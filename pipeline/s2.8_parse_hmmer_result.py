"""
TODO
Overall
-------
Inputs:
CSV chunks TODO: parquet
protien pair list

outputs:
output dictionary
"""
# system dependecies
import logging
import os
import sys
import tempfile
import time
from typing import Union

# library dependencies
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from codecarbon import OfflineEmissionsTracker
import duckdb as ddb
from joblib import Parallel, delayed
import pandas as pd
import pyhmmer
from yaml import dump as yaml_dump
from yaml import safe_load as yaml_load

# local dependencies
import learn2therm.database
import learn2therm.utils
import learn2therm.hmmer

## Paths
OUTPUT_DIR = './data/validation/hmmer/hmmer_labels/'

## get environmental variables
if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
    LOGLEVEL = getattr(logging, LOGLEVEL)
else:
    LOGLEVEL = logging.INFO
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'

# start logger/connect to log file
logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')

def create_accession_table():
    """
    Reads multiple CSV files and creates a DuckDB table from them.
    """
    # Creating temporary dir
    tmpdir = tempfile.mkdtemp(dir='./tmp', )
    # Establishing a connection with Duck DB
    conn = ddb.connect(tmpdir+'/pfam_mapping.db', read_only=False)
    # Create pairs table
    learn2therm.database.L2TDatabase._create_protein_pairs_table(conn, './data/')
    # Create the proteins_from_pairs table
    conn.execute("""
        CREATE TABLE pfam_mapping AS
        SELECT query_id AS pid, accession_id AS accession
        FROM read_parquet('./data/validation/hmmer/hmmer_outputs/*.parquet')
    """)
    # Committing DB
    unnaccounted = conn.execute("""SELECT COUNT(*) FROM (
        (SELECT meso_pid AS pid FROM pairs WHERE pid NOT IN (SELECT pid FROM pfam_mapping))
        UNION (SELECT thermo_pid AS pid FROM pairs WHERE pid NOT IN (SELECT pid FROM pfam_mapping))
    )""").fetchone()
    logger.info(f'Unaccounted for PIDs: {unnaccounted[0]}')
    conn.commit()
    conn.close()
    return tmpdir, tmpdir+'/pfam_mapping.db'


def process_pairs_table(dbpath, chunk_size:int, jaccard_threshold):
    """
    Processes the pairs table, calculates Jaccard similarity, and generates output CSV.
    """
    # Establish a connection with DuckDB
    conn = ddb.connect(dbpath, read_only=False)

    

    # Perform a join to get relevant information from the two tables
    query_string = """
        SELECT pairs.meso_pid, pairs.thermo_pid, mapping_meso.accession AS meso_accession, mapping_thermo.accession AS thermo_accession
        FROM pairs
        INNER JOIN pfam_mapping AS mapping_meso ON (pairs.meso_pid = mapping_meso.pid)
        INNER JOIN pfam_mapping AS mapping_thermo ON (pairs.thermo_pid = mapping_thermo.pid)
    """
    query_link = conn.execute(query_string)
        
    # Generate output parquet file
    # Execute the query
    data_remaining = True
    chunk_counter = 0  # Initialize the chunk counter
    while data_remaining:
        # Fetch the query result in chunks
        query_chunk = query_link.fetch_df_chunk(vectors_per_chunk=chunk_size)

        # Check if there is data remaining
        if len(query_chunk) == 0:
            data_remaining = False
            break

        # Calculate Jaccard similarity and determine functional status using apply function
        query_chunk[['score', 'bool_confirmed']] = query_chunk.apply(
            learn2therm.hmmer.evaluation_function,
            axis=1,
            args=(jaccard_threshold,),
            result_type='expand'
        )

        # Write DataFrame to parquet
        chunk_counter += 1  # Increment the chunk counter
        query_chunk = query_chunk[['meso_pid', 'thermo_pid', 'score']]
        query_chunk.to_parquet(f'{OUTPUT_DIR}{chunk_counter}_output.parquet')
        logger.info(f'Chunk {chunk_counter} of size {len(query_chunk)} written to parquet')

    # Commit the changes to the database
    conn.execute(f"CREATE TABLE labels AS SELECT * FROM read_parquet('{OUTPUT_DIR}/*.parquet')")
    logger.info("Saved labels to temporary database")

    conn.close()



if __name__ == '__main__':

    with open("./params.yaml", "r") as stream:
        params = yaml_load(stream)['run_hmmer']

    tracker = OfflineEmissionsTracker(
        project_name=f"s2.8",
        output_dir='./data/',
        country_iso_code='USA',
        region='Washington'
    ) 
    tracker.start()

    logger.info(f"Running script:{__file__}")

    # setup the database and get some pairs to run
    tmpdir_database, db_path = create_accession_table()

    # prepare output file
    try:
        os.makedirs(OUTPUT_DIR, exist_ok=True)
    except OSError as e:
        logger.error(f'Error creating directory: {e}')

    logger.info(f"Directory of output: {OUTPUT_DIR}, path to database {db_path}")

    logger.info('Creating pfam table') 
    threshold = params['jaccard_threshold']  # Set your own threshold
    chunksize = params['chunk_size'] # Vector size (2048 by default) * vector_multiple (we specify this).
    process_pairs_table(db_path, chunksize, threshold)
    logger.info('Pairs table processing completed')


    conn = ddb.connect(db_path, read_only=True)
    logger.info("Connected to DB to get statistics")

    # get some metrics
    metrics = {}
    count_pairs = conn.execute("SELECT COUNT(*) FROM pairs").fetchone()[0]
    count_labels = conn.execute("SELECT COUNT(*) FROM labels").fetchone()[0]
    assert count_pairs == count_labels, "Number of pairs and labels do not match"
    count_found = conn.execute("SELECT COUNT(*) FROM labels WHERE score IS NOT NULL").fetchone()[0]
    fraction_found = float(count_found/count_pairs)

    emissions = float(tracker.stop())

    std_jaccard = float(conn.execute("SELECT STDDEV(score) FROM labels").fetchone()[0])
    mean_jaccard = float(conn.execute("SELECT AVG(score) FROM labels").fetchone()[0])
    metrics = {
        'fraction_found': fraction_found,
        'std_jaccard': std_jaccard,
        'mean_jaccard': mean_jaccard,
        'co2': emissions
    }

    with open('./data/validation/hmmer/s2.8_metrics.yaml', 'w') as f:
        yaml_dump(metrics, f)

    
