"""
This script is just for testing stuff

TODO:
You could use df.apply() as a way to create the two columns. It is also vectorized; much quicker.
Just connect to tmp table -> debug strat

TODO:
Modularize the jacard similarity functoin. Inputs should only be a list of accessions
-> make function of splitting of accessions
-> make function used for the apply for the functional? score (True, False, Nan) 
-> make function
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
import csv
from collections import defaultdict
import duckdb as ddb
import glob
from joblib import Parallel, delayed
import pandas as pd
import pyhmmer

# local dependencies
import learn2therm.database
import learn2therm.utils


## Paths
OUTPUT_DIR = './tmp/results/'

## get environmental variables
if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
    LOGLEVEL = getattr(logging, LOGLEVEL)
else:
    LOGLEVEL = logging.INFO
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'


def create_accession_table():
    """
    Reads multiple CSV files and creates a DuckDB table from them.
    """
    # Creating temporary dir
    tmpdir = tempfile.mkdtemp(dir='./tmp', )
    # Establishing a connection with Duck DB
    conn = ddb.connect(tmpdir+'/proteins_from_pairs.db', read_only=False)
    # Create pairs table
    learn2therm.database.L2TDatabase._create_protein_pairs_table(conn, './data/')
    # Create the proteins_from_pairs table
    conn.execute("""
        CREATE TABLE proteins_from_pairs AS
        SELECT query_id AS pid, accession_id AS accession
        FROM read_csv_auto('./data/protein_pairs/protein_pair_targets/*.csv', HEADER=TRUE)
    """)
    # Committing DB
    conn.commit()
    conn.close()
    return tmpdir, tmpdir+'/proteins_from_pairs.db'



def calculate_jaccard_similarity(pair):
    """
    Calculates Jaccard similarity between meso_pid and thermo_pid pairs via their accessions.
    """
    meso_pid, thermo_pid, meso_accession, thermo_accession = pair

    # Print input pairs for debugging
    logger.debug("Meso PID:", meso_pid)
    logger.debug("Thermo PID:", thermo_pid)
    logger.debug("Meso Accession:", meso_accession)
    logger.debug("Thermo Accession:", thermo_accession)

    # Check if both pairs have "No accession Information"
    if meso_accession == "No accession Information" and thermo_accession == "No accession Information":
        return 0

    # Convert accessions to sets
    if isinstance(meso_accession, str):
        meso_accession_set = set(meso_accession.split(';'))
    else:
        meso_accession_set = set()

    if isinstance(thermo_accession, str):
        thermo_accession_set = set(thermo_accession.split(';'))
    else:
        thermo_accession_set = set()

    # print accession sets
    logger.debug("Meso Accession set:", meso_accession_set)
    logger.debug("Thermo Accession set:", thermo_accession_set)

    # Calculate Jaccard similarity
    intersection = len(meso_accession_set.intersection(thermo_accession_set))
    union = len(meso_accession_set.union(thermo_accession_set))
    jaccard_similarity = intersection / union if union > 0 else 0

    return jaccard_similarity




def process_pairs_table(dbpath, chunk_size:int, jaccard_threshold):
    """
    Processes the pairs table, calculates Jaccard similarity, and generates output CSV.
    """
    # Establish a connection with DuckDB
    conn = ddb.connect(dbpath, read_only=False)

    # Perform a join to get relevant information from the two tables
    query1 = """
        CREATE OR REPLACE TABLE joined_pairs AS 
        SELECT p.meso_pid, p.thermo_pid, pr.accession AS meso_accession, pr2.accession AS thermo_accession
        FROM pairs AS p
        RIGHT JOIN proteins_from_pairs AS pr ON (p.meso_pid = pr.pid)
        RIGHT JOIN proteins_from_pairs AS pr2 ON (p.thermo_pid = pr2.pid)
    """
    conn.execute(query1)

    # Generate output CSV file
    try:
        data_remaining = True
        chunk_counter = 0  # Initialize the chunk counter
        while data_remaining:
            # Fetch the query result in chunks
            query_result = conn.execute("SELECT * FROM joined_pairs").fetch_df_chunk(chunk_size)

            # Check if there is data remaining
            if query_result.empty:
                data_remaining = False
                break

            # Calculate Jaccard similarity and determine functional status using apply function
            query_result['score'] = query_result.apply(calculate_jaccard_similarity, axis=1)
            query_result['functional?'] = query_result['score'].apply(lambda x: 'yes' if x >= jaccard_threshold else 'no')

            # Write DataFrame to CSV
            chunk_counter += 1  # Increment the chunk counter
            query_result.to_csv(f'{OUTPUT_DIR}{chunk_counter}_output.csv', index=False, columns=['meso_pid', 'thermo_pid', 'functional?', 'score'])

    except IOError as e:
        logger.warning(f"Error writing to CSV file: {e}")

    conn.close()




if __name__ == '__main__':
    # start logger/connect to log file
    logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')

    logger.info("TEST LOG")

    # setup the database and get some pairs to run
    tmpdir_database, db_path = create_accession_table()

    # prepare output file
    try:
        os.makedirs(OUTPUT_DIR, exist_ok=True)
    except OSError as e:
        logger.error(f'Error creating directory: {e}')

    logger.info(f"Directory of output: {OUTPUT_DIR}, path to database {db_path}")



    logger.info('Creating process pair table') 
    threshold = 0.5  # Set your own threshold
    chunksize = 5 # Vector size (2048 by default) * vector_multiple (we specify).
    process_pairs_table(db_path, chunksize ,threshold)
    logger.info('Script done') 

    conn = ddb.connect(db_path, read_only=False)
    logger.info("Connected to DB")

    # Check if all pids in proteins table exist in pairs table
    query2 = """
        SELECT COUNT(*) 
        FROM proteins_from_pairs
        WHERE pid NOT IN (
        SELECT DISTINCT(meso_pid) FROM pairs
        UNION
        SELECT DISTINCT(thermo_pid) FROM pairs
        )
        """
    result = conn.execute(query2).fetchone()
    missing_count = result[0]

    logger.debug(f"Number of missing pids: {missing_count}")