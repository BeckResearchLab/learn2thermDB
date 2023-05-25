"""
This script is just for testing stuff

TODO:
You could use df.apply() as a way to create the two columns. It is also vectorized; much quicker.
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
    meso_accession_set = set(meso_accession.split(';'))
    thermo_accession_set = set(thermo_accession.split(';'))

    # print accession sets
    logger.debug("Meso Accession set:", meso_accession_set)
    logger.debug("Thermo Accession set:", thermo_accession_set)

    # Calculate Jaccard similarity
    intersection = len(meso_accession_set.intersection(thermo_accession_set))
    union = len(meso_accession_set.union(thermo_accession_set))
    jaccard_similarity = intersection / union if union > 0 else 0

    return jaccard_similarity



def process_pairs_table(dbpath, jaccard_threshold):
    """
    Processes the pairs table, calculates Jaccard similarity, and generates output CSV.
    """
    # Establish a connection with DuckDB
    conn = ddb.connect(dbpath, read_only=False)

    # Perform inner join to get relevant information
    query1 = """
        CREATE OR REPLACE TABLE joined_pairs AS 
        SELECT p.meso_pid, p.thermo_pid, pr.accession AS meso_accession, pr2.accession AS thermo_accession
        FROM pairs AS p
        RIGHT JOIN proteins_from_pairs AS pr ON (p.meso_pid = pr.pid)
        RIGHT JOIN proteins_from_pairs AS pr2 ON (p.thermo_pid = pr2.pid)
    """
    conn.execute(query1)

    # Fetch all rows from the query result
    number_of_joined_pairs = conn.execute("SELECT COUNT(*) FROM joined_pairs").df()
    logger.debug(f'number of rows inside the joined_pairs table:{number_of_joined_pairs}')
    query_result = conn.execute("SELECT * FROM joined_pairs")
    rows = query_result.fetchall()

    # Generate output CSV file
    try:
        with open(f'{OUTPUT_DIR}output.csv', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['meso_pid', 'thermo_pid', 'functional?', 'score'])

            # Iterate over joined pairs
            for row in rows:
                meso_pid = row[0]
                thermo_pid = row[1]
                meso_accession = row[2]
                thermo_accession = row[3]

                # Calculate Jaccard similarity
                jaccard_similarity = calculate_jaccard_similarity(row)

                # Determine if it meets the Jaccard threshold
                meets_threshold = jaccard_similarity >= jaccard_threshold

                # Write row to CSV
                writer.writerow([meso_pid, thermo_pid, 'yes' if meets_threshold else 'no', jaccard_similarity])
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
    process_pairs_table(db_path, threshold)
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