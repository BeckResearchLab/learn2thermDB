"""
This script is just for testing stuff
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



def find_jaccard_similarity(set1: set, set2: set) -> float:
    """
    Calculates the Jaccard similarity score between two sets.
    """
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return 0.0 if union == 0 else intersection / union



def calculate_similarity(pid_pairs, accessions1, accessions2, threshold):
    """
    Calculates the Jaccard similarity score for a pair of proteins.
    """
    print("pid_pairs:", pid_pairs)
    print("accessions1:", accessions1)
    print("accessions2:", accessions2)
    set1 = set(accessions1.split(';')) if accessions1 and accessions1 != 'No Accession Information' else set()
    set2 = set(accessions2.split(';')) if accessions2 and accessions2 != 'No Accession Information' else set()
    score = find_jaccard_similarity(set1, set2)
    
    # if isinstance(pid_pairs, str):
    #     # If pid_pairs is a string, convert it to a tuple
    #     pid_pairs = tuple(pid_pairs.split(','))
    
    return (pid_pairs[0], pid_pairs[1], score)




def worker_function(db_path, threshold):
    """
    Query the database in chunks and calculate Jaccard similarity using joblib.
    """
    conn = ddb.connect(db_path, read_only=True)
    cursor = conn.cursor()

    cursor.execute("""
        SELECT 
            pairs.meso_pid, 
            pairs.thermo_pid, 
            group_concat(proteins_from_pairs_meso.accession) as meso_accessions,
            group_concat(proteins_from_pairs_thermo.accession) as thermo_accessions
        FROM pairs
        LEFT JOIN proteins_from_pairs AS proteins_from_pairs_meso ON pairs.meso_pid = proteins_from_pairs_meso.pid
        LEFT JOIN proteins_from_pairs AS proteins_from_pairs_thermo ON pairs.thermo_pid = proteins_from_pairs_thermo.pid
        GROUP BY pairs.meso_pid, pairs.thermo_pid
    """)

    chunk_size = 10000  # You can adjust the chunk size
    njobs = 4  # You can adjust the number of jobs
    scores = {}

    while True:
        rows = cursor.fetchmany(chunk_size)
        if not rows:
            break
        results = Parallel(n_jobs=njobs)(delayed(calculate_similarity)(meso_pid, thermo_pid, meso_accessions, thermo_accessions) for meso_pid, thermo_pid, meso_accessions, thermo_accessions in rows)
        for result in results:
            (pid_pairs, thermo_pid, score) = result
            functional = 'Yes' if score >= threshold else 'No'
            scores[(pid_pairs[0], pid_pairs[1])] = (functional, score)

    return scores




def write_to_csv(scores, filename):
    """
    Writes the scores to a CSV file.
    """
    try:
        with open(filename, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['meso_pid', 'thermo_pid', 'functional?', 'Jaccard Score'])
            for (meso_pid, thermo_pid), (functional, score) in scores.items():
                writer.writerow([meso_pid, thermo_pid, functional, score])
    except IOError as e:
        print(f"Error writing to CSV file: {e}")


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

    threshold = 0.5  # Set your own threshold
    scores = worker_function(db_path, threshold)
    logger.info(f"Jaccard scores: {scores}")  

    # Write scores to a CSV file
    csv_filename = os.path.join(OUTPUT_DIR, 'scores.csv')
    write_to_csv(scores, csv_filename)
    logger.info(f"Scores written to {csv_filename}")