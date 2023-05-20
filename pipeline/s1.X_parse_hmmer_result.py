"""
TODO
Overall
-------
Inputs:
CSV chunks
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
from codecarbon import OfflineEmissionsTracker
import duckdb as ddb
from joblib import Parallel, delayed
import pandas as pd
import pyhmmer

# local dependencies
import learn2therm.database
import learn2therm.utils

## Paths
OUTPUT_DIR = './data/protein_pairs/protein_pair_HMM_label'

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
    TODO
    """
    # Creating temporary dir
    tmpdir = tempfile.mkdtemp(dir='./tmp', )
    # Establishing a connection with Duck DB
    conn = ddb.connect(tmpdir+'/proteins_from_pairs.db', read_only=False)
    learn2therm.database.L2TDatabase._create_protein_pairs_table(conn, './data/')
    conn.execute("CREATE TABLE proteins_from_pairs AS SELECT query_id AS pid, accession_id AS accession_id FROM read_parquet('./data/protein_pairs/protein_pair_targets/uniprot_chunk_0.parquet')")
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
    if union == 0:
        return 0.0
    else:
        return intersection / union
    
def calculate_similarity():
    """
    TODO
    """
    # Create a dictionary to store the Jaccard similarity scores
    scores = defaultdict(float)
    
    # Calculate the Jaccard similarity score between each paired protein 
    
    # Create a dictionary to store the functional tuple values
    functional = {}
    
     # Set the functional tuple value based on the Jaccard similarity score threshold
    for (query1, query2), score in scores.items():
        if score >= threshold:
            functional[(query1, query2)] = ('Yes', score)
        else:
            functional[(query1, query2)] = ('No', score)
    
    return functional


def worker_function():
    """
    A wrapping function that is able to link the chunked pairs table with the protein_from_pairs table
    To create a dictionary that has meso_pid, thermo_pid, Functional?, and Jaccard score, and 
    write that dictionary to file as CSV.

    The Jaccard score is calculated from the meso_pid and thermo_pid which are paired by looking at their accessions.
    """
    pass



if __name__ == '__main__':
    # start logger/connect to log file
    logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')

    logger.info("TEST LOG")

    # setup the database and get some pairs to run
    tmpdir_database, db_path = create_accession_table()

    # prepare output file
    if not os.path.exists(OUTPUT_DIR):
        os.mkdir(OUTPUT_DIR)

    logger.info(f"Directory of output: {OUTPUT_DIR}, path to database {db_path}")
