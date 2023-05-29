"""
This script aims to test out the evaluation_function logic as well as the Jaccard similarity calculation logic.
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

### from twoground.py
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


def preprocess_accessions(meso_acession, thermo_accession):
    """
    TODO
    """
    # Convert accessions to sets
    meso_accession_set = set(meso_acession.split(';'))
    thermo_accession_set = set(thermo_accession.split(';'))
    return meso_accession_set, thermo_accession_set


def calculate_jaccard_similarity(meso_accession_set, thermo_accession_set):
    """
    Calculates Jaccard similarity between meso_pid and thermo_pid pairs via their accessions.
    """
    # Calculate Jaccard similarity
    intersection = len(meso_accession_set.intersection(thermo_accession_set))
    union = len(meso_accession_set.union(thermo_accession_set))
    jaccard_similarity = intersection / union if union > 0 else 0

    return jaccard_similarity


# Define the evaluation function for the apply function
def evaluation_function(row, jaccard_threshold):
    """TODO
    """
    # Get the accessions
    meso_acc = row['meso_accession']
    thermo_acc = row['thermo_accession']

    # parsing accessions logic
    # Handle NULL values
    if pd.isnull(meso_acc) and pd.isnull(thermo_acc):
        score = None
        functional = None
    elif pd.notnull(meso_acc) and pd.notnull(thermo_acc):
        # Preprocess the accessions
        meso_acc_set, thermo_acc_set = preprocess_accessions(meso_acc, thermo_acc)

        if meso_acc_set == {"No accession Information"} and thermo_acc_set == {"No accession Information"}:
            score = None
            functional = None
        else:
            score = calculate_jaccard_similarity(meso_acc_set, thermo_acc_set)
            functional = score > jaccard_threshold
    else:
        # Handle unmatched rows
        score = None
        functional = False
    
    return {'functional?': functional, 'score': score}

if __name__ == '__main__':
    # Initialize logger
    logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')
    logger.info(f"Running {__file__}")

    # Define a sample row of data
    
    # only slighty matching case
    # Expected output: {'functional?': False, 'score': 0.3333333333333333}
    # sample_row = {
    #     'meso_accession': 'ABC;DEF',
    #     'thermo_accession': 'DEF;GHI'
    #     }

    # matching case
    # Expected output: {'functional?': True, 'score': 1.0}
    # sample_row = {
    # 'meso_accession': 'ABC;DEF;GHI',
    # 'thermo_accession': 'DEF;GHI;JKL'
    # }

    # case with one missing accession
    # Expected output: {'functional?': None, 'score': None}
    # sample_row = {
    # 'meso_accession': 'ABC;DEF',
    # 'thermo_accession': None
    # }

    # case with empty accessions (this one ended up being a bug)
    # Expected output: {'functional?': None, 'score': None}
    sample_row = {
    'meso_accession': '',
    'thermo_accession': ''
    }


    # Define the Jaccard threshold for functional label
    jaccard_threshold = 0.5

    # Call the evaluation_function with the sample row and threshold
    result = evaluation_function(sample_row, jaccard_threshold)

    # Print the output to check if it matches your expectations
    logger.info(result)
