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
    pass


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
    # Read the CSV files and create dictionaries with query IDs and accession IDs
    dict1 = parse_function_csv(file1)
    dict2 = parse_function_csv(file2)
    
    # Create a dictionary to store the Jaccard similarity scores
    scores = defaultdict(float)
    
    # Calculate the Jaccard similarity score between each protein in file1 and file2
    for query1, accs1 in dict1.items():
        for query2, accs2 in dict2.items():
            if query1 == query2:
                score = find_jaccard_similarity(set(accs1), set(accs2))
                scores[(query1, query2)] = score
    
    # Create a dictionary to store the functional tuple values
    functional = {}
    
     # Set the functional tuple value based on the Jaccard similarity score threshold
    for (query1, query2), score in scores.items():
        if score >= threshold:
            functional[(query1, query2)] = ('Yes', score)
        else:
            functional[(query1, query2)] = ('No', score)
    
    return functional


def write_function_output():
    """
    Writes a dictionary of protein query IDs and functional tuple values to a CSV file.

    Parameters
    ----------
    output_dict : Dict[str, Tuple[str, float]]
        A dictionary of protein query IDs and functional tuple values
    output_file : str
        File path to write the output CSV file
    """
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['File1', 'File2', 'Functional?', 'Jaccard Score']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        for query, (functional, score) in output_dict.items():
            writer.writerow({
                'File1': query[0],
                'File2': query[1],
                'Functional?': functional,
                'Jaccard Score': score
            })


if __name__ == '__main__':
    # start logger/connect to log file
    logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')
