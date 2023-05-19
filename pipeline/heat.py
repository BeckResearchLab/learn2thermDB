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
from collections import Counter
import duckdb as ddb
from joblib import Parallel, delayed
import pandas as pd
import pyhmmer

# local dependencies
import learn2therm.database
import learn2therm.utils


def load_protein_data():
    """
    Load protein data into a temporary database.
    Returns
    -------
    tuple
        A tuple containing the path to the temporary directory and the path to the temporary database file.
    Notes
    -----
    This function performs the following steps:
    1. Creates a temporary directory.
    2. Establishes a connection with Duck DB using the temporary database file.
    3. Creates SQL tables for proteins and protein_pairs based on Parquet files.
    4. Commits the changes to the database.
    5. Closes the database connection.
    The temporary directory and database are returned as a tuple.
    """
    # Creating temporary dir
    tmpdir = tempfile.mkdtemp(dir='./tmp', )
    # Establishing a connection with Duck DB
    conn = ddb.connect(tmpdir+'/proteins.db', read_only=False)
    # Making a SQL table of proteins and protein_pairs
    learn2therm.database.L2TDatabase._create_protein_pairs_table(conn, './data/')
    conn.execute("CREATE TABLE proteins AS SELECT pid AS pid, protein_seq AS protein_seq FROM read_parquet('./data/proteins/uniprot_chunk_0.parquet')") #purely for testing
    # learn2therm.database.L2TDatabase._create_proteins_table(conn, './data/')
    # Committing DB
    conn.commit()
    conn.close()

    return tmpdir, tmpdir+'/proteins.db'


if __name__== "__main__":
    # logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    fh = logging.FileHandler('./logs/heat.log', mode='w')
    formatter = logging.Formatter(
        '%(filename)-12s %(asctime)s;%(funcName)-12s: %(levelname)-8s %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    logger.info("TEST LOG")

    # setup the database and get some pairs to run
    tmpdir_database, db_path = load_protein_data()

    logger.info(f"Directory of output: TODO, path to database {db_path}")

    conn = ddb.connect(db_path, read_only=True)
    logger.info("Connected to DB")

    protein_pair_pids = conn.execute("SELECT meso_pid, thermo_pid FROM pairs ORDER BY RANDOM()").df() # think
    logger.info("Create PID dataframe from learn2therm database")
    logger.debug(f"Total number of protein pairs: {len(protein_pair_pids)} in pipeline")

    # Check if all pids in proteins table exist in pairs table
    query = """
        SELECT COUNT(*) 
        FROM proteins
        WHERE pid NOT IN (
        SELECT DISTINCT(meso_pid) FROM pairs
        UNION
        SELECT DISTINCT(thermo_pid) FROM pairs
        )
        """
    result = conn.execute(query).fetchone()
    missing_count = result[0]

    logger.info(f"Number of missing pids: {missing_count}")

    

    # get the unique pids from the chunked_pid_inputs
    pids = set(protein_pair_pids["meso_pid"]).union(protein_pair_pids["thermo_pid"])

    logger.info(f"how many total protein pairs: {len(pids)}")
    logger.info(f"how many unique total protein pairs: {len(Counter(pids).keys())}")

    # print(protein_pair_pids_df, "\n\n\n", pids)

