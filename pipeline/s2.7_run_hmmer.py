"""
{HA note: this needs updating}
TODO:
1. Create a PFAM downloading script maybe press it in that script (data/HMM) [x]
2. Figure out the best way to deal with paths in light of above [x]
3. Update run_pyhmmer function after step 1 & 2
    a. Figure out output file in light of output dir [x]
4. Switch everything to parque [x]
5. Re-run after protein PIDs are all available in pairs? [~]

Old things:
should contain:
* parameters
* metrics
* worker state thing
"""
# system dependecies
import logging
import os
import shutil
import sys
import tempfile
import time
from typing import Union
from yaml import dump as yaml_dump
from yaml import safe_load as yaml_load

# library dependencies
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from codecarbon import OfflineEmissionsTracker
import duckdb as ddb
from joblib import Parallel, delayed
import pandas as pd
import numpy as np

# local dependencies
import learn2therm.database
import learn2therm.utils
import learn2therm.hmmer
import pyhmmer.plan7

## Paths
HMM_PATH = './data/validation/hmmer/Pfam-A.hmm'  # ./Pfam-A.hmm
PRESS_PATH = './tmp/'
OUTPUT_DIR = './data/validation/hmmer/hmmer_outputs/'
DB_PATH = './data/database.ddb'
## get environmental variables
if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
    LOGLEVEL = getattr(logging, LOGLEVEL)
else:
    LOGLEVEL = logging.INFO
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'

def worker_function(chunk_index, chunked_inputs, e_value: float=1e-6, prefetch=True, cpu=1, wakeup=None, scan=True, **kwargs):
    """
    A wrapping function that runs and parses pyhmmer in chunks.

    Parameters
    ----------
    chunk_index : int
        Number of sequence chunks.

    dbpath : str
        Path to the database.

    chunked_pid_inputs : pandas.DataFrame
        Dataframe containing chunked PID and sequences inputs.

    wakeup : int or None, optional
        Delay in seconds before starting the execution, by default None.
    Returns
    -------
    None

    Notes
    -----
    This function performs the following steps:
    1. Converts string sequences to pyhmmer digital blocks.
    2. Runs HMMER via pyhmmer with the provided sequences.
    3. Parses the pyhmmer output and saves it to a CSV file.

    The parsed pyhmmer output is saved in the directory specified by OUTPUT_DIR,
    with each chunk having its own separate output file named '{chunk_index}_output.csv'.

    If the wakeup parameter is specified, the function will wait for the specified
    number of seconds before starting the execution.
    """
    # we want to wait for execution to see if this worker is actually being used
    # or if it is in the process of being killed
    if wakeup is not None:
        time.sleep(wakeup)

    # convert string sequences to pyhmmer digital blocks
    sequences = learn2therm.hmmer.save_to_digital_sequences(chunked_inputs)

    # run HMMER via pyhmmer
    if prefetch:
        hits = learn2therm.hmmer.run_pyhmmer(
            seqs=sequences,
            pressed_path=PRESS_PATH,
            prefetch=prefetch,
            cpu=cpu,
            eval_con=e_value,
            scan=scan,
            **kwargs
        )
    else:
        hits = learn2therm.hmmer.run_pyhmmer(
            seqs=sequences,
            hmms_path=HMM_PATH,
            prefetch=False,
            cpu=cpu,
            eval_con=e_value,
            scan=scan,
            **kwargs
        )

    # Parse pyhmmer output and save to CSV file
    accessions_parsed = learn2therm.hmmer.parse_pyhmmer(all_hits=hits, chunk_query_ids=chunked_inputs['pid'].tolist(), scanned=scan)
    accessions_parsed.to_parquet(
        f'{OUTPUT_DIR}/{chunk_index}_output.parquet',
        index=False)

if __name__== "__main__":

    # load dvc params
    # load params
    with open("./params.yaml", "r") as stream:
        params = yaml_load(stream)['run_hmmer']
    # start logger/connect to log file
    logger = learn2therm.utils.start_logger_if_necessary('root', LOGFILE, LOGLEVEL, filemode='w')
    logger.info(f"Loaded params: {params}")
    tracker = OfflineEmissionsTracker(
        project_name=f"s1.1",
        output_dir='./data/',
        country_iso_code='USA',
        region='Washington'
    ) 
    tracker.start()
    logger.info(f"Running script:{__file__}")

    # create pfam HMM directory (this was before HMM download script)
    shutil.rmtree(OUTPUT_DIR, ignore_errors=True)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Set up parallel processing and parsing
    chunk_size = params['chunk_size']
    njobs = params['njobs']
    logger.info(f'Parallel processing parameters obtained: chunk_size: {chunk_size}, njobs: {njobs}')

    # prepare output file
    try:
        os.makedirs(OUTPUT_DIR, exist_ok=True)
    except OSError as e:
        logger.error(f'Error creating directory: {e}')
    logger.info(f"Directory of output: {OUTPUT_DIR}, path to database {DB_PATH}")

    # start a clock
    start_time = time.time()

    conn = ddb.connect(DB_PATH, read_only=True)
    logger.info("Connected to learn2therm database")
    if params['dev_sample_data']:
        proteins_q = conn.execute(f"SELECT pid, protein_seq FROM proteins WHERE proteins.pid IN (SELECT DISTINCT(pairs.meso_pid) FROM pairs) OR proteins.pid IN (SELECT DISTINCT(pairs.thermo_pid) FROM pairs) LIMIT {params['dev_sample_data']}")
    else:
        proteins_q = conn.execute("SELECT pid, protein_seq FROM proteins WHERE proteins.pid IN (SELECT DISTINCT(pairs.meso_pid) FROM pairs) OR proteins.pid IN (SELECT DISTINCT(pairs.thermo_pid) FROM pairs)")
    
    # get number of hmms for evalue calc
    profiles = list(pyhmmer.plan7.HMMFile(HMM_PATH))
    n_hmms = len(profiles)
    del profiles
    logger.info(f"Number of HMMs: {n_hmms}")

    # here we need to prefetch or not manually and only need to do it once
    if params['prefetch']:
        learn2therm.hmmer.hmmpress_hmms(HMM_PATH, PRESS_PATH)
        logger.info(f'Pressed HMM DB: {PRESS_PATH}')
        targets = learn2therm.hmmer.prefetch_targets(PRESS_PATH)
        logger.info(f"prefetched targets, will use same block for all chunks")
        wrapper = lambda chunk_index, pid_chunk: worker_function(
            chunk_index, pid_chunk, cpu=njobs, prefetch=targets, e_value=params['e_value'], scan=params['scan'], Z=n_hmms)
    else:
        wrapper = lambda chunk_index, pid_chunk: worker_function(
            chunk_index, pid_chunk, cpu=njobs, prefetch=False, e_value=params['e_value'], scan=params['scan'], Z=n_hmms)
    
    complete = False
    chunk_index = 0
    total_processed = 0
    while not complete:
        pid_chunk = proteins_q.fetch_df_chunk(vectors_per_chunk=params['chunk_size'])
        logger.info(f"Loaded chunk of size {len(pid_chunk)}")
        if len(pid_chunk) == 0:
            complete = True
            break
        wrapper(chunk_index, pid_chunk)
        logger.info(f"Ran chunk, validating results")

        # do a check on the output
        df = pd.read_parquet(f'{OUTPUT_DIR}/{chunk_index}_output.parquet')
        assert set(list(pid_chunk['pid'].values)) == set(list(df['query_id'].values)), "Not all query ids are in the output file"

        logger.info(f"Completed chunk {chunk_index} with size {len(pid_chunk)}")
        total_processed += len(pid_chunk)
        chunk_index += 1
        
    # compute some metrics
    co2 = float(tracker.stop())
    logger.info(f"Total CO2 emissions: {co2} kg")

    # get the total number of proteins we labeled with accessions
    con = ddb.connect(DB_PATH, read_only=True)
    con.execute(f"CREATE TEMP TABLE results AS SELECT * FROM read_parquet('{OUTPUT_DIR}/*.parquet')")
    total_proteins = con.execute("SELECT COUNT(*) FROM proteins WHERE proteins.pid IN (SELECT DISTINCT(pairs.meso_pid) FROM pairs) OR proteins.pid IN (SELECT DISTINCT(pairs.thermo_pid) FROM pairs)").fetchone()[0]
    total_proteins_check = con.execute("SELECT COUNT(*) FROM ((SELECT DISTINCT(meso_pid) AS pid FROM pairs) UNION (SELECT DISTINCT(thermo_pid) AS pid FROM pairs))").fetchone()[0]
    completed_proteins = con.execute("SELECT COUNT(DISTINCT(query_id)) FROM results").fetchone()[0]
    labeled_proteins = con.execute("SELECT COUNT(DISTINCT(query_id)) FROM results WHERE accession_id != ''").fetchone()[0]

    logger.info(f"Processed {total_processed} proteins of the expected {total_proteins} (check: {total_proteins_check}) proteins in pairs in the database, {completed_proteins} of which were are in the mapping.")

    metrics = {
        "n_proteins_in_pairs": int(total_proteins),
        "n_proteins_labeled": int(labeled_proteins),
        "n_proteins_in_mapping": int(completed_proteins),
        "co2_emissions": co2,
        'execution_time': float((time.time() - start_time)/60)
    }

    with open(f"./data/validation/hmmer/s2.7_metrics.yaml", "w") as f:
        yaml_dump(metrics, f)
