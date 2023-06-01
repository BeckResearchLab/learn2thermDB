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

# local dependencies
import learn2therm.database
import learn2therm.utils
import learn2therm.hmmer

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

def worker_function(chunk_index, chunked_pid_inputs, e_value: float=1e-6, prefetch=True, wakeup=None):
    """
    A wrapping function that runs and parses pyhmmer in chunks.

    Parameters
    ----------
    chunk_index : int
        Number of sequence chunks.

    dbpath : str
        Path to the database.

    chunked_pid_inputs : pandas.DataFrame
        Dataframe containing chunked PID inputs.

    wakeup : int or None, optional
        Delay in seconds before starting the execution, by default None.
    Returns
    -------
    None

    Notes
    -----
    This function performs the following steps:
    1. Queries the database to get sequences only from chunked_pid_inputs.
    2. Converts the query result to a dataframe.
    3. Converts string sequences to pyhmmer digital blocks.
    4. Runs HMMER via pyhmmer with the provided sequences.
    5. Parses the pyhmmer output and saves it to a CSV file.

    The parsed pyhmmer output is saved in the directory specified by OUTPUT_DIR,
    with each chunk having its own separate output file named '{chunk_index}_output.csv'.

    If the wakeup parameter is specified, the function will wait for the specified
    number of seconds before starting the execution.
    """
    # we want to wait for execution to see if this worker is actually being used
    # or if it is in the process of being killed
    if wakeup is not None:
        time.sleep(wakeup)
    
    # define paths for input and output files
    # output_file_path = f'./tmp/results/{chunk_index}_output'

    # query the database to get sequences only from chunked_pid_inputs
    conn = ddb.connect(DB_PATH, read_only=True)
    
    # get the unique pids from the chunked_pid_inputs
    pids = set(chunked_pid_inputs)

    # Only extract protein_seqs from the list of PID inputs
    placeholders = ','.join(['?'] * len(pids))
    query = f"SELECT pid, protein_seq FROM proteins WHERE pid IN ({placeholders})"
    query_db = conn.execute(query, list(pids)).fetchall()

    # close db connection
    conn.close()

    # convert the query db to a dataframe
    result_df = pd.DataFrame(query_db, columns=['pid', 'protein_seq'])

    # convert string sequences to pyhmmer digital blocks
    sequences = learn2therm.hmmer.save_to_digital_sequences(result_df)

    # run HMMER via pyhmmer
    if prefetch:
        hits = learn2therm.hmmer.run_pyhmmer(
            seqs=sequences,
            pressed_path=PRESS_PATH,
            prefetch=True,
            cpu=1,
            eval_con=e_value)
    else:
        hits = learn2therm.hmmer.run_pyhmmer(
            seqs=sequences,
            hmm_path=HMM_PATH,
            prefetch=False,
            cpu=1,
            eval_con=e_value)

    # Parse pyhmmer output and save to CSV file
    accessions_parsed = learn2therm.hmmer.parse_pyhmmer(all_hits=hits, chunk_query_ids=pids)
    accessions_parsed.to_parquet(
        f'{OUTPUT_DIR}/{chunk_index}_output.parquet',
        index=False)

if __name__== "__main__":

    # load dvc params
    # load params
    with open("./params.yaml", "r") as stream:
        params = yaml_load(stream)['run_hmmer']
    # start logger/connect to log file
    logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')
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
    logger.info("Create PID dataframe from learn2therm database")
    if params['dev_sample_data']:
        proteins_in_pair_pids = conn.execute(f"SELECT DISTINCT(pid) FROM proteins WHERE proteins.pid IN (SELECT DISTINCT(pairs.meso_pid) FROM pairs) OR proteins.pid IN (SELECT DISTINCT(pairs.thermo_pid) FROM pairs) LIMIT {params['dev_sample_data']}").df()
        logger.info(f"Sampling data: {len(proteins_in_pair_pids)} total proteins")
    else:
        proteins_in_pair_pids = conn.execute("SELECT DISTINCT(pid) FROM proteins WHERE proteins.pid IN (SELECT DISTINCT(pairs.meso_pid) FROM pairs) OR proteins.pid IN (SELECT DISTINCT(pairs.thermo_pid) FROM pairs)").df()
        logger.info(f"Total number of protein in pairs: {len(proteins_in_pair_pids)} in pipeline")
    conn.close()
    proteins_in_pair_pids = list(proteins_in_pair_pids.values.reshape(-1))

    # chunking the PID so the worker function queries
    protein_pair_pid_chunks = [proteins_in_pair_pids[i:i + chunk_size]
                   for i in range(0, len(proteins_in_pair_pids), chunk_size)]
    
    # press the HMM db
    if params['prefetch']:
        learn2therm.hmmer.hmmpress_hmms(HMM_PATH, PRESS_PATH)
        logger.info(f'Pressed HMM DB: {PRESS_PATH}')

        wrapper = lambda chink_index, pid_chunks: worker_function(chink_index, pid_chunks, prefetch=True, e_value=params['e_value'])
    else:
        wrapper = lambda chink_index, pid_chunks: worker_function(chink_index, pid_chunks, prefetch=False, e_value=params['e_value'])

    # parallel computing on how many CPUs (n_jobs=)
    logger.info(f'Running pyhmmer in parallel on {len(protein_pair_pid_chunks)} chunks')

    Parallel(
        n_jobs=njobs)(
        delayed(wrapper)(
            chunk_index,
            protein_pair_pid_chunks
        ) for chunk_index, protein_pair_pid_chunks in enumerate(protein_pair_pid_chunks))

    # compute some metrics
    co2 = float(tracker.stop())
    logger.info(f"Total CO2 emissions: {co2} kg")

    # get the total number of proteins we labeled with accessions
    con = ddb.connect()
    con.execute(f"CREATE TEMP TABLE results AS SELECT * FROM read_parquet('{OUTPUT_DIR}/*.parquet')")
    total_proteins = con.execute("SELECT COUNT(DISTINCT(query_id)) FROM results").fetchone()[0]
    labeled_proteins = con.execute("SELECT COUNT(DISTINCT(query_id)) FROM results WHERE accession_id IS NOT NULL").fetchone()[0]

    metrics = {
        "n_proteins_in_pairs": int(total_proteins),
        "n_proteins_labeled": int(labeled_proteins),
        "co2_emissions": co2,
        'execution_time': float((time.time() - start_time)/60)
    }

    with open(f"./data/validation/hmmer/s2.7_metrics.yaml:", "w") as f:
        json.dump(metrics, f)
