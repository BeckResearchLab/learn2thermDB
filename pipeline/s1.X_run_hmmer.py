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
HMM_PATH = './data/HMM/Pfam-A.hmm'  # ./Pfam-A.hmm
PRESS_PATH = './data/HMM'
OUTPUT_DIR = './data/protein_pairs/protein_pair_targets'
WORKER_WAKE_UP_TIME = 25 # this is to ensure that if a worker that is about to be shut down due to previous task completetion doesn't actually start running

## get environmental variables
if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
    LOGLEVEL = getattr(logging, LOGLEVEL)
else:
    LOGLEVEL = logging.INFO
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'

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


def hmmpress_hmms(hmms_path, pfam_data_folder):
    """
    Presses the HMMs in the given HMM database and stores the resulting files in a specified directory.

    Parameters
    ----------
    hmmdb_path : str
        Path to the HMM database.

    pfam_data_folder : str, optional
        Path to the directory where the HMMs should be stored.

    Returns
    -------
    None

    Notes
    -----
    This function uses HMMER's hmmpress program to compress the HMMs in the given HMM database and
    stores the resulting files in the specified directory for faster access during future HMMER runs.
    If the specified directory does not exist, it will be created.
    """
    hmms = pyhmmer.plan7.HMMFile(hmms_path)
    pyhmmer.hmmer.hmmpress(hmms, pfam_data_folder)


def prefetch_targets(hmms_path: str):
    """
    Prefetch HMM profiles from a given HMM database.

    Parameters
    ----------
    hmms_path : str
        Path to the pressed HMM database.

    Returns
    -------
    targets : pyhmmer.plan7.OptimizedProfileBlock
        The HMM profiles loaded from the database.
    """
    # amino acid alphabet and prefetched inputs
    amino_acids = pyhmmer.easel.Alphabet.amino()
    optimized_profiles = list(pyhmmer.plan7.HMMPressedFile(hmms_path))
    targets = pyhmmer.plan7.OptimizedProfileBlock(
        amino_acids, optimized_profiles)
    return targets


def save_to_digital_sequences(dataframe: pd.DataFrame):
    """
    Save protein sequences from a DataFrame to a digital sequence block.

    Parameters
    ----------
    dataframe : pd.DataFrame
        DataFrame containing PIDs (Protein IDs) and sequences.

    Returns
    -------
    DigitalSequenceBlock
        A digital sequence block containing the converted sequences.
    """
    # Create empty list
    seqlist = []

    # Establish pyhmmer alphabet
    amino_acids = pyhmmer.easel.Alphabet.amino()

    # Convert proteins in dataframe to suitable format
    for _, row in dataframe.iterrows():
        pid = bytes(row['pid'], encoding='utf-8')
        seq_str = row['protein_seq']
        sequences = pyhmmer.easel.TextSequence(name=pid, sequence= seq_str)
        sequences = sequences.digitize(amino_acids)
        seqlist.append(sequences)
    
    # Convert so SequenceBlocks
    seqblock = pyhmmer.easel.DigitalSequenceBlock(amino_acids, seqlist)

    return seqblock


def run_pyhmmer(
        seqs: pyhmmer.easel.DigitalSequenceBlock,
        hmms_path: str,
        prefetch: bool = False,
        output_file: str = None,
        cpu: int = 4,
        eval_con: float = 1e-10):
    """
    Run HMMER's hmmscan program on a set of input sequences using with HMMs from a database.
    Parameters
    ----------
    seqs : pyhmmer.easel.DigitalSequenceBlock
        Path to the input sequence file.

    hmms_path : str
        Path to the HMM database.

    prefetch : bool, optional
        Specifies how the HMM are stored in meomry.

    output_file : str, optional
        Path to the output file if the users wants to write the file.

    cpu : int, optional
        The number of CPUs to use. Default is 4.

    eval_con : float, optional
        E-value threshold for domain reporting. Default is 1e-10.

    Returns
    -------
    all_hits : pyhmmer.plan7.TopHits or domtblout file
        If the output_file has a name, it will be written to a domtblout file.
        Otherwise, the user will get a list of pyhmmeer TopHits objects.

    Notes
    -----
    This function runs HMMER's hmmscan program on a set of input sequences
    using HMMs from a given database.
    The function supports two modes: normal mode and prefetching mode.
    In normal mode, the HMMs are pressed and stored in a directory before execution.
    In prefetching mode, the HMMs are kept in memory for faster search.
    """
    # ensure output_file has .domtblout extension
    if output_file is not None and not output_file.endswith('.domtblout'):
        output_file = f"{os.path.splitext(output_file)[0]}.domtblout"


    # HMM profile modes
    if prefetch:
        targets = prefetch_targets(hmms_path)
    else:
        targets = pyhmmer.plan7.HMMFile(PRESS_PATH)

    # HMMscan execution with or without saving output to file
    all_hits = list(pyhmmer.hmmer.hmmscan(seqs, targets, cpus=cpu, E=eval_con))
    # check if we should save the output
    if output_file is not None:
        with open(output_file, "wb") as dst:
            for i, hits in enumerate(all_hits):
                hits.write(dst, format="domains", header=i == 0)
    return all_hits


def parse_pyhmmer(all_hits, chunk_query_ids):
    """
    Parses the TopHit pyhmmer object getting the query and accession IDs and saves to a DataFrame

    Parameters
    ----------
    all_hits : list
        A list of TopHit objects from pyhmmer.
    chunk_query_ids : list
        A list of query IDs from the chunk.

    Returns
    -------
    pandas.DataFrame
        A dataframe containing the query and accession IDs.
    """
    # initialize an empty dictionary to store the data
    parsed_hits = {}

    # iterate over each protein hit
    for top_hits in all_hits:
        for hit in top_hits:
            # extract the query and accession IDs and decode the query ID
            query_id = hit.hits.query_name.decode('utf-8')
            accession_id = hit.accession.decode('utf-8')

            # if the query_id already exists in the dictionary, append the accession_id
            # to the existing value
            if query_id in parsed_hits:
                parsed_hits[query_id].append(accession_id)
            # otherwise, create a new key-value pair in the dictionary
            else:
                parsed_hits[query_id] = [accession_id]

    # find the query IDs that are missing from the parsed hits
    missing_query_ids = set(chunk_query_ids) - set(parsed_hits.keys())

    # add the missing query IDs with a placeholder value to indicate no accession information
    for missing_query_id in missing_query_ids:
        parsed_hits[missing_query_id] = ['No Accession Information']

    # create the DataFrame from the dictionary and convert list of accession IDs to string
    df = pd.DataFrame(parsed_hits.items(), columns=["query_id", "accession_id"])
    df["accession_id"] = df["accession_id"].apply(lambda x: ';'.join(x))

    return df


def worker_function(chunk_index, dbpath, chunked_pid_inputs, wakeup=None):
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
    conn = ddb.connect(dbpath, read_only=True)
    
    # get the unique pids from the chunked_pid_inputs
    pids = set(chunked_pid_inputs["pid"])

    # Only extract protein_seqs from the list of PID inputs
    placeholders = ', '.join(['?'] * len(pids))
    query = f"SELECT pid, protein_seq FROM proteins WHERE pid IN ({placeholders})"
    query_db = conn.execute(query, list(pids)).fetchall()

    # close db connection
    conn.close()

    # convert the query db to a dataframe
    result_df = pd.DataFrame(query_db, columns=['pid', 'protein_seq'])

    # convert string sequences to pyhmmer digital blocks
    sequences = save_to_digital_sequences(result_df)

    # run HMMER via pyhmmer
    hits = run_pyhmmer(
        seqs=sequences,
        hmms_path=PRESS_PATH,
        prefetch=True,
        cpu=1,
        eval_con=1e-5)
    
    # get the query IDs from the chunked_pid_inputs
    chunk_query_ids = chunked_pid_inputs["pid"].tolist()

    # Parse pyhmmer output and save to CSV file
    accessions_parsed = parse_pyhmmer(all_hits=hits, chunk_query_ids=chunk_query_ids)
    accessions_parsed.to_parquet(
        f'{OUTPUT_DIR}/{chunk_index}_output.parquet',
        index=False)

if __name__== "__main__":
    # start logger/connect to log file
    logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')
    
    logger.info("TEST LOG")

    # create pfam HMM directory (this was before HMM download script)
    try:
        os.makedirs('./data/HMM', exist_ok=True)
    except OSError as e:
        logger.error(f'Error creating directory: {e}')

    # press the HMM db
    hmmpress_hmms(HMM_PATH, PRESS_PATH)

    logger.info(f'Pressed HMM DB: {PRESS_PATH}')

    # Set up parallel processing and parsing
    chunk_size = int(sys.argv[1]) # Number of sequences to process in each chunk
    njobs = int(sys.argv[2])  # Number of parallel processes to use

    logger.info('Parallel processing parameters obtained')

    # setup the database and get some pairs to run
    tmpdir_database, db_path = load_protein_data()

    # prepare output file
    try:
        os.makedirs(OUTPUT_DIR, exist_ok=True)
    except OSError as e:
        logger.error(f'Error creating directory: {e}')

    logger.info(f"Directory of output: {OUTPUT_DIR}, path to database {db_path}")

    conn = ddb.connect(db_path, read_only=True)
    logger.info("Create PID dataframe from learn2therm database")
    proteins_in_pair_pids = conn.execute("SELECT pid FROM proteins WHERE proteins.pid IN (SELECT DISTINCT(pairs.meso_pid) FROM pairs) OR proteins.pid IN (SELECT DISTINCT(pairs.thermo_pid) FROM pairs)").df()
    logger.debug(f"Total number of protein in pairs: {len(proteins_in_pair_pids)} in pipeline")

    # chunking the PID so the worker function queries
    protein_pair_pid_chunks = [proteins_in_pair_pids[i:i + chunk_size]
                   for i in range(0, len(proteins_in_pair_pids), chunk_size)]
    
    logger.debug(f"number of protein chunks: {len(protein_pair_pid_chunks)}")
    
    # parallel computing on how many CPUs (n_jobs=)
    logger.info('Running pyhmmer in parallel on all chunks')

    Parallel(
        n_jobs=njobs)(
        delayed(worker_function)(
            chunk_index,
            db_path,
            protein_pair_pid_chunks,
            None) for chunk_index,
        protein_pair_pid_chunks in enumerate(protein_pair_pid_chunks))
    
    logger.info('Parallelization complete')