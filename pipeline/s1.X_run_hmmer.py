"""
{HA note: this needs updating}
TODO:
1. Create a PFAM downloading script maybe press it in that script
2. Figure out the best way to deal with paths in light of above
3. Update run_pyhmmer function after step 1 & 2

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
import learn2therm.utils

## Paths
HMM_PATH = './data/Pfam-A.hmm'  # ./Pfam-A.hmm
PROTEIN_SEQ_DIR = './data/taxa/proteins/'
OUTPUT_DIR = './data/taxa_pairs/protein_pair_targets'
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
    TODO:
    check s1.3 for inspiration
    """
    # stolen from Evan's t1.0 scripts
    with tempfile.TemporaryDirectory(dir='./tmp') as tmpdir:
      # Establishing a connection with Duck DB
      conn = ddb.connect(tmpdir+'proteins.db', read_only=False)
      # Making a SQL table
      conn.execute("CREATE TABLE proteins AS SELECT taxid AS taxid, pid AS pid, protein_seq AS sequence FROM read_parquet('./data/proteins/uniprot_chunk_0.parquet')") # changes
      # Committing DB
      conn.commit()

      # Create table of taxa pairs
      conn.execute("CREATE TEMP TABLE pair_labels AS SELECT * FROM read_parquet('./data/taxa_pairs/pair_labels/*.parquet')") # change
      conn.execute("CREATE TEMP TABLE pair_scores AS SELECT * FROM read_parquet('./data/taxa_pairs/alignment/*.parquet')") # change
      conn.execute("CREATE TABLE pairs AS SELECT * FROM pair_labels INNER JOIN pair_scores ON (pair_labels.__index_level_0__ = pair_scores.__index_level_0__) WHERE pair_labels.is_pair=True")
      conn.commit()
      conn.execute("CREATE INDEX meso_index ON pairs (subject_id)")
      conn.execute("CREATE INDEX thermo_index ON pairs (query_id)")
      conn.commit()
      conn.close()

    return tmpdir, tmpdir+'/proteins.db'

def prefetch_targets(hmms_path: str):
    """
    Prefetch HMM profiles from a given HMM database.
    Parameters
    ----------
    hmms_path : str
        Path to the HMM database.
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
    TODO

    Parameters
    ----------
    dataframe : pd.DataFrame
        _description_

    Returns
    -------
    _type_
        _description_
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
        seqs: str,
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
    if not output_file.endswith('.domtblout'):
        output_file = f"{os.path.splitext(output_file)[0]}.domtblout"

    # HMM profile modes
    if prefetch:
        targets = prefetch_targets(hmms_path)
    else:
        targets = pyhmmer.plan7.HMMFile("../data/HMM")

    # HMMscan execution with or without saving output to file
    all_hits = list(pyhmmer.hmmer.hmmscan(seqs, targets, cpus=cpu, E=eval_con))
    # check if we should save the output
    if output_file is not None:
        with open(output_file, "wb") as dst:
            for i, hits in enumerate(all_hits):
                hits.write(dst, format="domains", header=i == 0)
    return all_hits


def parse_pyhmmer(all_hits):
    """
    Parses the TopHit pyhmmer object getting the query and accession IDs and saves to a DataFrame
    Parameters
    ----------
    all_hits : list
        A list of TopHit objects from pyhmmer.
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

    # create the DataFrame from the dictionary and convert list of accession IDs to string
    df = pd.DataFrame(parsed_hits.items(), columns=["query_id", "accession_id"])
    df["accession_id"] = df["accession_id"].apply(lambda x: ';'.join(x))

    return df


def worker_function(chunk_index, dbpath, pid_inputs, wakeup=None):
    """
    A wrapping function that runs and parses pyhmmer in chunks
    Parameters
    ----------
    chunk_index : int
        number of sequences chunks
    sequences : str
        a list of dataframe containing protein sequences
    """
    # we want to wait for execution to see if this worker is actually being used
    # or if it is in the process of being killed
    if wakeup is not None:
        time.sleep(wakeup)
    
    # define paths for input and output files
    input_file_path = f'../tmp/results/{chunk_index}_input'
    output_file_path = f'../tmp/results/{chunk_index}_output'

    # query ID
    conn = ddb.connect(dbpath, read_only=True)
    query_db = conn.execute("SELECT pid, protein_seq FROM proteins WHERE pid = {pid_inputs}").df()



    # convert string sequences to pyhmmer digital blocks
    save_to_digital_sequences(query_db)

    # run HMMER via pyhmmer
    hits = run_pyhmmer(
        input_file=input_file_path,
        hmms_path=HMM_PATH,
        prefetch=True,
        output_file=output_file_path,
        cpu=1,
        eval_con=1e-5)

    # Parse pyhmmer output and save to CSV file
    accessions_parsed = parse_pyhmmer(all_hits=hits)
    accessions_parsed.to_csv(
        f'../tmp/results/{chunk_index}_output.csv',
        index=False)

if __name__== "__main__":
    # start logger/connect to log file
    logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')
    
    logger.info("TEST LOG")

    # create pfam HMM directory
    if not os.path.exists('./data/HMM'):
        os.mkdir('./data/HMM')

    # # Set up parallel processing and parsing
    # total_size = int(sys.argv[1])  # Number of total sequences read
    # # Number of sequences to process in each chunk
    # chunk_size = int(sys.argv[2])
    # njobs = int(sys.argv[3])  # Number of parallel processes to use

    logger.info('Parallel processing parameters obtained')

    # setup the database and get some pairs to run
    tmpdir_database, db_path = load_protein_data()

    # prepare output file
    if not os.path.exists(OUTPUT_DIR):
        os.mkdir(OUTPUT_DIR)

    logger.info(f"Directory of output: {OUTPUT_DIR}, path to database {db_path}")

    conn = ddb.connect(db_path, read_only=True)
    pairs = conn.execute("SELECT query_id, subject_id FROM pairs ORDER BY RANDOM()").fetchall()
    logger.info(f"Total number of taxa pairs: {len(pairs)} in pipeline")

    # # processing query proteins
    # query_test = load_protein_data()
    # logger.info('Sample loaded')
    
    # query_db = query_test[["pid", "protein_seq"]]
    # protein_list = query_test.set_index("pid").iloc[:total_size]
    # protein_list.index.name = None
    # logger.info('Sample preprocessed')


    # chunking the data to chunk_size sequence bits (change if sample or all
    # proteins)
    test_chunks = [protein_list[i:i + chunk_size]
                   for i in range(0, len(protein_list), chunk_size)]
    
    # parallel computing on how many CPUs (n_jobs=)
    logger.info('Running pyhmmer in parallel on all chunks')

    Parallel(
        n_jobs=njobs)(
        delayed(worker_function)(
            chunk_index,
            sequences,
            "test") for chunk_index,
        sequences in enumerate(test_chunks))
    
    logger.info('Parallelization complete')