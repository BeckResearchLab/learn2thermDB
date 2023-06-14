"""Tools for running hmmer and saving the output."""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from codecarbon import OfflineEmissionsTracker
import duckdb as ddb
from joblib import Parallel, delayed
import pandas as pd
import pyhmmer

from typing import Union
import logging
logger = logging.getLogger(__name__)

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

def prefetch_targets(pressed_path: str):
    """
    Prefetch HMM profiles from a given HMM database.

    Parameters
    ----------
    pressed_path : str
        Path to the pressed HMM database.

    Returns
    -------
    targets : pyhmmer.plan7.OptimizedProfileBlock
        The HMM profiles loaded from the database.
    """
    # amino acid alphabet and prefetched inputs
    amino_acids = pyhmmer.easel.Alphabet.amino()
    optimized_profiles = list(pyhmmer.plan7.HMMPressedFile(pressed_path))
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
        columns: 'pid', 'protein_seq'

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
        seqs: Union[pyhmmer.easel.DigitalSequenceBlock, str],
        hmms_path: str = None,
        pressed_path: str = None,
        prefetch: Union[bool, pyhmmer.plan7.OptimizedProfileBlock] = False,
        output_file: str = None,
        cpu: int = 4,
        scan: bool=True,
        eval_con: float = 1e-10,
        **kwargs
    ):
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
        Also, can be a pyhmmer.plan7.OptimizedProfileBlock object.


    output_file : str, optional
        Path to the output file if the users wants to write the file.

    cpu : int, optional
        The number of CPUs to use. Default is 4.
    
    scan: bool, optional
        Whether to run hmmscan or hmmsearch. Default is True (hmmscan).

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
    if not hmms_path and not pressed_path:
        raise ValueError("Must specifity one of hmm path (to a .hmm file) or pressed_path (containing .h3m, etc.)")

    # ensure output_file has .domtblout extension
    if output_file is not None and not output_file.endswith('.domtblout'):
        output_file = f"{os.path.splitext(output_file)[0]}.domtblout"

    # HMM profile modes
    if prefetch:
        if isinstance(prefetch, pyhmmer.plan7.OptimizedProfileBlock):
            targets = prefetch
        elif pressed_path is None:
            raise ValueError("Spcified prefetch but did not pass a path to pressef iles")
        else:
            targets = prefetch_targets(pressed_path)
    else:
        if hmms_path is None:
            raise ValueError("Spcified prefetch but did not pass a path to the .hmm file")
        targets = pyhmmer.plan7.HMMFile(hmms_path)

    # are the sequences preloaded?
    if isinstance(seqs, str):
        seqs = pyhmmer.easel.SequenceFile(seqs, format='fasta', digital=True, alphabet=pyhmmer.easel.Alphabet.amino())
    else:
        pass

    # HMMscan execution with or without saving output to file
    if hasattr(seqs, '__len__'):
        seqs_size = len(seqs)
    else:
        seqs_size = "In file, unknown length"

    if hasattr(targets, '__len__'):
        targets_size = len(targets)
    else:
        targets_size = "In file, unknown length"

    if scan:
        logger.info(f"Running hmmscan... {seqs_size} sequences against {targets_size} HMMs, using {cpu} CPUs, additional kwargs: {kwargs}")
        all_hits = pyhmmer.hmmer.hmmscan(seqs, targets, cpus=cpu, incE=eval_con, **kwargs)
    else:
        logger.info(f"Running hmmsearch... {targets_size} HMMs against {seqs_size} seqs, using {cpu} CPUs, additional kwargs: {kwargs}")
        all_hits = pyhmmer.hmmer.hmmsearch(targets, seqs, cpus=cpu, incE=eval_con, **kwargs)
    # check if we should save the output
    if output_file is not None:
        with open(output_file, "wb") as dst:
            for i, hits in enumerate(all_hits):
                hits.write(dst, format="domains", header=i == 0)
    return all_hits

def parse_pyhmmer(all_hits, chunk_query_ids, scanned: bool = True):
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
            # check e value
            if hit.evalue > top_hits.incE:
                continue

            # extract the query and accession IDs and decode the query ID
            if scanned:
                query_id = hit.hits.query_name.decode('utf-8')
                accession_id = hit.accession.decode('utf-8')
            else:
                query_id = hit.name.decode('utf-8')
                accession_id = hit.hits.query_accession.decode('utf-8')

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
        parsed_hits[missing_query_id] = [""]

    # create the DataFrame from the dictionary 
    df = pd.DataFrame(parsed_hits.items(), columns=["query_id", "accession_id"])

    # convert list of accession IDs to string
    df["accession_id"] = df["accession_id"].apply(lambda x: ";".join(x) if x else "")

    return df

def preprocess_accessions(meso_acession, thermo_accession):
    """
    TODO
    """
    # Convert accessions to sets
    if meso_acession == '':
        meso_accession_set = set()
    else:
        meso_accession_set = set(meso_acession.split(';'))
    if thermo_accession == '':
        thermo_accession_set = set()
    else:
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
    if type(meso_acc) == str:
        meso_acc_set, thermo_acc_set = preprocess_accessions(meso_acc, thermo_acc)
    elif type(meso_acc) == list:
        meso_acc_set = set(meso_acc)
        thermo_acc_set = set(thermo_acc)
    else:
        raise ValueError("meso_acc must be either a string or a list")

    # parsing accessions logic
    if not meso_acc_set and not thermo_acc_set:
        score = None
        functional = None
    elif meso_acc_set and thermo_acc_set:
        # Preprocess the accessions
        score = calculate_jaccard_similarity(meso_acc_set, thermo_acc_set)
        functional = score > jaccard_threshold
    else:
        # Handle unmatched rows
        score = 0.0
        functional = False
    return score, functional