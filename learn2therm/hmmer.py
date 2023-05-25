"""Functions for HMMER
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
        press_path: str,
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

    press_path : str
        Path to the pressed HMM database

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
        targets = pyhmmer.plan7.HMMFile(press_path)

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