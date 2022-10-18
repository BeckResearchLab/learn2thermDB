"""Tools for File IO"""
from Bio import SeqIO
import gzip
import tempfile
import shutil
import os

import pandas as pd

from typing import Collection

import logging
logger = logging.getLogger(__name__)

def seq_io_gnuzipped(filepath: str, filetype: str):
    """Open gnuzipped sequence files.
    
    Wrapper around biopython SeqIO.

    Parameters
    ----------
    filepath : str
        path to zipped data file
    filetype : str
        type of sequence data file within the zipped file. Passed directly to SeqIO
    
    Returns
    -------
    list of records
    """
    with gzip.open(filepath, 'rb') as f_in:
        tmp = tempfile.NamedTemporaryFile(delete=False)
        tmp.close()
        with open(tmp.name, 'w+b') as f_out:
            shutil.copyfileobj(f_in, f_out)
            f_out.close()
        with open(tmp.name, 'rt') as f_out:
            seq = SeqIO.parse(f_out, filetype)
            records = [r for r in seq]
            f_out.close()
        os.remove(tmp.name)
    return records

def csv_id_seq_iterator(csv_filepath: str, seq_col: str, id_filter: Collection = None, chunksize: int = 512, **kwargs):
    """Returns a one by one iterator of seq ids and sequences to avoid OOM.
    
    Parameters
    ----------
    csv_filepath : str
        path to file containing data
    seq_col : str
        name of column containing sequences
    id_filter : Collection, Optional
        If given, only return sequences with the provided indexes
    chunksize : int, default 512
        Number of sequences that will be stored in memory at once.
    **kwargs passed to pandas read csv 
    """
    for i, df_chunk in enumerate(pd.read_csv(csv_filepath, chunksize=chunksize, **kwargs)):
        chunk = df_chunk[seq_col]
        logger.debug(f'Iterating chunk {i} seq in {csv_filepath}')
        # filter down indexes
        if id_filter is not None:
            mask = chunk.index.isin(id_filter)
            chunk = chunk[chunk.index[mask]]
            logger.debug(f"{len(chunk)} viable sequences in chunk")
        for id_, seq in chunk.items():
            yield id_, seq
