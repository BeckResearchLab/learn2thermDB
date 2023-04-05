"""Tools for File IO"""
from Bio import SeqIO
import gzip
import tempfile
import shutil
import os

import pandas as pd
import numpy as np

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
            records = SeqIO.parse(f_out, filetype)
            try:
                for r in records:
                    yield r
            except StopIteration:
                f_out.close()
                os.remove(tmp.name)
                raise


def csv_id_seq_iterator(csv_filepath: str, seq_col: str, index_col: str=None, id_filter: Collection = None, chunksize: int = 512, max_seq_length: int=None, **kwargs):
    """Returns a one by one iterator of seq ids and sequences to avoid OOM.
    
    Parameters
    ----------
    csv_filepath : str
        path to file containing data
    seq_col : str
        name of column containing sequences
    index_col: str, default None
        which column name is associated with the index, otherwise the 0th column will be used
    id_filter : Collection, Optional
        If given, only return sequences with the provided indexes
    chunksize : int, default 512
        Number of sequences that will be stored in memory at once.
    max_seq_length : int, default None
        Maximum length of sequence to return
    **kwargs passed to pandas read csv 
    """
    # first get the row numbers to consider
    if id_filter is not None:
        if index_col is not None:
            # get column positions to figure out which full col to load into memory
            columns = pd.read_csv(csv_filepath, nrows=1, **kwargs).columns
            index_col_position = np.argwhere(columns==index_col)[0][0]
            # load only that column
            row_indexes = pd.read_csv(csv_filepath, usecols=[index_col_position], **kwargs)
            row_indexes = pd.Series(row_indexes.set_index(index_col, drop=True).index)
        else:
            row_indexes = pd.read_csv(csv_filepath, usecols=[0]).index # just take the first column becuase we only need the indexes
            row_indexes = pd.Series(index=row_indexes, data=row_indexes)
        row_indexes_to_keep_mask = row_indexes.isin(id_filter)
        skiprows = lambda row_num: False if row_num==0 else not row_indexes_to_keep_mask.loc[row_num-1]
        logger.debug(f"{row_indexes_to_keep_mask.sum()} viable sequences in in file to iterate")
        seq_index_iterator = iter(list(row_indexes[row_indexes_to_keep_mask].values))
    else:
        skiprows=None

    for i, df_chunk in enumerate(pd.read_csv(csv_filepath, chunksize=chunksize, skiprows=skiprows, dtype=str, **kwargs)):
        if index_col is not None:
            df_chunk = df_chunk.set_index(index_col, drop=True)
        chunk = df_chunk[seq_col]
        logger.debug(f'Iterating chunk {i} seq in {csv_filepath}')
        for id_, seq in chunk.items():
            # skip long sequences
            if max_seq_length and len(seq) > max_seq_length:
                continue
            # in the case that there were no id filters, the id in the chunk corresponds to the correct sequence id
            # but if there was a filter, many rows were skipped and the indexes got jumbled, so we have to recapitulate
            # the correct seq index
            if id_filter is not None:
                yield next(seq_index_iterator), seq
            else:
                yield id_, seq
