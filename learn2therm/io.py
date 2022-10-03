"""Tools for File IO"""
from Bio import SeqIO
import gzip
import tempfile
import shutil
import os

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

