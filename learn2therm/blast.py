"""Wraps up use of the NCBI blast software for use with this package.

We are not using a known database, we are comparing lists of sequences, so we need tools
to support this.

Note BioPython provides a python API to the command line
"""
import tempfile
import os


class BlastFiles:
    """Temporary files for use with BLAST CLI.
    
    Blast expects two input FASTA and produces an XML. The FASTA are redundant to CSV
    we already have. These files are created for the context and removed after completion.

    Parameters
    ----------
    query_iterator : iterator of (seq id, sequence)
        sequences to be used as query
    subject_iterator : iterator of (seq id, sequence)
        sequences to be used as the "database"

    Returns
    -------
    query_filename : str, name of fasta file with query sequences
    subject_filename : str, name of fasta file with subject sequences
    output_filename : str, name of output file for blast to save reuslts, will be deleted out of context
    """
    def __init__(self, query_iterator, subject_iterator):
        # we have to create the temporary fasta files
        query_temp = tempfile.NamedTemporaryFile('w', delete=False)
        self.qt = query_temp.name
        for id_, seq in query_iterator:
            if seq == 'None' or seq is None:
                continue
            query_temp.write(f">{id_}\n{seq}\n")
        query_temp.close()

        subject_temp = tempfile.NamedTemporaryFile('w', delete=False)
        self.st = subject_temp.name
        for id_, seq in subject_iterator:
            if seq == 'None' or seq is None:
                continue
            subject_temp.write(f">{id_}\n{seq}\n")
        subject_temp.close()

        # create the output xml file
        out_temp = tempfile.NamedTemporaryFile('w', delete=False)
        self.ot = out_temp.name

    def __enter__(self):
        return self.qt, self.st, self.ot

    def __exit__(self, type, value, traceback):
        os.remove(self.qt)
        os.remove(self.st)
        os.remove(self.ot)
        