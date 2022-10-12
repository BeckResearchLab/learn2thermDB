"""Wraps up use of the NCBI blast software for use with this package.

We are not using a known database, we are comparing lists of sequences, so we need tools
to support this.

Note BioPython provides a python API to the command line
"""
import tempfile
import os
import re

import numpy as np
import pandas as pd

import logging
logger = logging.getLogger(__name__)

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
        logger.info("Creating temporary files to deposit blast inputs and outputs.")
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
        logger.info("Removing temporary files used by blast")
        os.remove(self.qt)
        os.remove(self.st)
        os.remove(self.ot)
        

class BlastMetrics:
    """Handles computation of metrics for each alignment in a blast record.
    
    Parameters
    ----------
    blast_record : the record containing all hits for a query

    Notes
    -----
    metric names:
        local_X - metric uses only information from a single HSP
        scaled_local_X - metric uses single HSP information normalized by global sequence length. 
            This is an attempt to filter out small but accurate chunks without needing the global alignment/qcoverage.
            Eg if match is |, mismatch is X and gap is -
                |||||||||||||||X|||||||||-------------------------------------------------------------------------------------------------
            Could score very well using local alignment score even though there is a huge gap tail.
        global_X - metric uses global alignment, currently cannot do without HSP tiling
    """
    def __init__(self, blast_record):
        self.record = blast_record
        self.qid = self.record.query.split(' ')[0]
        logger.info(f"Query {self.qid} with {len(self.record.alignments)} alignments with ids {[a.hit_id for a in self.record.alignments]}.")

    def compute_metric(self, metric_name: str):
        """Compute the metric with specified name for each alignment"""
        if not hasattr(self, metric_name):
            raise ValueError(f"No metric found with name : {metric_name}")
        else:
            metric = getattr(self, metric_name)
        
        logger.info(f"Computing metric `{metric_name}` for all alignments in query {self.qid}")

        outputs = []
        for alignment in self.record.alignments:
            outputs.append((self.qid, alignment.hit_id, metric(alignment)))
        return pd.DataFrame(data=outputs, columns=['query_id', 'subject_id', metric_name])

    @staticmethod
    def raw_gap_excluding_percent_id(n_matches, n_gaps, n_columns):
        """Percent matches in sequence, excluding gaps.
        
        Parameters
        ----------
        n_matches : int, number of matches in match columns
        n_gaps : number of gaps in match columns
        n_columns : total number of alignment match columns
        """
        return n_matches / (n_columns - n_gaps)

    @staticmethod
    def raw_gap_including_percent_id(n_matches, n_columns):
        """Percent matches in sequence, including gaps gaps.
        
        Parameters
        ----------
        n_matches : int, number of matches in match columns
        n_columns : total number of alignment match columns
        """
        return n_matches / (n_columns)

    @staticmethod
    def raw_gap_compressed_percent_id(n_matches, n_gaps, n_columns, n_compressed_gaps):
        """Percent matches in sequence, including but compressing gaps.
        
        Parameters
        ----------
        n_matches : int, number of matches in match columns
        n_gaps : number of gaps in match columns
        n_columns : total number of alignment match columns
        n_compressed_gaps : number of compressed gaps in match columns
        """
        return n_matches / (n_columns - n_gaps + n_compressed_gaps)

    def local_gap_compressed_percent_id(self, alignment):
        """Percent matches in match sequence, including but compressing gaps.
        
        The largest local HSP score is used
        """
        scores = []
        for hsp in alignment.hsps:
            n_matches = hsp.identities
            n_gaps = hsp.gaps
            n_columns = len(hsp.query)
            n_compressed_gaps = len(re.findall('-+', hsp.query))+len(re.findall('-+', hsp.sbjct))
            scores.append(self.raw_gap_compressed_percent_id(n_matches, n_gaps, n_columns, n_compressed_gaps))
        return max(scores)

    def scaled_local_query_percent_id(self, alignment):
        """Percent matches in query sequence based on best HSP."""
        ids = [hsp.identities for hsp in alignment.hsps]
        return max(ids)/self.record.query_length

    def scaled_local_symmetric_percent_id(self, alignment):
        """Percent matches compared to average seq length of query and subject based on best HSP"""
        ids = [hsp.identities for hsp in alignment.hsps]
        return 2*max(ids)/(self.record.query_length + alignment.length)

    def local_E_value(self, alignment):
        """E value of HSP with most identities."""
        ids = []
        Es = []
        for hsp in alignment.hsps:
            ids.append(hsp.identities)
            Es.append(hsp.expect)
        return Es[np.argmax(ids)]
