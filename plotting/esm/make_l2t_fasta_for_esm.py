"""Select a sample of l2t proteins that occur in pairs equal to the fraction of Atlas used,
and compute and save the ESM embeddings.
"""
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import torch
import numpy as np
import duckdb as ddb
import pandas as pd
from esm import Alphabet, FastaBatchedDataset, pretrained, MSATransformer
import pathlib
from typing import List
import tempfile
import os

import logging
logging.basicConfig(filename='l2t_esm_fasta.log',
                    filemode='w',
                    format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                    datefmt='%H:%M:%S',
                    level=logging.INFO)

logger = logging.getLogger('l2t_esm')

FRAC = 0.06010428815245909 # fraction of original dataset to use


if __name__ == '__main__':

    # get data from duck db
    con = ddb.connect('../../data/database.ddb', read_only=True)
    count = con.execute("SELECT COUNT(pid) FROM proteins WHERE proteins.pid IN (SELECT DISTINCT(pairs.meso_pid) FROM pairs) OR proteins.pid IN (SELECT DISTINCT(pairs.thermo_pid) FROM pairs)").fetchone()[0]
    n_to_keep = int(FRAC * count)

    dataframe = con.execute(f"SELECT pid, protein_seq FROM proteins WHERE proteins.pid IN (SELECT DISTINCT(pairs.meso_pid) FROM pairs) OR proteins.pid IN (SELECT DISTINCT(pairs.thermo_pid) FROM pairs) ORDER BY RANDOM() LIMIT {n_to_keep}").df()
    logger.info(f'Got {len(dataframe)} proteins')

    # write fasta file
    fasta_file_loc = './data/esm_fasta_input.fasta'
    if os.path.exists(fasta_file_loc):
        fasta_file = fasta_file_loc
    else:
        records = []
        ids_ = []
        for i, row in dataframe.iterrows():
            if row['pid'] in ids_:
                continue
            else:
                record = SeqRecord(Seq(row['protein_seq']), id=row['pid'], description='')
                records.append(record)
                ids_.append(row['pid'])
        fasta_file = open(fasta_file_loc, 'w')
        SeqIO.write(records, fasta_file, 'fasta')
        fasta_file.close()
        fasta_file = fasta_file.name
        del records
    fp = pathlib.Path(fasta_file)
    logger.info(f"Using fasta file {fasta_file}")

    
