"""Get some fasta files to play around with from the DB"""
import os
import pandas as pd
import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

SOURCE = '../../data/taxa/proteins/'

n_taxa = len(os.listdir(SOURCE)) - 1

# proteins from 100 random organisms for training
def seq_iter(indexes):
    n = 0
    for i in indexes:
        df = pd.read_csv(SOURCE+f'taxa_index_{i}.csv', sep=';')
        seqs = df['protein_seq'].values
        for s in seqs:
            seq = Seq(s)
            seq = SeqRecord(
                seq,
                id=str(n),
            )
            n += 1
            yield seq

train_indexes = np.random.choice(n_taxa, size=100, replace=False)
train_seqs = seq_iter(train_indexes)
with open("./data/train.fasta", "w") as output_handle:
    SeqIO.write(train_seqs, output_handle, "fasta")
    
test_indexes = []
while len(test_indexes) < 3:
    index = np.random.choice(n_taxa)
    if index not in train_indexes:
        test_indexes.append(index)
test_seqs = seq_iter(test_indexes)
with open("./data/test.fasta", "w") as output_handle:
    SeqIO.write(test_seqs, output_handle, "fasta")