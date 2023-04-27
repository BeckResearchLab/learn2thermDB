"""
{HA note: this needs updating}
Script that labels putative BLAST protein pairs using domains from HMMER.

HMMER is a tool for searching sequence databases for homologs of protein.
HMM databases include Pfam, BFD.

This script HMMER to search all proteins against the HMM database. It then
labels putative BLAST protein pairs as homologous if they share some
amount of similarity on identified domains.

Eg. if a protein is labeled with 4 domains, and three are shared with
the other protein which has 5 domains, then the two proteins have 
Jaccard similarity of 3/9 = 0.33. Some threshold is set for this
to be validated as a protein pair.

Inputs to script:
-----------------
1 For DB version 1.0 using refseq and bacdive: proteins in ./data/taxa/proteins
  in the form of csv files of prien sequences and associated taxa ids.
1 For DB version 1.1 using uniprot: proteins in ./data/proteins in the form
  parquet files
2 parquet files at ./data/taxa_pairs/protein_alignment of the form "taxa_pair_XX-YY.parquet"
  Where XX and YY are the taxa ids of the two taxa in the pair. Contents are protein IDs and
  BLAST metics.

Outputs of script:
------------------
1 parquet files at ./data/taxa_pairs/hmmer_val of the form "taxa_pair_XX-YY.parquet"
  that are 1:1 mirrors of the blast files, with domains IDs identified for each protein,
  and jaccard score

Expected steps, approximately:
------------------------------
1. Load the protein sequences and associated taxonomic information from the input CSV files
  for DB version 1.0, or from the input Parquet files for DB version 1.1.
2. Use HMMER to search all proteins against the HMM database, which includes Pfam and BFD.
3. Identify the domains for each protein using HMMER output, and calculate the Jaccard similarity
  between the domains of each pair of proteins.
    In accomplishing this step, you may need to save the HMMER domains to file, and then load and 
    compute the Jaccard similarity in a separate step.
    Otherwise, do both at runtime 
4. Save the results, including protein IDs in parquet files at 
  ./data/taxa_pairs/hmmer_val of the form "taxa_pair_XX-YY.parquet", which are 1:1
  mirrors of the blast files. Each parquet file should contain the domains IDs identified for each protein
  and the calculated Jaccard score.

Make sure to handle errors and exceptions appropriately and document any assumptions or limitations of the script.
"""
# system dependecies
import logging
import os
import sys
from typing import Union


# library dependencies
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from joblib import Parallel, delayed
import pyhmmer


# local dependencies
import learn2therm.utils

# get environmental variables
if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
    LOGLEVEL = getattr(logging, LOGLEVEL)
else:
    LOGLEVEL = logging.INFO
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'

if __name__== "__main__":
    # start logger/connect to log file
    logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')
