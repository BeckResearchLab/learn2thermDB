"""Extract all protein sequences paired to the taxa that produced them."""
import logging
import os
import time
import fcntl

from joblib import delayed, Parallel
import pandas as pd
from yaml import dump as yaml_dump
from yaml import safe_load as yaml_load

import learn2therm.io
import learn2therm.utils

# get environmental variables
if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
    LOGLEVEL = getattr(logging, LOGLEVEL)
else:
    LOGLEVEL = logging.INFO
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'

# where we will deposit_data
OUTPUT_FILE_PROTEINS = './data/taxa/proteins.csv'
OUTPUT_FILE_16srRNA = './data/taxa/16s_rRNA.csv'

def extract_proteins_from_one(taxa_index: int, filepath: str):
    """Open a zipped GBFF file and extract all protein sequences.
    
    Parameters
    ----------
    filepath : path to gbff.gz file

    Returns
    -------
    DataFrame of (sequences, descriptions, seq length), 16s rRNA sequence
    """
    # get the logger in subprocesses
    logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL)

    logger.debug(f"Extracting file {filepath}")
    refseq_records = learn2therm.io.seq_io_gnuzipped(filepath, 'genbank')

    # Each record has features, and some of those features are CDS. Get those translations
    protein_sequences = {
        'sequence': [],
        'desc': [],
    }
    seq_16srRNA = None
    for record in refseq_records:
        for feature in record.features:
            # check for 16s
            if feature.type == 'rRNA' and feature.qualifiers['product'][0] == '16S ribosomal RNA':
                seq_16srRNA = feature.extract(record.seq)
            # check for protein with a translation
            if feature.type == 'CDS' and 'translation' in feature.qualifiers:
                if len(feature.qualifiers['translation']) > 1:
                    raise ValueError(f'Multiple translations for feature')
                protein_sequences['sequence'].append(feature.qualifiers['translation'][0])
                if 'product' in feature.qualifiers:
                    protein_sequences['desc'].append(feature.qualifiers['product'][0])
                else:
                    protein_sequences['desc'].append(None)
    if seq_16srRNA is None:
        logger.info(f"Could not find 16s rRNA in {filepath}")
    logger.info(f"Found {len(protein_sequences['sequence'])} proteins")
    protein_sequences['seq_len'] = [len(s) for s in protein_sequences['sequence']]
    
    # save to file
    with open(OUTPUT_FILE_PROTEINS, "a") as g:
        fcntl.flock(g, fcntl.LOCK_EX)
        for i in range(len(protein_sequences['seq_len'])):
            g.write(f"{taxa_index};{protein_sequences['sequence'][i]};{protein_sequences['desc'][i].replace(';', ',')};{protein_sequences['seq_len'][i]}\n")
        fcntl.flock(g, fcntl.LOCK_UN)
    with open(OUTPUT_FILE_16srRNA, "a") as g:
        fcntl.flock(g, fcntl.LOCK_EX)
        g.write(f"{taxa_index},{seq_16srRNA}\n")
        fcntl.flock(g, fcntl.LOCK_UN)
    
    num_proteins = len(protein_sequences['sequence'])
    has_16srRNA = seq_16srRNA is not None
    return num_proteins, has_16srRNA

if __name__ == '__main__':
    # connect to log file
    logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')

    # DVC tracked parameters
    with open("./pipeline/s1_data_processing_params.yaml", "r") as stream:
        params = yaml_load(stream)['get_protein_sequences']
    logger.info(f"Loaded parameters: {params}")

    # load the filepaths and their identities
    taxa = pd.read_csv('./data/taxa/taxa_info_and_ogt.csv', index_col=0,  usecols=[0,3])['filepath']
    # select a sample if specified
    if params['n_sample']:
        taxa = taxa.iloc[:params['n_sample']]
        logger.info(f"Using only the first {params['n_sample']} files.")

    # start the files, we have to write intermittedly because of memory issues
    with open(OUTPUT_FILE_PROTEINS, "w") as g:
        g.write("taxa_index;protein_seq;protein_desc;protein_len\n")
    with open(OUTPUT_FILE_16srRNA, "w") as g:
        g.write("taxa_index,seq_16srRNA\n")

    # get the info in parallel
    outputs = Parallel(n_jobs=params['n_jobs'])(delayed(extract_proteins_from_one)(taxa_index, file) for taxa_index, file in taxa.items())

    # process into two files, one for proteins, one for 16s rRNA
    proteins_counts, has_16srRNA = list(zip(*outputs))
    # proteins is a list of dataframe, taxa_16srRNA is a list of string or none if not found
    assert len(proteins_counts) == len(taxa)
    assert len(has_16srRNA) == len(taxa)

    # save metrics
    metrics = {}
    metrics['n_total_sequences'] = int(sum(proteins_counts))
    metrics['n_taxa_with_16srRNA'] = int(sum(has_16srRNA))
    with open('./data/metrics/s1.1_metrics.yaml', "w") as stream:
        yaml_dump(metrics, stream)