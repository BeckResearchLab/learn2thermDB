"""Extract all protein sequences paired to the taxa that produced them."""
import logging
import os
import time

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

def extract_proteins_from_one(filepath: str):
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
    logger.info(f"Found {len(protein_sequences)} proteins")
    protein_sequences['seq_len'] = [len(s) for s in protein_sequences['sequence']]
    protein_sequences = pd.DataFrame(protein_sequences)
    return protein_sequences, seq_16srRNA

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

    # get the info in parallel
    outputs = Parallel(n_jobs=params['n_jobs'])(delayed(extract_proteins_from_one)(file) for file in taxa.values)

    # process into two files, one for proteins, one for 16s rRNA
    proteins, taxa_16srRNA = list(zip(*outputs))
    # proteins is a list of dataframe, taxa_16srRNA is a list of string or none if not found
    assert len(proteins) == len(taxa)
    assert len(taxa_16srRNA) == len(taxa)

    # assign the taxa index to every protein from that organism
    for i, protein_df in enumerate(proteins):
        protein_df['taxa_index'] = taxa.index[i]
    proteins = pd.concat(proteins, ignore_index=True)
    proteins.to_csv('./data/taxa/proteins.csv')

    # now save 16s
    taxa_16srRNA = pd.Series(index=taxa.index, data=taxa_16srRNA)
    taxa_16srRNA.to_csv('./data/taxa/16s_rRNA.csv')

    # save metrics
    metrics = {}
    metrics['n_total_sequences'] = int(len(proteins))
    metrics['n_taxa_with_16srRNA'] = int((~taxa_16srRNA.isna()).sum())
    with open('./data/metrics/s1.1_metrics.yaml', "w") as stream:
        yaml_dump(metrics, stream)