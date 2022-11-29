"""Extract all protein sequences paired to the taxa that produced them."""
import logging
import os
import time
import fcntl
import shutil

from joblib import delayed, Parallel
from codecarbon import OfflineEmissionsTracker
import pandas as pd
import numpy as np
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
OUTPUT_DIR_PROTEINS = './data/taxa/proteins/'
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
    logger.info(f"Found {len(protein_sequences['sequence'])} proteins for taxa {taxa_index}")
    protein_sequences['seq_len'] = [len(s) for s in protein_sequences['sequence']]
    
    # save to file
    # first we will have to create protein file for this taxa
    with open(OUTPUT_DIR_PROTEINS+f'taxa_index_{taxa_index}.csv', "w") as g:
        g.write("seq_id;protein_seq;protein_desc;protein_len\n")
    logger.debug(f"Created new proteins file for taxa {taxa_index}")

    # save proteins to the protein file for this taxa and 16s to the global file fir 16s
    with open(OUTPUT_DIR_PROTEINS+f'taxa_index_{taxa_index}.csv', "a") as g:
        fcntl.flock(g, fcntl.LOCK_EX)
        for i in range(len(protein_sequences['seq_len'])):
            g.write(f"{taxa_index}.{i};{protein_sequences['sequence'][i]};{protein_sequences['desc'][i].replace(';', ',')};{protein_sequences['seq_len'][i]}\n")
        fcntl.flock(g, fcntl.LOCK_UN)
    with open(OUTPUT_FILE_16srRNA, "a") as g:
        fcntl.flock(g, fcntl.LOCK_EX)
        g.write(f"{taxa_index},{seq_16srRNA}\n")
        fcntl.flock(g, fcntl.LOCK_UN)
    
    num_proteins = len(protein_sequences['sequence'])
    has_16srRNA = seq_16srRNA is not None
    logger.debug(f"Completed protein deposit for taxa {taxa_index}")
    return num_proteins, has_16srRNA

if __name__ == '__main__':
    # connect to log file
    logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')

    # DVC tracked parameters
    with open("./params.yaml", "r") as stream:
        params = yaml_load(stream)['get_protein_sequences']
    logger.info(f"Loaded parameters: {params}")

    # load the filepaths and their identities
    taxa = pd.read_csv('./data/taxa/taxa_info_and_ogt.csv', index_col=0,  usecols=[0,3])['filepath']
    # select a sample if specified
    if params['n_sample']:
        taxa = taxa.iloc[:params['n_sample']]
        logger.info(f"Using only the first {params['n_sample']} files.")

    # start the files, we have to write intermittedly because of memory issues
    with open(OUTPUT_FILE_16srRNA, "w") as g:
        g.write("taxa_index,seq_16srRNA\n")
    # start a dir for proteins
    shutil.rmtree(OUTPUT_DIR_PROTEINS, ignore_errors=True)
    os.makedirs(OUTPUT_DIR_PROTEINS)

    # get the info in parallel
    with OfflineEmissionsTracker(
        project_name=f"s1.1",
        output_dir='./logs/',
        country_iso_code='USA',
        region='Washington'
    ) as tracker:
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

    # save length plot
    length_list = []
    for i, f in enumerate(os.listdir(OUTPUT_DIR_PROTEINS)):
        lengths = pd.read_csv(OUTPUT_DIR_PROTEINS+f, usecols=[3], sep=';').values
        length_list.append(lengths)
    import seaborn as sns
    import matplotlib.pyplot as plt
    sns.set_context('talk')
    fig, ax = plt.subplots()
    sns.histplot(np.vstack(length_list), ax=ax)
    plt.savefig(f'./data/plots/protein_length_hist.png', bbox_inches='tight', dpi=250)