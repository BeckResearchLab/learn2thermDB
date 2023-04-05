"""Assigns mesophilic, thermophilic labels.

A method for dealing with taxa with data but not OGT should be passed.
"""
import ast
from locale import normalize
import logging
import os

import matplotlib.pyplot as plt
import pandas as pd
from yaml import dump as yaml_dump
from yaml import safe_load as yaml_load

import learn2therm.utils

# get environmental variables
if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
    LOGLEVEL = getattr(logging, LOGLEVEL)
else:
    LOGLEVEL = logging.INFO
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'


if __name__ == '__main__':
    # connect to log file
    logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')

    # DVC tracked parameters
    with open("./params.yaml", "r") as stream:
        params = yaml_load(stream)['label_taxa']
    logger.info(f"Loaded parameters: {params}")

    # load the OGT and taxa information
    taxa = pd.read_parquet('./data/taxa.parquet', usecols=['taxid', 'ogt'])
    if taxa.isna().sum() > 0:
        raise ValueError("Found taxa with NaN values")
    else:
        pass
    
    # extract meso or thermo based on threshold
    logger.info("Labeling thermophiles")
    thermo_bools = taxa['ogt'].apply(lambda ogt: ogt > params['ogt_threshold'])
    taxa['thermophile_label'] = thermo_bools

    # make plot of OGT and remove it as a column
    # save plot
    fig, ax = plt.subplots(figsize=(5,5))
    ax.set_xlabel('OGT [C]')
    taxa['ogt'].plot.hist(bins=15, ax=ax)
    ax.vlines(x=[params['ogt_threshold']], ymin=0, ymax=ax.get_ylim()[1], colors=['r'])
    plt.savefig('./data/plots/ogt_hist.png', bbox_inches='tight', dpi=250)
    taxa = taxa.drop(columns=['ogt'])
    
    # save to file
    taxa.to_parquet('./data/taxa_thermophile_labels.parquet')

    # save metrics
    metrics = {}
    metrics['n_meso'] = int((taxa['thermophile_label'] == False).sum())
    metrics['n_therm'] = int((taxa['thermophile_label'] == True).sum())
    with open('./data/metrics/s1.0_metrics.yaml', "w") as stream:
        yaml_dump(metrics, stream)
    
    
