"""Sample some pairs in order to run structural alignment on.

Alignment is too expensive to run on millions of pairs, so we must take samples from the dataset.
To be condusive to downtream analysis, this sample is not taken randomly but rather
taken to appriximate a uniform distribution in some metric such as coverage. 
"""
import os
import shutil
import time
import io
from yaml import safe_load as yaml_load
from yaml import dump as yaml_dump
import duckdb as ddb

import numpy as np
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('whitegrid')
sns.set_context('paper')

from typing import List
import logging

import learn2therm.utils

if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
    LOGLEVEL = getattr(logging, LOGLEVEL)
else:
    LOGLEVEL = logging.INFO
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'
# get the logger in subprocesses
logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')

def uniform_sample(df, column, M, num_bins=10):
    """Sample M data points from dataframe to approximate a uniform distribution in column."""
    # Compute the histogram
    counts, bins = np.histogram(df[column], bins=num_bins)
    
    # Compute the number of samples to draw from each bin
    samples_per_bin = M // num_bins
    
    # Initialize the sample DataFrame
    sample_df = pd.DataFrame()
    
    # For each bin...
    for i in range(num_bins):
        # Get the data points in this bin
        bin_data = df[(df[column] >= bins[i]) & (df[column] < bins[i+1])]
        
        # If there are not enough data points in this bin, take all of them
        if len(bin_data) <= samples_per_bin:
            sample_df = pd.concat([sample_df, bin_data])
        else:
            # Otherwise, sample without replacement
            bin_sample = bin_data.sample(n=samples_per_bin, replace=False)
            sample_df = pd.concat([sample_df, bin_sample])
    
    # If we didn't get enough samples (because some bins had too few data points), 
    # sample the remaining from the entire data
    if len(sample_df) < M:
        remaining = M - len(sample_df)
        remaining_sample = df.sample(n=remaining, replace=False)
        sample_df = pd.concat([sample_df, remaining_sample])
    
    return sample_df

if __name__ == '__main__':

    # load the params
    with open('./params.yaml') as f:
        params = yaml_load(f)['sample_data_for_structure']

    # get the per metric sample sizes
    sub_sample_size = int(params['sample_size']/len(params['metrics']))

    # connect to the database
    con = ddb.connect('./data/database.ddb', read_only=True)

    samples = []
    # for each metric, get a sample of pairs
    df = con.execute(f"""
        SELECT thermo_pid, meso_pid, p2.pdb_id AS meso_pdb, p1.pdb_id AS thermo_pdb, {','.join(params['metrics'])} 
        FROM pairs INNER JOIN proteins AS p1 ON (p1.pid=pairs.thermo_pid)
        INNER JOIN proteins AS p2 ON (p2.pid=pairs.meso_pid)""").df()
    new_columns = ['thermo_pid', 'meso_pid', 'meso_pdb', 'thermo_pdb'] + params['metrics']
    df.columns = new_columns
    for metric in params['metrics']:
        logger.info(f"Sampling {sub_sample_size} pairs by {metric}")
        df = uniform_sample(df, metric, sub_sample_size)
        samples.append(df)
    sample = pd.concat(samples, ignore_index=True)
    sample = sample.drop_duplicates(subset=['thermo_pid', 'meso_pid'])
    logger.info(f"Sampled {len(sample)} pairs total")

    # save the sample
    if not os.path.exists('./data/validation/structure'):
        os.makedirs('./data/validation/structure')
    sample.to_csv('./data/validation/structure/sample_l2t_data.csv')

