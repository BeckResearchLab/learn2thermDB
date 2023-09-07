"""Load Engqvist's data of enzymes to OGT and check the overlap

Datasets:
 - https://zenodo.org/record/2539114
    - https://zenodo.org/record/2539114/files/enzyme_ogt_topt.tsv?download=1

Donwload the dataset, use uniprot ID to map to our dataset, then
get R2 of OGTs
"""
import os
import shutil
import time
import io
from yaml import safe_load as yaml_load
from yaml import dump as yaml_dump
import tempfile
import duckdb as ddb

import numpy as np
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt

from typing import List
import requests
import logging

import learn2therm.database
import learn2therm.utils
import learn2therm.blast
import learn2therm.uniprot

if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
    LOGLEVEL = getattr(logging, LOGLEVEL)
else:
    LOGLEVEL = logging.INFO
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'
# get the logger in subprocesses
logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')

ENGQVIST_DL = "https://zenodo.org/record/2539114/files/enzyme_ogt_topt.tsv?download=1"

if __name__ == "__main__":

    # load up our database
    con = ddb.connect('./data/database.ddb', read_only=True)
    
    # download enqvist's data to temporary directory
    logger.info(f"Downloading Engqvist's data from {ENGQVIST_DL}")
    
    response = requests.get(ENGQVIST_DL)
    if response.status_code != 200:
        logger.error(f"Failed to download data from {ENGQVIST_DL}")
        exit(1)

    with open('./tmp/enzyme_ogt_topt.tsv', 'wb') as file:
        file.write(response.content)

    logger.info(f"Engqvist's data.")

    # Add the data as a temporary table
    logger.info("Adding Engqvist's data as a temporary table")
    con.execute("CREATE TEMP TABLE engqvist AS SELECT * FROM read_csv_auto('./tmp/enzyme_ogt_topt.tsv', sep='\t')")
    logger.info(f'{con.execute("SELECT COUNT(*) FROM engqvist").fetchone()} examples from enqvist, some predicted')

    # inner join our OGTs with engqvist's OGTs
    logger.info("Joining our OGTs with Engqvist's OGTs")
    df = con.execute(
        """SELECT engqvist.ogt AS ogt_e, taxa.temperature AS ogt_u FROM taxa
            INNER JOIN proteins ON (taxa.taxid=proteins.taxid)
            INNER JOIN engqvist ON (engqvist.uniprot_id=proteins.pid)
            WHERE engqvist.ogt_source='experimental'
    """).df().dropna()
    logger.info(f"Data overlap: {len(df)}")
    # Compute R2 and make a plot

    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(df['ogt_e'], df['ogt_u'])

    # Calculate R^2
    r2_value = r_value**2
    logger.info(f"Computed R^2 value: {r2_value}")

    # Making a plot
    plt.scatter(df['ogt_e'], df['ogt_u'], label='Data')
    plt.plot(df['ogt_e'], intercept + slope * df['ogt_e'], color='red', label=f'Linear Fit R^2={r2_value:.2f}')
    plt.xlabel('Engqvist OGT')
    plt.ylabel('Our OGT')
    plt.legend()
    plt.title('Comparison of OGTs')

    if not os.path.exists('./data/validation/engqvist'):
        os.makedirs('./data/validation/engqvist')
    plt.savefig('./data/validation/engqvist/ogt_comparison.png', bbox_inches='tight', dpi=600)

    metrics = {
        'r2_value': float(r2_value),
        'slope': float(slope),
        'intercept': float(intercept),
        'p_value': float(p_value),
        'std_err': float(std_err)
    }
    with open('./data/validation/engqvist/metrics.yaml', 'w') as f:
        yaml_dump(metrics, f)
    
