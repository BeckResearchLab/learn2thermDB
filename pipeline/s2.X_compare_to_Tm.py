"""Load a dataset of melting temperatures, search for those proteins in the dataset
and see if they correspond to higher OGT on average.
"""
import os
import shutil
import io
from yaml import safe_load as yaml_load
from yaml import dump as yaml_dump
import pathlib
import tempfile
import torch
import duckdb as ddb

import numpy as np
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('whitegrid')
sns.set_context('paper')

from typing import List
import requests
import logging

import learn2therm.database
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

MELTOME_URL = "https://meltomeatlas.proteomics.wzw.tum.de/master_meltomeatlasapp/session/11bb3ca3b855652d324584b78fe7dfc7/download/downloadData?w="

if __name__ == "__main__":

    # download the data
    # s = requests.get(MELTOME_URL).content
    # content = s.decode('utf-8')
    # print(content)
    # meltome_df = pd.read_csv(io.StringIO(content))[['ProteinID', 'meltPoint']]
    meltome_df = pd.read_csv('./tmp/cross-species.csv')[['Protein_ID', 'meltPoint']]
    logger.info(f"{len(meltome_df)} raw data points from meltome atlas")

    # clean up the data  abit
    meltome_df.dropna(inplace=True)
    logger.info(f"{len(meltome_df)} data points after dropping NaNs")
    # get uniprot id
    def split_uniprot_id(id_):
        return id_.split("_")[0]
    meltome_df['Protein_ID'] = meltome_df['Protein_ID'].apply(split_uniprot_id)

    # spin up database
    tmpfile = tempfile.NamedTemporaryFile(dir='./tmp', mode='w+b', delete = True)
    tmpfile.close()
    con = ddb.connect(database=tmpfile.name, read_only=False)
    learn2therm.database.L2TDatabase._create_taxa_table(con, './data/')
    learn2therm.database.L2TDatabase._create_proteins_table(con, './data/')
    logger.info("Created learn2therm database")
    # add meltome data
    con.execute("CREATE OR REPLACE TABLE meltome AS SELECT * FROM meltome_df")
    logger.info("Added meltome data")
    con.execute("CREATE UNIQUE INDEX taxa_primary ON taxa (taxid)")
    con.execute("CREATE UNIQUE INDEX prot_primary ON proteins (pid)")
    con.execute("CREATE INDEX prot_taxa_foreign ON proteins (taxid)")
    con.execute("CREATE INDEX meltome_primary ON meltome (Protein_ID)")
    logger.info("Created indexes")

    # run a join on data to get ogt vs meting temp
    data = con.execute("""
        SELECT meltome.meltPoint, taxa.temperature FROM taxa
        INNER JOIN proteins ON proteins.taxid = taxa.taxid
        INNER JOIN meltome ON meltome.Protein_ID = proteins.pid
    """).df()
    con.close()
    os.remove(tmpfile.name)
    logger.info(f"Got {len(data)} data points with both OGT and Tm")

    # make a plot of data
    fig, ax = plt.subplots(figsize=(5,5))
    sns.regplot(data=data, x='temperature', y='meltPoint', ax=ax)
    ax.set_xlabel("OGT [C]")
    ax.set_ylabel("Tm [C]")

    # save the plot
    fig.savefig('./data/plots/ogt_vs_tm_check.png', dpi=300, bbox_inches='tight')

    # save spearmans
    r, p = scipy.stats.spearmanr(data['temperature'], data['meltPoint'])
    metrics = {
        'Tm_OGT_spearman_r': float(r),
        'Tm_OGT_spearman_p': float(p),
    }
    with open('./data/metrics/s2.Tm_metrics.yaml', 'w') as f:
        yaml_dump(metrics, f)
    
