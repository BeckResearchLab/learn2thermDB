"""Load a dataset of melting temperatures, search for those proteins in the dataset
and see if they correspond to higher OGT on average.
"""
import os
import shutil
import time
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
import learn2therm.blast

if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
    LOGLEVEL = getattr(logging, LOGLEVEL)
else:
    LOGLEVEL = logging.INFO
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'
# get the logger in subprocesses
logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')

FIREPROTDB_URL = "https://loschmidt.chemi.muni.cz/fireprotdb/v1/export?searchType=advanced&type=csv"

if __name__ == "__main__":

    # download the data
    logger.info("Downloading data")
    payload = {"searchData":{"type":"expr","key":"all","value":" ","checkOptions":[]},"filter":{"filterKey":"ddG","order":"asc"}}
    response = requests.post(FIREPROTDB_URL, json=payload)
    content = response.content.decode('utf-8')

    tm_df = pd.read_csv(io.StringIO(content))[['uniprot_id', 'tm', 'dTm', 'sequence']]
    logger.info(f"{len(tm_df)} raw data points from FireProtDB")

    # clean up the data  abit
    tm_df.dropna(inplace=True)
    logger.info(f"{len(tm_df)} data points after dropping NaNs")
    tm_df = tm_df.drop_duplicates(subset=['uniprot_id'])
    logger.info(f"{len(tm_df)} data points after dropping duplicate uniprot ids")
    tm_df['wTm'] = tm_df['tm'] - tm_df['dTm']
    logger.info(f"Tm data: {tm_df['wTm'].describe()}")
    tm_df.rename(columns={'uniprot_id': 'pid'}, inplace=True)

    # spin up database
    tmpfile = tempfile.NamedTemporaryFile(dir='./tmp', mode='w+b', delete = True)
    tmpfile.close()
    logger.info("Creating database  at %s", tmpfile.name)
    con = ddb.connect(database=tmpfile.name, read_only=False)
    learn2therm.database.L2TDatabase._create_taxa_table(con, './data/')
    learn2therm.database.L2TDatabase._create_proteins_table(con, './data/')
    logger.info("Created learn2therm database")

    # need to search out dataset of proteins for matches
    # since uniprot ids are not necessarily non redundant for >99% identity
    # we need to search by sequence
    our_proteins = con.execute("SELECT pid, protein_seq AS sequence FROM proteins ORDER BY RANDOM()").df()
    logger.info(f"Got {len(our_proteins)} proteins from our dataset")
    
    handler = learn2therm.blast.DiamondAlignmentHandler(
        seqs_A = tm_df,
        seqs_B = our_proteins,
        metrics = ['scaled_local_query_percent_id'],
        alignment_params = {
            'num_threads': 6,
            'sensitivity': 'fast',
            'iterate': False,              
            'global_ranking': None,           
            'gapopen': 11,
            'gapextend': 2,
            'matrix': 'BLOSUM62',
            'evalue': .0000001,
            'hsp_cov': 99,  
        }
    )
    time0 = time.time()
    logger.info("Running blast")
    results, metadata = handler.run()
    logger.info(f"Blast took {(time.time() - time0)/60} minutes")
    logger.info(f"Got {len(results)} results from blast")

    # aggregate for the best result for each query
    results = results.sort_values('scaled_local_query_percent_id', ascending=False).drop_duplicates(subset=['query_id'])
    logger.info(f"Got {len(results)} results after taking best match for each query")
    logger.info(f"Results description: {results['scaled_local_query_percent_id'].describe()}")

    # merge with tm data
    tm_df = tm_df.merge(results, left_on='pid', right_on='query_id', how='inner')
    # drop those without 99% identity
    tm_df = tm_df[tm_df['scaled_local_query_percent_id'] >= .99][['subject_id', 'wTm']]
    logger.info(f"Protein map: {tm_df[['query_id', 'subject_id']]}")
    logger.info(f"Got {len(tm_df)} results after dropping those without 99% identity")
    # clear some memory
    del our_proteins

    # add wTm data
    con.execute("CREATE OR REPLACE TABLE tm AS SELECT * FROM tm_df")
    logger.info("Added tm data")
    con.execute("CREATE UNIQUE INDEX taxa_primary ON taxa (taxid)")
    con.execute("CREATE UNIQUE INDEX prot_primary ON proteins (pid)")
    con.execute("CREATE INDEX prot_taxa_foreign ON proteins (taxid)")
    con.execute("CREATE INDEX tm_primary ON tm (subject_id)")
    logger.info("Created indexes")

    # run a join on data to get ogt vs meting temp
    data = con.execute("""
        SELECT tm.wTm, taxa.temperature FROM taxa
        INNER JOIN proteins ON proteins.taxid = taxa.taxid
        INNER JOIN tm ON tm.subject_id = proteins.pid
    """).df()
    con.close()
    os.remove(tmpfile.name)
    logger.info(f"Got {len(data)} data points with both OGT and Tm")

    # make a plot of data
    fig, ax = plt.subplots(figsize=(5,5))
    sns.regplot(data=data, x='temperature', y='wTm', ax=ax)
    ax.set_xlabel("OGT [C]")
    ax.set_ylabel("Tm [C]")

    # save the plot
    fig.savefig('./data/plots/ogt_vs_tm_check.png', dpi=300, bbox_inches='tight')

    # save spearmans
    r, p = scipy.stats.spearmanr(data['temperature'], data['wTm'])
    metrics = {
        'Tm_OGT_spearman_r': float(r),
        'Tm_OGT_spearman_p': float(p),
    }
    with open('./data/metrics/s2.Tm_metrics.yaml', 'w') as f:
        yaml_dump(metrics, f)
    
