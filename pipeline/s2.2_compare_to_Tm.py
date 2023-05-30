"""Load datasets of melting temperatures, search for those proteins in the dataset
and see if they correspond to higher OGT on average.

Datasets:
 - FIREPROTDB: https://loschmidt.chemi.muni.cz/fireprotdb/
 - Meltome Atlas: https://meltomeatlas.proteomics.wzw.tum.de/

Steps for each dataset:
    - download data
    - get amino acid sequences and Tm
    - blast against our database of proteins to find near identical sequences
    - plot ogt vs Tm
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
from codecarbon import OfflineEmissionsTracker

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

FIREPROTDB_URL = "https://loschmidt.chemi.muni.cz/fireprotdb/v1/export?searchType=advanced&type=csv"
MELTOME_URL = "https://meltomeatlas.proteomics.wzw.tum.de/master_meltomeatlasapp/cross-species.csv"

RELOAD = False

# functions for uniprot api
def get_next_link(headers):
        if "Link" in headers:
            match = re_next_link.match(headers["Link"])
            if match:
                return match.group(1)

def get_batch(batch_url):
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(response.headers)

if __name__ == "__main__":

    # co2 tracking
    tracker = OfflineEmissionsTracker(
        project_name=f"s2.X",
        output_dir='./data/',
        country_iso_code='USA',
        region='Washington'
    ) 
    tracker.start()

    # spin up database
    tmpfile = tempfile.NamedTemporaryFile(dir='./tmp', mode='w+b', delete = True)
    tmpfile.close()
    logger.info("Creating database  at %s", tmpfile.name)
    con = ddb.connect(database=tmpfile.name, read_only=False)
    learn2therm.database.L2TDatabase._create_taxa_table(con, './data/')
    learn2therm.database.L2TDatabase._create_proteins_table(con, './data/')
    con.execute("CREATE UNIQUE INDEX taxa_primary ON taxa (taxid)")
    con.execute("CREATE UNIQUE INDEX prot_primary ON proteins (pid)")
    con.execute("CREATE INDEX prot_taxa_foreign ON proteins (taxid)")
    logger.info("Created learn2therm database")

    # FIREPROT DB
    ########################################################################
    logger.info("Downloading fireprot data")
    payload = {"searchData":{"type":"expr","key":"all","value":" ","checkOptions":[]},"filter":{"filterKey":"ddG","order":"asc"}}
    response = requests.post(FIREPROTDB_URL, json=payload)
    content = response.content.decode('utf-8')

    fireprot_df = pd.read_csv(io.StringIO(content))[['uniprot_id', 'tm', 'dTm', 'sequence']]
    logger.info(f"{len(fireprot_df)} raw data points from FireProtDB")

    # clean up the data  abit
    fireprot_df.dropna(inplace=True)
    logger.info(f"{len(fireprot_df)} fireprot data points after dropping NaNs")
    fireprot_df = fireprot_df.drop_duplicates(subset=['uniprot_id'])
    logger.info(f"{len(fireprot_df)} fireprot data points after dropping duplicate uniprot ids")
    fireprot_df['wTm'] = fireprot_df['tm'] - fireprot_df['dTm']
    logger.info(f"Tm data: {fireprot_df['wTm'].describe()}")
    fireprot_df.rename(columns={'uniprot_id': 'pid'}, inplace=True)

    # need to search out dataset of proteins for matches
    # since uniprot ids are not necessarily non redundant for >99% identity
    # we need to search by sequence
    our_proteins = con.execute("SELECT pid, protein_seq AS sequence FROM proteins").df()
    logger.info(f"Got {len(our_proteins)} proteins from our dataset")
    
    handler = learn2therm.blast.DiamondAlignmentHandler(
        seqs_A = fireprot_df,
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
    logger.info("Running blast on fireprot")
    results, metadata = handler.run()
    logger.info(f"Blast took {(time.time() - time0)/60} minutes")
    logger.info(f"Got {len(results)} results from blast for fireprot")

    # aggregate for the best result for each query
    results = results.sort_values('scaled_local_query_percent_id', ascending=False).drop_duplicates(subset=['query_id'])
    logger.info(f"Got {len(results)} results after taking best match for each query fireprot")
    logger.info(f"Results description: {results['scaled_local_query_percent_id'].describe()}")

    # merge with tm data
    fireprot_df = fireprot_df.merge(results, left_on='pid', right_on='query_id', how='inner')
    # drop those without 99% identity
    fireprot_df = fireprot_df[fireprot_df['scaled_local_query_percent_id'] >= .99]
    logger.info(f"Protein map: {fireprot_df[['query_id', 'subject_id']]}")
    fireprot_df = fireprot_df[['subject_id', 'wTm']]
    logger.info(f"Got {len(fireprot_df)} results after dropping those without 99% identity")
    # clear some memory
    del our_proteins

    # add wTm data
    con.execute("CREATE OR REPLACE TABLE tm_fpdb AS SELECT * FROM fireprot_df")
    logger.info("Added tm data")
    con.execute("CREATE INDEX tm_primary ON tm_fpdb (subject_id)")
    logger.info("added fpdb to database")

    # run a join on data to get ogt vs meting temp
    fpdb_data = con.execute("""
        SELECT tm_fpdb.wTm, taxa.temperature FROM taxa
        INNER JOIN proteins ON proteins.taxid = taxa.taxid
        INNER JOIN tm_fpdb ON tm_fpdb.subject_id = proteins.pid
    """).df()
    fpdb_data = fpdb_data.rename(columns={'temperature': 'OGT', 'wTm': 'Tm'})[['OGT', 'Tm']]
    fpdb_data['from'] = 'fireprot'
    logger.info(f"Got {len(fpdb_data)} data points with both OGT and Tm from FireProt")

    # Meltome Atlas
    ########################################################################
    if not RELOAD:
        meltome_df = pd.read_csv(MELTOME_URL)[['Protein_ID', 'meltPoint']]
        meltome_df.to_csv('./tmp/meltome_raw.csv')
    else:
        meltome_df = pd.read_csv('./tmp/meltome_raw.csv', index_col=0)
    
    logger.info(f"{len(meltome_df)} raw data points from Meltome Atlas")
    meltome_df = meltome_df.dropna()
    logger.info(f"{len(meltome_df)} meltome data points after dropping NaNs")
    # clean up the data  abit
    # get uniprot id
    def split_uniprot_id(id_):
        return id_.split("_")[0]
    meltome_df['Protein_ID'] = meltome_df['Protein_ID'].apply(split_uniprot_id)
    logger.info(f'Fixed uniprot ids: {meltome_df.head()}')
    logger.info(f"Total meltome data count, unique pids: {len(meltome_df)}, {len(meltome_df['Protein_ID'].unique())}")
    # drop duplicate uniprots
    meltome_groups = meltome_df.groupby('Protein_ID')
    logger.info(f"Variability in melting temperatures for uniprot ids in Meltome: {meltome_groups.std()['meltPoint'].describe()}")
    meltome_df = meltome_groups.mean().reset_index()
    logger.info(f"{len(meltome_df)} meltome data points after aggregating duplicate PIDs")
    logger.info(f"{meltome_df.head()}")
    # subset data by only proteins we do not have an exact match for
    unknown_meltome_pids = con.execute(
        """SELECT DISTINCT(Protein_ID) 
        FROM meltome_df 
        WHERE Protein_ID NOT IN (SELECT pid FROM proteins)""").df()['Protein_ID'].values
    known_meltome_pids = con.execute(
        """SELECT COUNT(DISTINCT(Protein_ID))
        FROM meltome_df
        WHERE Protein_ID IN (SELECT pid FROM proteins)""").df()
    logger.info(f"{known_meltome_pids} pids from meltone in our dataset, {len(unknown_meltome_pids)} unknown meltome proteins")
    # get the amino acid sequence to balst vs our database
    logger.info("Getting amino acid sequences for unkown meltome data by querying uniprot")

    if not RELOAD:
        # retrieve sequences
        unknown_meltome_seqs = learn2therm.uniprot.get_uniprot_sequences(
            list(unknown_meltome_pids)).rename(columns={'id': 'pid', 'seq': 'sequence'})
        unknown_meltome_seqs.to_csv("./tmp/unknown_meltome_seqs.csv")
    else:
        unknown_meltome_seqs = pd.read_csv("./tmp/unknown_meltome_seqs.csv")

    # run diamond to find near identical sequences
    our_proteins = con.execute("SELECT pid, protein_seq AS sequence FROM proteins").df()
    logger.info(f"Got {len(our_proteins)} proteins from our dataset")
    
    handler = learn2therm.blast.DiamondAlignmentHandler(
        seqs_A = unknown_meltome_seqs,
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
    logger.info("Running blast on meltome")
    results, metadata = handler.run()
    logger.info(f"Blast took {(time.time() - time0)/60} minutes")
    logger.info(f"Got {len(results)} results from blast for meltome")
    # clear some memory
    del our_proteins

    # aggregate for the best result for each query
    results = results.sort_values('scaled_local_query_percent_id', ascending=False).drop_duplicates(subset=['query_id'])
    logger.info(f"Got {len(results)} results after taking best match for each query for meltome unknown sequences")
    logger.info(f"Results description: {results['scaled_local_query_percent_id'].describe()}")

    # drop those without 99% identity
    results = results[results['scaled_local_query_percent_id'] >= .99]
    logger.info(f"Protein map: {results[['query_id', 'subject_id']]}")
    logger.info(f"Got {len(results)} meltome PID maps after dropping those without 99% identity")
    # replace pids with near identities
    # turn the map into a dict
    meltome_pid_map = dict(zip(results['query_id'], results['subject_id']))
    def mapper(meltome_pid):
        return meltome_pid_map.get(meltome_pid, meltome_pid)
    meltome_df['Protein_ID'] = meltome_df['Protein_ID'].apply(mapper)
    logger.info('Mapped some meltome PIDs to our PIDs')

    # add meltome data
    con.execute("CREATE OR REPLACE TABLE tm_meltome AS SELECT * FROM meltome_df")
    logger.info("Added tm data")
    con.execute("CREATE INDEX tm_primary_meltome ON tm_meltome (Protein_ID)")

    # run a join on data to get ogt vs meting temp
    meltome_data = con.execute("""
        SELECT tm_meltome.meltPoint, taxa.temperature FROM taxa
        INNER JOIN proteins ON proteins.taxid = taxa.taxid
        INNER JOIN tm_meltome ON tm_meltome.Protein_ID = proteins.pid
    """).df()
    meltome_data = meltome_data.rename(columns={'temperature': 'OGT', 'meltPoint': 'Tm'})[['OGT', 'Tm']]
    meltome_data['from'] = 'meltome'
    logger.info(f"Got {len(meltome_data)} data points with both OGT and T from meltome")
    
    # now get XXX et al data

    # make a plot of data
    data = pd.concat([meltome_data, fpdb_data], ignore_index=True)
    fig, ax = plt.subplots(figsize=(5,5))
    min_ = min(data['OGT'].min(), data['Tm'].min())
    max_ = max(data['OGT'].max(), data['Tm'].max())
    ax.plot([min_, max_], [min_, max_], color='black', linestyle='--')
    ax.set_xlim(min_, max_)
    ax.set_ylim(min_, max_)
    sns.scatterplot(data=data, x='OGT', y='Tm', hue='from', ax=ax)
    ax.set_xlabel("OGT [C]")
    ax.set_ylabel("Tm [C]")
    data.to_csv('./data/validation/tm/ogt_vs_tm.csv')

    # save the plot
    fig.savefig('./data/validation/tm/ogt_vs_tm_check.png', dpi=300, bbox_inches='tight')

    # save spearmans
    r, p = scipy.stats.spearmanr(data['OGT'], data['Tm'])

    # metrics
    co2 = float(tracker.stop())
    metrics = {
        'Tm_OGT_spearman_r': float(r),
        'Tm_OGT_spearman_p': float(p),
        'Tm_OGT_co2': co2,
    }
    with open('./data/validation/tm/metrics.yaml', 'w') as f:
        yaml_dump(metrics, f)
    
