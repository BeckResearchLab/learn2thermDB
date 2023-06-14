"""Compare out alignment scores to haits alignment scores."""

from yaml import safe_load as yaml_load
from yaml import dump as yaml_dump
import os
import sys
import time
import tempfile
import pandas as pd
import numpy as np
import duckdb as ddb
import logging
import tqdm
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_style('whitegrid')
sns.set_context('paper')

import learn2therm.utils

if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
    LOGLEVEL = getattr(logging, LOGLEVEL)
else:
    LOGLEVEL = logging.INFO
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'
logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')

if __name__ == '__main__':
    # load hait alignment data
    hait_data = pd.read_csv('./data/validation/hait_aligned_scores.csv', index_col=0)

    # connect to out database
    con = ddb.connect(database='./data/database.ddb', read_only=True)

    # get the output ready
    os.makedirs('./data/validation/hait_alignment', exist_ok=True)

    # make plot for each alignment metric as a comparison
    # compute the percentiles in haits pairs
    percentiles = (5, 10, 25, 50)
    hait_percentiles = {}
    remaining_data = {}
    queries = {
        'normalized_bit_score': "(2 * bit_score::FLOAT/(query_align_len::FLOAT + subject_align_len::FLOAT))::FLOAT",
        'scaled_local_symmetric_percent_id': "scaled_local_symmetric_percent_id",
        'coverage': "(query_align_cov::FLOAT + subject_align_cov::FLOAT)/2"
    }
    for metric, query in queries.items():
        if metric == 'normalized_bit_score':
            hait_data[metric] = hait_data.apply(lambda row: 2 * row['bit_score'] / (row['query_align_len'] + row['subject_align_len']), axis=1)
        elif metric == 'coverage':
            hait_data['coverage'] = hait_data.apply(lambda row: (row['query_align_cov'] + row['subject_align_cov']) / 2, axis=1)
        logger.info(f"Doing {metric}...")
        logger.info(f"Hait data for {metric}: {hait_data[metric]}")
        # now compute percentiles
        hait_percentiles[metric] = {p: hait_data[metric].quantile(p/100) for p in percentiles}
        logger.info(f"Hait percentiles for {metric}: {hait_percentiles[metric]}")

        our_values = con.execute(f"SELECT {query} FROM pairs").df().iloc[:,0].values
        logger.debug(f"len(our_values) = {len(our_values)}")
        
        # add aggregate dataframe so that we can plot it
        df1_ = pd.DataFrame({metric: hait_data[metric].values, 'source': 'Hait'})
        df2_ = pd.DataFrame({metric: our_values, 'source': 'learn2therm'})
        df = pd.concat([df1_, df2_], axis=0, ignore_index=True) 
        # plot
        fig, ax = plt.subplots(figsize=(5,5))
        sns.histplot(data=df, x=metric, bins=25, hue='source', ax=ax, stat='probability', common_norm=False, common_bins=True)
        ax.set_xlabel(metric)
        ax.set_ylabel('probability')
        fig.savefig(f'./data/validation/hait_alignment/hait_vs_learn2therm_{metric}.png', dpi=300, bbox_inches='tight')
        logger.info(f"Saved plot to ./validation/hait_alignment/hait_vs_learn2therm_{metric}.png")
        remaining_data[metric] = {
            p: con.execute(f"SELECT COUNT(*) FROM pairs WHERE {query}>{hait_percentiles[metric][p]}").df().iloc[0,0] for p in percentiles
        }
    remaining_data['all'] = {}
    for p in percentiles:
        aggregate_query = "SELECT COUNT(*) FROM pairs WHERE"
        for metric, query in queries.items():
            aggregate_query += f" {query} > {hait_percentiles[metric][p]} AND"
        aggregate_query = aggregate_query[:-4]
        remaining_data['all'][p] = con.execute(aggregate_query).df().iloc[0,0]
    
    # now make plots of the remaining data based on thresholds
    fig, axes = plt.subplots(1, len(queries)+1, figsize=(20,5), sharey=True)
    for i, (metric, remaining_data_subdict) in enumerate(remaining_data.items()):
        ax = axes[i]
        if metric == 'all':
            x = list(percentiles)
        else:
            x = list(hait_percentiles[metric].values())
        y = list(remaining_data_subdict.values())
        ax.set_ylim(0, 1.1*max(y))
        logger.info(f"{x}, {y}")
        sns.lineplot(x=x, y=y, ax=ax, label=metric)
        if metric == 'all':
            ax.set_xlabel('Hait percentile all metrics')
        else:
            ax.set_xlabel(metric)
        ax.set_ylabel('number of pairs remaining')
    fig.savefig(f'./data/validation/hait_alignment/remaining_data_after_hait_metric_filter.png', dpi=300, bbox_inches='tight')
    
    




        



    

    