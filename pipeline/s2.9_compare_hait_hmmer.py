"""Compare labels from hmmer on our data vs Hait et al.
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

if __name__ == "__main__":

    # get l2t data
    con = ddb.connect('./data/database.ddb', read_only=True)
    con.execute("CREATE TEMP TABLE labels AS SELECT * FROM read_parquet('./data/validation/hmmer/hmmer_labels/*.parquet')")
    
    l2t_data = con.execute("""SELECT (pairs.query_align_cov+pairs.subject_align_cov)/2.0, labels.score, labels.bool_confirmed
        FROM pairs INNER JOIN labels ON pairs.meso_pid=labels.meso_pid AND pairs.thermo_pid=labels.thermo_pid
    """).df().rename(
        columns={'score': 'jaccard_score', 'bool_confirmed': 'functional', '((pairs.query_align_cov + pairs.subject_align_cov) / 2.0)': 'coverage'})

    # get hait data
    hait_data = pd.read_csv('./data/validation/hmmer/hait_scores.csv', index_col=0)[['jaccard_score', 'functional']]
    logger.info(f'Loaded {hait_data.shape[0]} rows from Hait et al. and {l2t_data.shape[0]} rows from learn2therm.')

    # compute some metrics about finding at least one
    metrics = {}
    metrics['l2t_not_in_Pfam'] = float(l2t_data['functional'].isna().mean())
    metrics['hait_not_in_Pfam'] = float(hait_data['functional'].isna().mean())
    
    # add score metrics for those we did find
    metrics['l2t_jaccard_mean'] = float(l2t_data['jaccard_score'].mean())
    metrics['hait_jaccard_mean'] = float(hait_data['jaccard_score'].mean())
    # compute a t statistic
    t, p = scipy.stats.ttest_ind(l2t_data['jaccard_score'], hait_data['jaccard_score'], alternative='less',equal_var=False, nan_policy='omit')
    logger.info(f"t, p: {t}, {p} full")
    t, p = scipy.stats.ttest_ind(l2t_data[l2t_data['coverage']>0.95]['jaccard_score'], hait_data['jaccard_score'], alternative='less',equal_var=False, nan_policy='omit')
    logger.info(f"t, p: {t}, {p} 95+")
    metrics['t_pvalue_base'] = float(p)
    metrics['t_pvalue_95'] = float(p)
    # make historgram of jacccard for those we did find vs the 2 datasets
    fig, ax = plt.subplots()
    # first defined bins
    # bins = list(np.linspace(0, 1, 10))
    # bins.append(1.0)
    # plot l2t to the left and hait to the right so that bins to not overlap
    ax.hist(
        [l2t_data['jaccard_score'], l2t_data[l2t_data['coverage']>0.95]['jaccard_score'], hait_data['jaccard_score']],
        bins=10,
        alpha=1.0, 
        label=[f'learn2therm N={len(l2t_data)}',f'learn2therm 95+ N={int((l2t_data["coverage"]>0.95).sum())}', f'Hait et al. N={len(hait_data)}'], density=True)
    ax.set_xlabel('Jaccard score of Pfam annotations')
    ax.set_ylabel('Density')
    ax.legend()
    fig.savefig('./data/validation/hmmer/compare_jaccard_hist.png', dpi=300, bbox_inches='tight')

    # save metrics
    with open('./data/validation/hmmer/s2.9_metrics.yaml', 'w') as f:
        yaml_dump(metrics, f)

    
    