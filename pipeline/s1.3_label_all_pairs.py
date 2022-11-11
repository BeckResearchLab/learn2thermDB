"""Label the thermo-meso pairs that we will run BLASTP on
"""
import ast
from curses import intrflush
from distutils.errors import DistutilsClassError
import logging
import os

from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from joblib import Parallel, delayed
import pandas as pd
from yaml import dump as yaml_dump
from yaml import safe_load as yaml_load

import learn2therm.blast
import learn2therm.io
import learn2therm.blast
import learn2therm.utils

from typing import List

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
        params = yaml_load(stream)['label_all_pairs']
    logger.info(f"Loaded parameters: {params}")

    # get the columns in the 16s blast file. This will tell us what metrics we actually have available to blast
    pairwise_score_cols = pd.read_csv('./data/taxa_pairs/pairwise_16s_blast.csv', nrows=1).columns

    # the form of metric threshold params is very specific, check and parse them
    if type(params['blast_metric_thresholds']) != dict:
        raise ValueError(f"`blast_metric_thresholds` params should be a dict of (name, dict(thres: float, greater: bool))")
    for mname, vals in params['blast_metric_thresholds'].items():
        try:
            assert type(mname) == str
            assert type(vals) == dict
            assert 'thresh' in vals and 'greater' in vals
            assert type(vals['thresh']) == float
            assert type(vals['greater']) == bool
        except AssertionError:
            raise ValueError(f"`blast_metric_thresholds` params should be a dict of (name, dict(thres: float, greater: bool))")
        if mname not in pairwise_score_cols:
            raise ValueError(f"Specified threshold on {mname} but that metric was not computed.")

    # for each specified threshold, create a boolean mask
    masks = []
    for mname, vals in params['blast_metric_thresholds'].items():
        column_index = pairwise_score_cols.get_loc(mname)
        if type(column_index) != int:
            raise ValueError(f'Something went wrong! The position of column {mname} in the blast metrics file is not an integer')
        metric_vals = pd.read_csv('./data/taxa_pairs/pairwise_16s_blast.csv', usecols=[column_index])[mname]
        if metric_vals.isna().any():
            raise ValueError(f'Some metrics came back Nan, this should not be')
        if vals['greater'] == True:
            mask = metric_vals > vals['thresh']
        else:
            mask = metric_vals < vals['thresh']
        logger.info(f"Only considering pairs with {mname} {'>' if vals['greater'] else '<'} {vals['thresh']}")
        masks.append(mask)
    mask = pd.DataFrame(masks).all(axis=0)
    mask.name='is_pair'
    
    # save the data
    mask.to_csv('./data/taxa_pairs/pair_labels.csv')
    
    # save metrics
    metrics = {
        'num_pairs_conservative': int(mask.sum()),
        'true_pair_ratio': float(mask.sum()/len(mask))
    }
    with open('./data/metrics/s1.3_metrics.yaml', "w") as stream:
        yaml_dump(metrics, stream)