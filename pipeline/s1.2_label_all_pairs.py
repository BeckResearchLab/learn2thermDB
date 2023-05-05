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

ALIGNMENT_SCORE_DIR = './data/taxa_pairs/alignment/'
LABELS_DIR = './data/taxa_pairs/pair_labels/'

def process_file(input_file, params, output_dir, output_filename):
    metric_vals = pd.read_parquet(input_file)
    pairwise_score_cols = list(metric_vals.columns)
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
    for mname, thresh in params['blast_metric_thresholds'].items():
        this_metric_vals = metric_vals[mname]
        if this_metric_vals.isna().any():
            raise ValueError(f'Some metrics came back Nan, this should not be')
        if thresh['greater'] == True:
            mask = this_metric_vals > thresh['thresh']
        else:
            mask = this_metric_vals < thresh['thresh']
        logger.info(f"Only considering pairs with {mname} {'>' if thresh['greater'] else '<'} {thresh['thresh']}")
        masks.append(mask)
    mask = pd.DataFrame(masks).all(axis=0)
    mask.name='is_pair'

    # save the data
    mask = pd.DataFrame(mask)
    mask.to_parquet(output_dir+output_filename, index=True)

    # save metrics
    metrics = {
        'num_taxa_pairs_conservative': int(mask['is_pair'].sum()),
        'taxa_pair_found_ratio': float(mask['is_pair'].sum()/len(mask))
    }
    return metrics

if __name__ == '__main__':
    # connect to log file
    logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')

    # DVC tracked parameters
    with open("./params.yaml", "r") as stream:
        params = yaml_load(stream)['label_all_pairs']
    logger.info(f"Loaded parameters: {params}")

    # check directory exists
    if not os.path.exists(LABELS_DIR):
        os.makedirs(LABELS_DIR)

    # execute function on all files in the directory
    files = os.listdir(ALIGNMENT_SCORE_DIR)

    all_metrics = []
    for file in files:
        if not file.endswith('.parquet'):
            continue
        logger.info(f"Processing {file}")
        metrics = process_file(ALIGNMENT_SCORE_DIR+file, params, LABELS_DIR, file)
        all_metrics.append(metrics)
        if params['dev_only_one_file']:
            break
    
    agg_metrics = {
        'num_taxa_pairs_conservative': sum(m['num_taxa_pairs_conservative'] for m in all_metrics),
        'taxa_pair_found_ratio': sum(m['taxa_pair_found_ratio'] for m in all_metrics) / len(all_metrics)
    }
   
    with open('./data/metrics/s1.2_metrics.yaml', "w") as stream:
        yaml_dump(agg_metrics, stream)