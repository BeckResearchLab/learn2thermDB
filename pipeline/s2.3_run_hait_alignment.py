"""Compute the alignment metrics for hait pairs using identical parameters as
the alignment metrics for the full dataset.

See how the resulting metrics compare.
"""

from yaml import safe_load as yaml_load
from yaml import dump as yaml_dump
import os
import sys
import time
import tempfile
import pandas as pd
import logging
import tqdm


import learn2therm.utils
import learn2therm.blast

if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
    LOGLEVEL = getattr(logging, LOGLEVEL)
else:
    LOGLEVEL = logging.INFO
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'
logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')


if __name__ == '__main__':
    # load params
    with open("./params.yaml", "r") as stream:
        params = yaml_load(stream)['get_protein_blast_scores']

    # get aligner handler class
    alignment_method = params['method']
    try:
        aligner_params = params[f"method_{alignment_method}_params"]
    except:
        raise ValueError(f"specified method {alignment_method} for alignment but found no parameters for it.")

    # manually override minimum criteria so that we get all of the results back
    if 'qcov_hsp_perc' in aligner_params:
        aligner_params['qcov_hsp_perc'] = 1
    elif 'hsp_cov' in aligner_params:
        aligner_params['hsp_cov'] = 1
    aligner_params['evalue'] = .01
    
    if alignment_method == 'blast':
        Aligner = getattr(learn2therm.blast, 'BlastAlignmentHandler')
    elif alignment_method == 'diamond':
        Aligner = getattr(learn2therm.blast, 'DiamondAlignmentHandler')
    else:
        raise ValueError(f"Unknown alignment method {alignment_method}")

    # Parse sequences from Hait into the correct format for the aligner
    hait_data_raw = pd.read_csv('./data/validation/hait_pairs.csv', index_col=0)
    assert hait_data_raw.dropna().shape[0] == hait_data_raw.shape[0]
    logger.info(f"Got {len(hait_data_raw)} pairs from Hait et al.")
    thermo_seqs = hait_data_raw[['thermo_pid', 'thermo_seq']].rename(columns={'thermo_pid':'pid', 'thermo_seq':'sequence'}).drop_duplicates(subset='pid')
    meso_seqs = hait_data_raw[['meso_pid', 'meso_seq']].rename(columns={'meso_pid':'pid', 'meso_seq':'sequence'}).drop_duplicates(subset='pid')
    logger.info(f"Unique sequences: {len(thermo_seqs)} thermo sequences and {len(meso_seqs)} meso sequences")

    # ready the aligner
    # thermo vs meso
    aligner = Aligner(thermo_seqs, meso_seqs, metrics=params['blast_metrics'],alignment_params=aligner_params)
    logger.info(f"Running alignment on {len(thermo_seqs)} pairs")
    align_results, meta_data = aligner.run()
    # the outputs contain alignments against off-pairs as well, filter back to only the actual Hait pairs
    align_results = align_results.rename(columns={'query_id':'thermo_pid', 'subject_id':'meso_pid'})
    logger.info(f"Got {len(align_results)} alignments, filtering to only Hait pairs")
    hait_align_results = hait_data_raw.merge(align_results, on=['thermo_pid', 'meso_pid'], how='left')
    logger.info(f"Got {len(hait_align_results)} alignments for Hait pairs")
    # save the results
    hait_align_results.to_csv('./data/validation/hait_aligned_scores.csv')