"""Runs blast and reports metrics for every pair of thermophiles, mesophiles
"""
import ast
import fcntl
import logging
import os

from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import pandas as pd
from yaml import dump as yaml_dump
from yaml import safe_load as yaml_load

import learn2therm.blast
import learn2therm.io
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

OUTFILENAME = './data/taxa_pairs/pairwise_16s_blast.csv'

def apply_blast_metric_and_append_pairwise(
    blast_record,
    metrics: List[str],
):
    """Compute metrics for a single blast record and save to file with locking.
    
    Metrics returns dataframe of (query index, subject_index, metric value) for each metric
    """
    logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')
    # compute all metrics
    metric_handler = learn2therm.blast.BlastMetrics(blast_record)
    df_outputs = [metric_handler.compute_metric(m) for m in metrics]
    # each df has query and subject id and one metric
    # join them all
    joined_df = df_outputs[0].set_index(['query_id', 'subject_id'])
    for df in df_outputs[1:]:
        joined_df = joined_df.join(df.set_index(['query_id', 'subject_id']))
    out = joined_df.reset_index()
    # ensure ordering
    out = out[['query_id', 'subject_id']+metrics]
    logger.debug(f"Adding {len(out)} alignments for thermophile {metric_handler.qid}")
    # save to file
    with open(OUTFILENAME, "a") as g:
        fcntl.flock(g, fcntl.LOCK_EX)
        out.to_csv(g, mode='a',index=False, header=False)
        fcntl.flock(g, fcntl.LOCK_UN)
    return

if __name__ == '__main__':
    # connect to log file
    logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')

    # DVC tracked parameters
    with open("./params.yaml", "r") as stream:
        params = yaml_load(stream)['get_16s_blast_scores']
    logger.info(f"Loaded parameters: {params}")

    # load the thermo meso labels
    labels = pd.read_csv('./data/taxa/labels.csv', index_col=0,  usecols=[0,2])['thermophile_label']

    # get the tax indexes of themophiles and mesophiles
    thermo_indexes = list(labels[labels == True].index)
    meso_indexes = list(labels[labels == False].index)
    logger.info(f"Found taxa indexes for {thermo_indexes} thermophiles and {meso_indexes} mesophiles")

    # run blast
    # ###########
    # these iterators just mean that we don't have to load all sequences into memory
    thermophile_iterator = learn2therm.io.csv_id_seq_iterator(
        './data/taxa/16s_rRNA.csv', seq_col='seq_16srRNA', id_filter=thermo_indexes, index_col=0)
    mesophile_iterator = learn2therm.io.csv_id_seq_iterator(
        './data/taxa/16s_rRNA.csv', seq_col='seq_16srRNA', id_filter=meso_indexes, index_col=0)

    # sample if we neeed to
    if params['n_sample'] != None:
        thermophile_iterator = [next(thermophile_iterator) for i in range(params['n_sample'])]
        mesophile_iterator = [next(mesophile_iterator) for i in range(params['n_sample'])]
        logger.info(f"Downsample to each of {params['n_sample']} thermo and meso")

    # this context manager creates temporary files for BLAST inputs and outputs, again to avoid OOM
    with learn2therm.blast.BlastFiles(thermophile_iterator, mesophile_iterator, dbtype='nucl') as (query_fname, subject_fname, out_fname):
        logger.info('Running blast...')
        NcbiblastnCommandline(
            query=query_fname,
            db=subject_fname,
            outfmt=5,
            out=out_fname,
            max_target_seqs=10000000, # very large so we do not throw out any pairs. will have to increase if there is more than this num of mesos
            perc_identity=0.0, # no cutoff to lose data
            evalue=10000000, # very large so we do not lose any hits
            # the rest are tunable params
            word_size=params['word_size'],
            gapopen=params['gapopen_penalty'],
            gapextend=params['gapextend_penalty'],
            reward=params['reward'],
            penalty=params['penalty'],
            ungapped=params['ungapped'],
            num_threads=params['num_threads']
        )()
        logger.info('Blast complete. Parsing and saving metrics.')

        # create empty file
        file = open(OUTFILENAME, 'w')
        file.write(f"thermo_index,meso_index,{','.join(params['blast_metrics'])}\n")
        file.close()

        # blast record io
        xml_f = open(out_fname, 'r')
        blast_result_records = NCBIXML.parse(xml_f)

        # compute and save
        for record in blast_result_records:
            apply_blast_metric_and_append_pairwise(record, params['blast_metrics'])


        # save some plots
        import matplotlib.pyplot as plt
        for i, m in enumerate(params['blast_metrics']):
            scores = pd.read_csv(OUTFILENAME, usecols=[i+2])
            scores = scores[m]

            fig, ax = plt.subplots(figsize=(5,5))
            ax.set_xlabel(m)
            scores.plot.hist(bins=50, ax=ax)
            plt.savefig(f'./data/plots/blast/blast_hist_{m}.png', bbox_inches='tight', dpi=250)
        
        # save metrics
        metrics = {'percent_full_pairwise_16s_blast': float(len(pd.read_csv(OUTFILENAME, usecols=[0]))/(len(thermo_indexes)*len(meso_indexes)))}
        with open('./data/metrics/s1.2_metrics.yaml', "w") as stream:
            yaml_dump(metrics, stream)
    
    
