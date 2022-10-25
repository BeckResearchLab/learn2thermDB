"""Produce a very large file of pairwise blast scores for each protein in meso and thermo pairs
"""
import ast
import fcntl
import logging
import time
import os

from codecarbon import EmissionsTracker
from Bio.Blast.Applications import NcbiblastpCommandline
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
LOGNAME = ''
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'

OUTFILENAME = './data/taxa_pairs/pairwise_protein_blast.csv'
PROTEIN_SEQ_FILE = './data/taxa/proteins.csv'

def apply_blast_metric_and_append_pairwise(
    blast_record,
    metrics: List[str],
):
    """Compute metrics for a single blast record and save to file with locking.
    
    Metrics returns dataframe of (query index, subject_index, metric value) for each metric

    Returns
    -------
    number of protein pairs with metrics
    """
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
    # save to file
    with open(OUTFILENAME, "a") as g:
        fcntl.flock(g, fcntl.LOCK_EX)
        out.to_csv(g, mode='a',index=False, header=False)
        fcntl.flock(g, fcntl.LOCK_UN)
    return len(out)

def run_blast_and_apply_blast_metric_and_append_pairwise_one_taxa_pair(
    thermo_taxa_index: int,
    meso_taxa_index: int,
    protein_index_map: pd.Series,
    metrics: List[str]
):
    """Filter proteins to only the proteins of a particular thermo-meso pair, and run blast on them.
    
    Blast metrics are applied to the results and stored with indexes of proteins and metrics.

    Parameters
    ----------
    *index : int
        index of taxa in the taxa list file
    protein_index_map : pd.Series
        index is protein index, values are associated taxa index
    metrics : list of str
        names of blast metrics to apply to blast hits
    """
    # start logger and carbon tracker
    logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='a')

    # get list of protein indexes we want to blast on
    meso_protein_indexes = protein_index_map[protein_index_map == meso_taxa_index].index
    thermo_protein_indexes = protein_index_map[protein_index_map == thermo_taxa_index].index
    logger.info(f"BLASTing {len(thermo_protein_indexes)} thermo proteins against {len(meso_protein_indexes)} meso sequences for pair {(thermo_taxa_index,meso_taxa_index)}")

    # create iterators of protein sequences to give to blast input files
    # this file does not have an index column, we have been assuming that the position is the index, so do not
    # pass index col to kwargs
    meso_iter = learn2therm.io.csv_id_seq_iterator(PROTEIN_SEQ_FILE, seq_col="protein_seq", id_filter=meso_protein_indexes, sep=';')
    thermo_iter = learn2therm.io.csv_id_seq_iterator(PROTEIN_SEQ_FILE, seq_col="protein_seq", id_filter=thermo_protein_indexes, sep=';')

    time0 = time.time()
    with learn2therm.blast.BlastFiles(thermo_iter, meso_iter, dbtype='prot') as (query_fname, subject_fname, out_fname):
        logger.debug(f"Pair {(thermo_taxa_index,meso_taxa_index)} file transfer and DB creation finished, took {(time.time()-time0)/60} minutes")
        logger.debug(f"Pair {(thermo_taxa_index,meso_taxa_index)} query and subject files: {query_fname, subject_fname}")
        NcbiblastpCommandline(
            query=query_fname,
            db=subject_fname,
            outfmt=5,
            out=out_fname,
            max_target_seqs=10000000, # very large so we do not throw out any pairs. will have to increase if there is more than this num of mesos
            evalue=10000000, # very large so we do not lose any hits
            # the rest are tunable params
            matrix=params['matrix'],
            word_size=params['word_size'],
            gapopen=params['gapopen_penalty'],
            gapextend=params['gapextend_penalty'],
            threshold=params['word_score_thresh'],
            ungapped=params['ungapped'],
        )()
        time1 = time.time()
        logger.info(f"BLASTing complete for pair {(thermo_taxa_index,meso_taxa_index)}, took {(time1-time0)/60} minutes")

        # blast record io
        xml_f = open(out_fname, 'r')
        blast_result_records = NCBIXML.parse(xml_f)

        # compute mtrics and track how many hits we have
        hits = 0
        for record in blast_result_records:
            hits += apply_blast_metric_and_append_pairwise(record, metrics)
        time2=time.time()
        logger.info(f"Computed metrics for pair {(thermo_taxa_index,meso_taxa_index)}, took {(time2-time1)/60} minutes")
    
    # record carbon emissions and output
    logger.info(f"Roundtrip execution for pair {(thermo_taxa_index,meso_taxa_index)} took {(time2-time0)/60}")
    return {'ids': (thermo_taxa_index,meso_taxa_index), 'pairwise_space': len(meso_protein_indexes)*len(thermo_protein_indexes), 'computed space': hits, 'time':(time2-time0)/60}

if __name__ == '__main__':
    # connect to log file and carbon tracker file
    logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')
    c_tracker = EmissionsTracker(project_name='blastp', output_dir='./logs/')
    c_tracker.start()
    time0= time.time()

    # DVC tracked parameters
    with open("./pipeline/s1_data_processing_params.yaml", "r") as stream:
        params = yaml_load(stream)['get_protein_blast_scores']
    logger.info(f"Loaded parameters: {params}")

    # get list if pairs to consider
    # this is a Series of boolean attached to index of pairs 
    logger.debug("Loading pair labels")
    pairs = pd.read_csv('./data/taxa_pairs/pair_labels.csv', index_col=0)['is_pair']
    # get the taxa pairs associated with each
    logger.debug("Loading taxa indexes for pairs")
    pair_indexes = pd.read_csv('./data/taxa_pairs/pairwise_16s_blast.csv', usecols=[0,1])
    pairs = pair_indexes[pairs]
    logger.debug(f"found {len(pairs)} pairs")
    # make it just a list of tuple for simplicity
    pairs = [(row['thermo_index'], row['meso_index']) for _, row in pairs.iterrows()]

    # sample if needed for testing
    if params['n_sample']:
        logger.info(f"Running for only {params['n_sample']} pairs")
        pairs = pairs[:params['n_sample']]
    
    # load the protein to taxa mapping
    logger.debug("Loading protein to taxa index map")
    protein_index_map = pd.read_csv(PROTEIN_SEQ_FILE, usecols=[0], sep=';')['taxa_index']
    logger.debug(f"Found {len(protein_index_map)} for {len(protein_index_map.unique())} taxa")

    # create the headers in the output file
    logger.info("Creating output file and beginning to BLASTP...")
    file = open(OUTFILENAME, 'w')
    file.write(f"thermo_protein_index,meso_protein_index,{','.join(params['blast_metrics'])}\n")
    file.close()
    
    #startup loading time
    logger.info(f"Took {(time.time() - time0)/60} m to prepare files.")

    # we need to run blast on proteins in each pair
    # this can be done in parallel
    outs = Parallel(n_jobs=params['n_jobs'])(delayed(
        lambda pair: run_blast_and_apply_blast_metric_and_append_pairwise_one_taxa_pair(*pair, protein_index_map=protein_index_map, metrics=params['blast_metrics'])
        )(pair) for pair in pairs)
    logger.debug(outs)
    emissions = c_tracker.stop()

    # get the outputs to record metrics
    outs = pd.DataFrame(outs)
    emissions_per_pair = float(emissions/len(outs))
    time_per_pair = float(outs['time'].mean())
    percent_full_pairwise = float(outs['computed space']/outs['pairwise_space'].sum())

    # save metrics
    metrics = {
        'blastp_kg_CO2_per_pair': emissions_per_pair,
        'percent_full_pairwise_blastp': percent_full_pairwise,
        'blastp_time_per_pair': time_per_pair
    }
    with open('./data/metrics/s1.4_metrics.yaml', "w") as stream:
        yaml_dump(metrics, stream)

    # save some plots
    import matplotlib.pyplot as plt
    for i, m in enumerate(params['blast_metrics']):
        scores = pd.read_csv(OUTFILENAME, usecols=[i+2])
        scores = scores[m]

        fig, ax = plt.subplots(figsize=(5,5))
        ax.set_xlabel(m)
        scores.plot.hist(bins=50, ax=ax)
        plt.savefig(f'./data/plots/blast/protein_blast_hist_{m}.png', bbox_inches='tight', dpi=250)

    
