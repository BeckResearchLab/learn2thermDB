"""Sample taxa pairs across delta GT space, run them with blast vs diamond

This script looks a lot like the true s.14 script, however it takes a small sample
of pairs.

Compares the distribution of hits between methods
"""
from yaml import safe_load as yaml_load
from yaml import dump as yaml_dump
import os
import shutil
import time
import tempfile
# set the dask config
os.environ['DASK_CONFIG'] = './.config/dask/'

import pandas as pd
import logging
import logging_tree

import dask
import dask_jobqueue
import distributed

import learn2therm.utils
import learn2therm.blast

if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
    LOGLEVEL = getattr(logging, LOGLEVEL)
else:
    LOGLEVEL = logging.INFO
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'

PROTEIN_SEQ_DIR = './data/taxa/proteins/'
NUM_SAMPLE = 32
WORKER_WAKE_UP_TIME = 25 # this is to ensure that if a worker that is about to be shut down due to previous task completetion doesn't actually start running

def worker_function(alignment_handler):
    """Run one taxa pair on a worker."""
    # we want to wait for execution to see if this worker is actually being used
    # or if it is in the process of being killed
    time.sleep(WORKER_WAKE_UP_TIME)
    # begin execution
    t0=time.time()
    worker_logfile = f'./logs/t1.5_protein_alignment_hits_workers/pair_{alignment_handler.pair_indexes}.log'
    logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, worker_logfile, LOGLEVEL, filemode='a', worker=True)
    learn2therm.blast.logger = learn2therm.utils.start_logger_if_necessary('learn2therm.blast', worker_logfile, LOGLEVEL, filemode='a', worker=True)
    learn2therm.io.logger = learn2therm.utils.start_logger_if_necessary('learn2therm.io', worker_logfile, LOGLEVEL, filemode='a', worker=True)
    # logger.info(logging_tree.tree())
    logger.info(f"recieved pair {alignment_handler.pair_indexes}")

    out_dic = alignment_handler.run()
    t1=time.time()
    logger.info(f"Completed pair {alignment_handler.pair_indexes}. Took {(t1-t0)/60}m")
    return out_dic

if __name__ == '__main__':
    # clear logs
    shutil.rmtree('./logs/t1.5_protein_alignment_hits_workers/', ignore_errors=True, onerror=None)
    os.mkdir('./logs/t1.5_protein_alignment_hits_workers/')
    # load params
    with open("./params.yaml", "r") as stream:
        params = yaml_load(stream)['get_protein_blast_scores']
    logger = learn2therm.utils.start_logger_if_necessary("", LOGFILE, LOGLEVEL, filemode='w')
    logger.info(f"Loaded parameters: {params}")
    
    # get list if pairs to consider
    # this is a Series of boolean attached to index of pairs 
    logger.debug("Loading pair labels")
    labels = pd.read_csv('./data/taxa_pairs/pair_labels.csv', index_col=0)['is_pair']
    # get the taxa pairs associated with each
    logger.debug("Loading taxa indexes for pairs")
    pair_indexes = pd.read_csv('./data/taxa_pairs/pairwise_16s_blast.csv')
    pairs = pair_indexes[labels]
    pairs[['thermo_index', 'meso_index']] = pairs[['thermo_index', 'meso_index']].astype(int)
    
    # get OGT information and create a join
    taxa_T_info = pd.read_csv('./data/taxa/labels.csv', index_col=0)
    meso_OGT = taxa_T_info.loc[pairs['meso_index']]['ogt'].values
    thermo_OGT = taxa_T_info.loc[pairs['thermo_index']]['ogt'].values
    delta_OGT = thermo_OGT - meso_OGT
    pairs['del_ogt'] = delta_OGT
    
    # build a distribution over the samples and pick from each bin
    values = pairs[pairs.columns[2:]]
    logger.debug(f'Determining density over {values.describe()}')
    cuts = pd.DataFrame([pd.cut(values[c], NUM_SAMPLE, labels=range(NUM_SAMPLE)) for c in values.columns])
    cuts = cuts.T.apply(tuple, axis=1) # combine into one column
    # this is a tuple indicating the bin for each example
    # shuffle the dataset with as seed and select the first from each bin
    pairs['cuts'] = cuts
    pairs = pairs.sample(frac=1, random_state=777).reset_index(drop=True)
    pairs.drop_duplicates(subset=['cuts'], keep='first', inplace=True)
    # take N
    pairs = pairs.sample(n=NUM_SAMPLE, random_state=777).drop(columns=['cuts'])
    # Now we have N pairs that are apprximately uniformly sampled over OGT and 16s metric space
    pairs_indexes = [(int(row['thermo_index']), int(row['meso_index'])) for _, row in pairs.iterrows()]

    # start a cluster object and run it
    t0 = time.time()
    Cluster = getattr(dask_jobqueue, params['dask_cluster_class'])
    cluster = Cluster(silence_logs=None)
    
    # determine if we have the job capacity to be killing workers after jobs
    # heuristically at three, the few is more likely that all jobs die
    # this option is to save compute when a task would normally start on a worker that just
    # finished a job
    if params['n_jobs'] > 3:
        minimum_jobs = 2
    else:
        minimum_jobs = 1
    cluster.adapt(minimum=minimum_jobs, maximum=params['n_jobs'], target_duration='5s')

    logger.info(f"{cluster.job_script()}")
    
    # create plugin to use for worker killing and start the client
    with distributed.Client(cluster) as client:
        # start with blast
        ##############
        Aligner = learn2therm.blast.BlastAlignmentHandler
        aligner_params = params["method_blast_params"]
        output_dir = tempfile.TemporaryDirectory(dir='./tmp/', prefix='t1.5_blast_dia_comparison')
        output_dir = output_dir.name
        results_blast = []
        with learn2therm.blast.AlignmentClusterFutures(
            pairs=pairs_indexes,
            client=client,
            worker_function=worker_function,
            max_seq_len=params['max_protein_length'],
            protein_deposit=PROTEIN_SEQ_DIR,
            aligner_class=Aligner,
            alignment_score_deposit=output_dir,
            metrics=params['blast_metrics'],
            alignment_params=aligner_params,
            restart=True
        ) as futures:
            for future in futures:
                results_blast.append(future)
        
        results_blast = pd.DataFrame(results_blast)
        results_blast['thermo_index'] = results_blast['pair'].apply(lambda s: s.split('-')[0]).astype(int)
        results_blast['meso_index'] = results_blast['pair'].apply(lambda s: s.split('-')[1]).astype(int)
        results_blast.set_index(['thermo_index', 'meso_index'], drop=True, inplace=True)

        # now diamond
        ##############
        Aligner = learn2therm.blast.DiamondAlignmentHandler
        aligner_params = params["method_diamond_params"]
        output_dir = tempfile.TemporaryDirectory(dir='./tmp/', prefix='t1.5_blast_dia_comparison')
        output_dir = output_dir.name
        results_dia = []
        with learn2therm.blast.AlignmentClusterFutures(
            pairs=pairs_indexes,
            client=client,
            worker_function=worker_function,
            max_seq_len=params['max_protein_length'],
            protein_deposit=PROTEIN_SEQ_DIR,
            aligner_class=Aligner,
            alignment_score_deposit=output_dir,
            metrics=params['blast_metrics'],
            alignment_params=aligner_params,
            restart=True
        ) as futures:
            for future in futures:
                results_dia.append(future)
                
        results_dia = pd.DataFrame(results_dia)
        results_dia['thermo_index'] = results_dia['pair'].apply(lambda s: s.split('-')[0]).astype(int)
        results_dia['meso_index'] = results_dia['pair'].apply(lambda s: s.split('-')[1]).astype(int)
        results_dia.set_index(['thermo_index', 'meso_index'], drop=True, inplace=True)
        
    # combine the dataframes from two methods
    results = pd.merge(results_blast, results_dia, left_index=True, right_index=True, suffixes=['_blast', '_dia'])
    pairs.set_index(['thermo_index', 'meso_index'], drop=True, inplace=True)
    results.to_csv('./tmp/results.csv')
    pairs.to_csv('./tmp/pairs.csv')
    
    # add to the raw data, eg ogt and 16s scores
    results = pd.concat([pairs, results], axis=1)
    results.to_csv('./data/analysis/blast_vs_diamond_hits.csv')
    t1 = time.time()
    logger.info(f'Completed all tasks, took {(t1-t0)/60} min')
    
    

