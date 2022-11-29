"""Run alignment on all taxa pairs

Measure the carbon cost, pairwise space, and distribution of scores
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

from codecarbon import OfflineEmissionsTracker
import dask
import dask_jobqueue
import distributed
from dvc.api import make_checkpoint

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
OUTPUT_DIR = './data/taxa_pairs/protein_alignment/'
WORKER_WAKE_UP_TIME = 25 # this is to ensure that if a worker that is about to be shut down due to previous task completetion doesn't actually start running

def try_again_read_csv(filepath: str, retries: int = 5, **kwargs):
    """Meant to circumvent weird error where pandas interacts with code_carbon"""
    i = 0
    while i < retries:
        try:
            return pd.read_csv(filepath, **kwargs)
        except pd.error.EmptyDataError:
            time.sleep(1)
            i += 1
    raise pd.error.EmptyDataError(f"Tried {retries} times to read file {filepath}")

def worker_function(alignment_handler, wakeup=None):
    """Run one taxa pair on a worker."""
    # we want to wait for execution to see if this worker is actually being used
    # or if it is in the process of being killed
    if wakeup is not None:
        time.sleep(wakeup)
    # begin execution
    t0=time.time()
    worker_logfile = f'./logs/s1.4_get_protein_blast_scores_workers/pair_{alignment_handler.pair_indexes}.log'
    logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, worker_logfile, LOGLEVEL, filemode='a', worker=True)
    learn2therm.blast.logger = learn2therm.utils.start_logger_if_necessary('learn2therm.blast', worker_logfile, LOGLEVEL, filemode='a', worker=True)
    learn2therm.io.logger = learn2therm.utils.start_logger_if_necessary('learn2therm.io', worker_logfile, LOGLEVEL, filemode='a', worker=True)

    logger.info(f"recieved pair {alignment_handler.pair_indexes}")
    
    with OfflineEmissionsTracker(
        project_name=f"s1.4_{alignment_handler.pair_indexes}",
        output_dir='./logs/',
        country_iso_code='USA',
        region='Washington'
    ) as tracker:
        out_dic = alignment_handler.run()
    t1=time.time()
    logger.info(f"Completed pair {alignment_handler.pair_indexes}. Took {(t1-t0)/60}m")
    return out_dic

if __name__ == '__main__':
    # load params
    with open("./params.yaml", "r") as stream:
        params = yaml_load(stream)['get_protein_blast_scores']
    if params['restart']:
        logger = learn2therm.utils.start_logger_if_necessary("", LOGFILE, LOGLEVEL, filemode='w')
        shutil.rmtree('./logs/s1.4_get_protein_blast_scores_workers/', ignore_errors=True, onerror=None)
        os.mkdir('./logs/s1.4_get_protein_blast_scores_workers/')
    else:
        logger = learn2therm.utils.start_logger_if_necessary("", LOGFILE, logging.INFO, filemode='a')
    logger.info(f"Loaded parameters: {params}")
    
    # get list if pairs to consider
    # this is a Series of boolean attached to index of pairs 
    logger.debug("Loading pair labels")
    pairs = pd.read_csv('./data/taxa_pairs/pair_labels.csv', index_col=0)['is_pair']
    # get the taxa pairs associated with each
    logger.debug("Loading taxa indexes for pairs")
    pair_indexes = pd.read_csv('./data/taxa_pairs/pairwise_16s_blast.csv', usecols=[0,1])
    pairs = pair_indexes[pairs]
    
    # list out the indexes to run on
    pairs = [(row['thermo_index'], row['meso_index']) for _, row in pairs.iterrows()]

    # create aligners
    alignment_method = params['method']
    try:
        aligner_params = params[f"method_{alignment_method}_params"]
    except:
        raise ValueError(f"specified method {alignment_method} for alignment but found no parameters for it.")
    
    if alignment_method == 'blast':
        Aligner = getattr(learn2therm.blast, 'BlastAlignmentHandler')
    elif alignment_method == 'diamond':
        Aligner = getattr(learn2therm.blast, 'DiamondAlignmentHandler')
    else:
        raise ValueError(f"Unknown alignment method {alignment_method}")

    # start a cluster object and run it
    t0 = time.time()
    Cluster = getattr(dask_jobqueue, params['dask_cluster_class'])
    cluster = Cluster(silence_logs=None)
    
    # determine if we have the job capacity to be killing workers after jobs
    # heuristically at three, the few is more likely that all jobs die
    # this option is to save compute when a task would normally start on a worker that just
    # finished a job
    if params['n_jobs'] > 3:
        minimum_jobs = params['n_jobs'] - 1
    else:
        minimum_jobs = 1
    cluster.adapt(minimum=minimum_jobs, maximum=params['n_jobs'], target_duration='5s')

    logger.info(f"{cluster.job_script()}")
    
    # create plugin to use for worker killing and start the client
    with distributed.Client(cluster) as client:
        if params['primary_sweep']:
            # run one without killer workers, faster option for 
            # fast tasks
            logger.info(f"Running primary fast sweep")
            with learn2therm.blast.AlignmentClusterState(
                pairs=pairs,
                client=client,
                worker_function=lambda a: worker_function(a, None),
                max_seq_len=params['max_protein_length'],
                protein_deposit=PROTEIN_SEQ_DIR,
                aligner_class=Aligner,
                alignment_score_deposit=OUTPUT_DIR,
                metrics=params['blast_metrics'],
                alignment_params=aligner_params,
                restart=params['restart'],
                killer_workers=False
            ) as futures:
                for i, future in enumerate(futures):
                    if params['checkpoint'] and (i+1) % params['checkpoint'] == 0:
                        # compute metrics
                        results = try_again_read_csv(OUTPUT_DIR+'/completion_state.metadat', index_col=0)
                        metrics = {}
                        metrics['perc_protein_pairwise'] = float((results['hits']/results['pw_space']).mean())
                        metrics['hits'] = float(results['hits'].sum())

                        with open('./data/metrics/s1.4_metrics.yaml', "w") as stream:
                            yaml_dump(metrics, stream)
                        make_checkpoint()

        # now run one with killer workers
        logger.info(f"Running safe sweep")
        with learn2therm.blast.AlignmentClusterState(
            pairs=pairs,
            client=client,
            worker_function=lambda a: worker_function(a, WORKER_WAKE_UP_TIME),
            max_seq_len=params['max_protein_length'],
            protein_deposit=PROTEIN_SEQ_DIR,
            aligner_class=Aligner,
            alignment_score_deposit=OUTPUT_DIR,
            metrics=params['blast_metrics'],
            alignment_params=aligner_params,
            restart=params['restart'] and not params['primary_sweep'],
            killer_workers=False
        ) as futures:
            for i, future in enumerate(futures):
                if params['checkpoint'] and (i+1) % params['checkpoint'] == 0:
                    # compute metrics
                    results = try_again_read_csv(OUTPUT_DIR+'/completion_state.metadat', index_col=0)
                    metrics = {}
                    metrics['perc_protein_pairwise'] = float((results['hits']/results['pw_space']).mean())
                    metrics['hits'] = float(results['hits'].sum())
                    
                    with open('./data/metrics/s1.4_metrics.yaml', "w") as stream:
                        yaml_dump(metrics, stream)
                    make_checkpoint()
            
    t1 = time.time()
    logger.info(f'Completed all tasks, took {(t1-t0)/60} min')
    
    # compute metrics
    results = try_again_read_csv(OUTPUT_DIR+'/completion_state.metadat', index_col=0)
    metrics = {}
    metrics['perc_protein_pairwise'] = float((results['hits']/results['pw_space']).mean())
    metrics['hits'] = float(results['hits'].sum())
    with open('./data/metrics/s1.4_metrics.yaml', "w") as stream:
        yaml_dump(metrics, stream)


