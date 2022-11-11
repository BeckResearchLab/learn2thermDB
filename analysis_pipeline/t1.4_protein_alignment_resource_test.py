"""Test a few taxa pair alignments using current parameters.

This script looks a lot like the true s.14 script, however it takes a small sample
of pairs to run and it does not handle restarting.

Measure the carbon and time requirements.
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

from codecarbon import EmissionsTracker
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

def worker_function(alignment_handler):
    """Run one taxa pair on a worker."""
    t0=time.time()
    worker_logfile = f'./logs/t1.4_protein_alignment_resource_test_workers/pair_{alignment_handler.pair_indexes}.log'
    logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, worker_logfile, LOGLEVEL, filemode='a', worker=True)
    learn2therm.blast.logger = learn2therm.utils.start_logger_if_necessary('learn2therm.blast', worker_logfile, LOGLEVEL, filemode='a', worker=True)
    learn2therm.io.logger = learn2therm.utils.start_logger_if_necessary('learn2therm.io', worker_logfile, LOGLEVEL, filemode='a', worker=True)
    # logger.info(logging_tree.tree())
    logger.info(f"recieved pair {alignment_handler.pair_indexes}")
    
    with EmissionsTracker(project_name=f"t1.4_{alignment_handler.pair_indexes}", output_dir='./logs/t1.4_protein_alignment_resource_test_workers/') as tracker:
        out_dic = alignment_handler.run()
    t1=time.time()
    logger.info(f"Completed pair {alignment_handler.pair_indexes}. Took {(t1-t0)/60}m")
    return out_dic

if __name__ == '__main__':
    # clear logs
    shutil.rmtree('./logs/t1.4_protein_alignment_resource_test_workers/', ignore_errors=True, onerror=None)
    os.mkdir('./logs/t1.4_protein_alignment_resource_test_workers/')
    # load params
    with open("./params.yaml", "r") as stream:
        params = yaml_load(stream)['get_protein_blast_scores']
    logger = learn2therm.utils.start_logger_if_necessary("", LOGFILE, logging.INFO, filemode='w')
    logger.info(f"Loaded parameters: {params}")
    
    # get list if pairs to consider
    # this is a Series of boolean attached to index of pairs 
    logger.debug("Loading pair labels")
    pairs = pd.read_csv('./data/taxa_pairs/pair_labels.csv', index_col=0)['is_pair']
    # get the taxa pairs associated with each
    logger.debug("Loading taxa indexes for pairs")
    pair_indexes = pd.read_csv('./data/taxa_pairs/pairwise_16s_blast.csv', usecols=[0,1])
    pairs = pair_indexes[pairs]
    
    # sample a small number for this test
    pairs = pairs.sample(n=NUM_SAMPLE, random_state=777)
    logger.debug(f"Using {len(pairs)} pairs for resource test")
    # make it just a list of tuple for simplicity
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

    # create place to deposit output. This will not be cleaned up but is clearly marked as
    # temporary for the user
    output_dir = tempfile.TemporaryDirectory(dir='./tmp/', prefix='t1.4_resource_test')
    output_dir = output_dir.name

    aligners = [Aligner(
        meso_index=mi,
        thermo_index=ti,
        max_seq_len=params['max_protein_length'],
        protein_deposit=PROTEIN_SEQ_DIR,
        alignment_score_deposit=output_dir,
        metrics=params['blast_metrics'],
        alignment_params=aligner_params,
        restart=True
    ) for (ti, mi) in pairs]

    # start a cluster object and run it
    t0 = time.time()
    Cluster = getattr(dask_jobqueue, params['dask_cluster_class'])
    cluster = Cluster(silence_logs=None)
    cluster.adapt(minimum=params['n_jobs'], maximum=params['n_jobs'], target_duration='45m')

    logger.info(f"{cluster.job_script()}")
    
    with distributed.Client(cluster) as client:
        while len(client.scheduler_info()["workers"]) < params['n_jobs']:
            time.sleep(3)
        logger.info(f"Workers ready to go... submitting tasks.")
        # run the job
        results = []
        for future in distributed.as_completed(client.map(
            worker_function, aligners
        )):
            results.append(future.result())

    t1 = time.time()

    # compute metrics
    results = pd.DataFrame(results)
    metrics = {}
    metrics['apx_minutes_per_pair_avg'] = float(results['execution_time'].mean())
    metrics['apx_minutes_per_pair_std'] = float(results['execution_time'].std())
    metrics['apx_perc_protein_pairwise'] = float((results['hits']/results['pw_space']).mean())
    metrics['apx_hits_per_pair_avg'] = float(results['hits'].mean())
    metrics['apx_hits_per_pair_std'] = float(results['hits'].std())
    
    # get the carbon cost
    co2 = pd.read_csv('./logs/t1.4_protein_alignment_resource_test_workers/emissions.csv')['emissions']
    metrics['apx_co2_per_pair_avg'] = float(co2.mean())
    metrics['apx_co2_per_pair_std'] = float(co2.std())
    with open('./data/metrics/t1.4_metrics.yaml', "w") as stream:
        yaml_dump(metrics, stream)


