"""Test a few taxa pair alignments using current parameters.

This script looks a lot like the true s.14 script, however it takes a small sample
of pairs to run and it does not handle restarting.

Measure the carbon and time requirements.
"""
from yaml import safe_load as yaml_load
from yaml import dump as yaml_dump
import os
import time
import tempfile
# set the dask config
os.environ['DASK_CONFIG'] = './.config/dask/'

import pandas as pd
import logging
from logging_tree import printout

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
NUM_SAMPLE = 4

def worker_function(alignment_handler):
    """Run one taxa pair on a worker."""
    logging.setLoggerClass(learn2therm.utils.WorkerLogger)
    logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='a')
    logging.setLoggerClass(learn2therm.utils.WorkerLogger)
    fh = logger.handlers[-1]
    fh.setFormatter(logging.Formatter("%(filename)s:%(worker)s-%(asctime)s %(levelname)-8s %(message)s"))
    # logger = logging.LoggerAdapter(logger, {'worker': worker.name})
    logger.info(f"Recieved pair {alignment_handler.pair_indexes}")
    logger.info(printout())
    
    tracker = EmissionsTracker(project_name=f"t1.4", output_dir='./logs/')
    tracker.start()

    out_dic = alignment_handler.run()
    emissions = tracker.stop()
    out_dic.update({'co2': emissions})
    return out_dic

if __name__ == '__main__':
    # load params
    with open("./params.yaml", "r") as stream:
        params = yaml_load(stream)['get_protein_blast_scores']
    
    logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, logging.INFO, filemode='w')
    logger.info(f"Loaded parameters: {params}")
    logger.debug(dask.config.config['jobqueue'])
    
    
    # get list if pairs to consider
    # this is a Series of boolean attached to index of pairs 
    logger.debug("Loading pair labels")
    pairs = pd.read_csv('./data/taxa_pairs/pair_labels.csv', index_col=0)['is_pair']
    # get the taxa pairs associated with each
    logger.debug("Loading taxa indexes for pairs")
    pair_indexes = pd.read_csv('./data/taxa_pairs/pairwise_16s_blast.csv', usecols=[0,1])
    pairs = pair_indexes[pairs]
    
    print(printout())
    
    # sample a small number for this test
    pairs = pairs.sample(n=NUM_SAMPLE)
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
    cluster.adapt(minimum=1, maximum=params['n_jobs'])
    client = distributed.Client(cluster)

    # run the job
    results = []
    for future in distributed.as_completed(client.map(
        worker_function, aligners, resources={'processes':1}
    )):
        results.append(future.result())
    cluster.close()
    t1 = time.time()

    # compute metrics
    results = pd.DataFrame(results)
    metrics = {}
    metrics['apx_minutes_per_pair'] = float((t1-t0)/60/NUM_SAMPLE)
    metrics['apx_co2_per_pair'] = float(results['co2'].mean())
    metrics['apx_perc_protein_pairwise'] = float((results['hits']/results['pw_space']).mean())
    with open('./data/metrics/t1.4_metrics.yaml', "w") as stream:
        yaml_dump(metrics, stream)


