"""Runs blast and reports metrics for every pair of thermophiles, mesophiles
"""
import logging
import os

from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import pandas as pd
from yaml import dump as yaml_dump
from yaml import safe_load as yaml_load
from codecarbon import OfflineEmissionsTracker
import duckdb as ddb

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

OUTDIR = "./data/taxa_pairs/"

def apply_blast_metric_make_files_chunks(
    records,
    metrics: List[str],
    chunksize: int=100000
):
    """Compute metrics for a single blast record and save to file with locking.
    
    Metrics returns dataframe of (query index, subject_index, metric value) for each metric
    """
    logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')
    dataframes = []
    current_size = 0
    total_files = 0
    for record in records:
        # compute all metrics
        metric_handler = learn2therm.blast.BlastMetrics(blast_record)
        metric_values_df = metric_handler.compute_metrics(metrics)
        logger.debug(f"Adding {len(out)} alignments for thermophile {metric_handler.qid}")
        dataframes.append(metric_values_df)
        current_size += len(metric_values_df)

        if current_size > chunksize:
            df = pd.concat(dataframes, ignore_index=True)
            df.to_parquet(OUTDIR+'taxa_pair_blast_chunk_{total_files}.parquet')
            logger.info(f"Saved file {total_files} with {current_size} pairwise comparisons")
            total_files += 1
            current_size = 0
            dataframes = []


if __name__ == '__main__':
    # connect to log file
    logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')

    # DVC tracked parameters
    with open("./params.yaml", "r") as stream:
        params = yaml_load(stream)['get_16s_blast_scores']
    logger.info(f"Loaded parameters: {params}")

    tracker = OfflineEmissionsTracker(
        project_name=f"s1.1",
        output_dir='./data/',
        country_iso_code='USA',
        region='Washington'
    ) 
    tracker.start()

    # load the thermo meso labels
    labels = pd.read_csv('./data/taxa_thermophile_labels.parquet').set_index('taxid', drop=True)['thermophile_label']

    # get the taxa indexes of themophiles and mesophiles
    thermo_indexes = list(labels[labels == True].index)
    meso_indexes = list(labels[labels == False].index)
    logger.info(f"Found taxa indexes for {thermo_indexes} thermophiles and {meso_indexes} mesophiles")

    # run blast
    # ###########
    # get the iterators of 16s sequence for mesophiles and thermophiles
    taxa = pd.read_parquet('./data/taxa.parquet', columns=['taxid', '16s_seq']).set_index('taxid', drop=True)
    thermophile_iterator = taxa.loc[thermo_indexes]
    mesophile_iterator = taxa.loc[meso_indexes]
    # need tuple of id, seq
    thermophile_iterator = thermophile_iterator.reset_index().intertuples(index=False)
    mesophile_iterator = mesophile_iterator.reset_index().intertuples(index=False)

    # sample if we neeed to
    if params['dev_n_sample'] != None:
        thermophile_iterator = [next(thermophile_iterator) for i in range(params['dev_n_sample'])]
        mesophile_iterator = [next(mesophile_iterator) for i in range(params['dev_n_sample'])]
        logger.info(f"Downsample to each of {params['dev_n_sample']} thermo and meso")
    
    # this context manager creates temporary files for BLAST inputs and outputs, again to avoid OOM
    with learn2therm.blast.BlastFiles(thermophile_iterator, mesophile_iterator, dbtype='nucl') as (query_fname, subject_fname, out_fname):
        logger.info('Running blast...')
        NcbiblastnCommandline(
            query=query_fname,
            db=subject_fname,
            outfmt=5,
            out=out_fname,
            max_target_seqs=10000000, # very large so we do not throw out any pairs. will have to increase if there is more than this num of mesos
            perc_identity=0.5, # low cutoff to lose no data
            evalue=1.0, # very large so we do not lose any hits
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

        # blast record io
        xml_f = open(out_fname, 'r')
        blast_result_records = NCBIXML.parse(xml_f)

        # compute metrics and save
        apply_blast_metric_make_files_chunks(blast_result_records, metrics=params['blast_metrics'], chunksize=100000)

    # save some plots
    con = ddb.connect()

    import matplotlib.pyplot as plt
    for i, m in enumerate(params['blast_metrics']):
        scores = con.execute(f"SELECT {m} FROM '{OUTDIR}*.parquet'").df()
        scores = scores[m]

        fig, ax = plt.subplots(figsize=(5,5))
        ax.set_xlabel(m)
        scores.plot.hist(bins=50, ax=ax)
        plt.savefig(f'./data/plots/blast/blast_16s_hist_{m}.png', bbox_inches='tight', dpi=250)

    # save metrics
    carbon = tracker.stop()
    total_hits = con.execute(f"SELECT COUNT(*) FROM '{OUTDIR}*.parquet'").fetchone()[0]
    metrics = {'percent_full_pairwise_16s_blast': float(total_hits/(len(thermo_indexes)*len(meso_indexes)))}
    metrics['s1.1_carbon'] = float(carbon)
    with open('./data/metrics/s1.1_metrics.yaml', "w") as stream:
        yaml_dump(metrics, stream)


