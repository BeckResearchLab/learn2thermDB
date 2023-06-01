"""Run and compute hmmer scores for Hait pairs."""
# system dependecies
import logging
import os
import sys
import tempfile
import time
from typing import Union
from yaml import dump as yaml_dump
from yaml import safe_load as yaml_load

# library dependencies
from codecarbon import OfflineEmissionsTracker
import duckdb as ddb
from joblib import Parallel, delayed
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import learn2therm.hmmer
import learn2therm.utils

## Paths
HMM_PATH = './data/validation/hmmer/Pfam-A.hmm'  # ./Pfam-A.hmm
DB_PATH = './data/database.ddb'

# get environmental variables
if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
    LOGLEVEL = getattr(logging, LOGLEVEL)
else:
    LOGLEVEL = logging.INFO
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'

if __name__ == '__main__':
    # load params
    with open("./params.yaml", "r") as stream:
        params = yaml_load(stream)['run_hmmer']
    # start logger/connect to log file
    logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')
    tracker = OfflineEmissionsTracker(
        project_name=f"s2.6",
        output_dir='./data/',
        country_iso_code='USA',
        region='Washington'
    ) 
    tracker.start()
    logger.info(f"Running script:{__file__}")

    # Read protein pairs data
    data_path = './data/validation/hait_pairs.csv'
    df = pd.read_csv(data_path)

    # Combine 'meso_seq' and 'thermo_seq' columns into a new dataframe
    unique_seqs_df = pd.concat([
        df[['meso_pid', 'meso_seq']].rename(columns={'meso_pid': 'pid', 'meso_seq': 'protein_seq'}),
        df[['thermo_pid', 'thermo_seq']].rename(columns={'thermo_pid': 'pid', 'thermo_seq': 'protein_seq'})
    ]).drop_duplicates()
    logger.info(f"Unique sequences data frame created, with number of proteins: {len(unique_seqs_df)}")

    # Save protein sequences to a digital sequence block
    seqblock = learn2therm.hmmer.save_to_digital_sequences(unique_seqs_df)
    logger.info("Protein sequences saved to a digital sequence block")

    # Run hmmscan
    all_hits = learn2therm.hmmer.run_pyhmmer(
        seqs=seqblock, hmms_path=HMM_PATH, prefetch=False, cpu=params['njobs'], eval_con=params['e_value'])
    logger.info("HMMER hmmscan completed")

    # Parse the output
    parsed_hits_df = learn2therm.hmmer.parse_pyhmmer(all_hits, unique_seqs_df['pid'].tolist())
    logger.info(f"{parsed_hits_df}")
    logger.info("Parsed the HMMER output")

    # Merge the parsed output with the original DataFrame
    df = df.merge(parsed_hits_df.rename(columns={'query_id': 'meso_pid', 'accession_id': 'meso_accession'}), on='meso_pid', how='left')
    df = df.merge(parsed_hits_df.rename(columns={'query_id': 'thermo_pid', 'accession_id': 'thermo_accession'}), on='thermo_pid', how='left')
    
    logger.info("Parsed hits merged with the original data")
    assert df['meso_accession'].isna().sum() == 0
    assert df['thermo_accession'].isna().sum() == 0

    # Apply evaluation function to each row
    jaccard_threshold = params['jaccard_threshold']
    df[['jaccard_score', 'functional']] = df.apply(
        lambda row: pd.Series(learn2therm.hmmer.evaluation_function(row, jaccard_threshold)),
        axis=1
    )
    df.drop(columns=['meso_seq', 'thermo_seq'], inplace=True)
    
    logger.info("Applied the evaluation function to each row")

    # Save the final dataframe
    df.to_csv('./data/validation/hmmer/hait_scores.csv', index=False)
    logger.info(f"{df[['jaccard_score', 'functional']].describe()}")

    fig, ax = plt.subplots()
    sns.histplot(data=df, x='jaccard_score', ax=ax)
    ax.set_xlabel('Jaccard score')
    logger.info("Final data saved")
    fig.savefig('./data/validation/hmmer/hait_jaccard.png', dpi=300, bbox_inches='tight')

    # Stop emissions tracker and save the output
    emissions = float(tracker.stop())

    # metrics:
    fraction_found = 1 - float(df['functional'].isna().sum()/len(df))
    min_jaccard = float(df['jaccard_score'].min())
    mean_jaccard = float(df['jaccard_score'].mean())
    fraction_functional = float(df['functional'].sum()/len(df))
    metrics = {
        'fraction_found': fraction_found,
        'min_jaccard': min_jaccard,
        'mean_jaccard': mean_jaccard,
        'fraction_functional': fraction_functional,
        'co2': emissions
    }
    with open('./data/validation/hmmer/s2.6_metrics.yaml', 'w') as f:
        yaml_dump(metrics, f)

    logger.info("Emissions tracker stopped and data saved")