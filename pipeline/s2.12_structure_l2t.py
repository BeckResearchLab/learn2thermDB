"""Run structural alignment on Our data

"""
import os
from yaml import safe_load as yaml_load
from yaml import dump as yaml_dump

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('whitegrid')
sns.set_context('paper')

from typing import List
import logging

import learn2therm.utils
import learn2therm.structure

if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
    LOGLEVEL = getattr(logging, LOGLEVEL)
else:
    LOGLEVEL = logging.INFO
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'
# get the logger in subprocesses
logger = learn2therm.utils.start_logger_if_necessary('', LOGFILE, LOGLEVEL, filemode='w')

if 'FATCAT_EXEC' not in os.environ:
    raise ValueError('Please set the FATCAT_EXEC environment variable to the path of the FATCAT executable.')
else:
    FATCAT_EXEC = os.environ['FATCAT_EXEC']

STRUCTURE_DIR = './tmp/structural_alignment/'

if __name__ == "__main__":

    if not os.path.exists(STRUCTURE_DIR):
        os.makedirs(STRUCTURE_DIR, exist_ok=True)

    # Read the dataframe
    df = pd.read_csv('./data/validation/structure/sample_l2t_data.csv', index_col=0)
    # Download structures from PDB or AlphaFold2
    learn2therm.structure.download_structures(df, pdb_column='meso_pdb', u_column='meso_pid', pdb_dir=STRUCTURE_DIR)
    learn2therm.structure.download_structures(df, pdb_column='thermo_pdb', u_column='thermo_pid', pdb_dir=STRUCTURE_DIR)
    # Run FATCAT on the structures
    logger.info(f"{len(df)} initial protein pairs")
    df_result =  learn2therm.structure.run_fatcat(df, pdb_dir=STRUCTURE_DIR, exec=FATCAT_EXEC)
    df_result = df_result.dropna(subset=['p_value'])
    logger.info(f"{len(df_result)} protein pairs with both structures found and aligned")

    # Save the results
    df_result.to_csv('./data/validation/structure/l2t_sample_fatcat.csv', index=False)