"""Extract AA sequences for Hait meso-thermo pairs to be used for
downstream validation

Outputs
-------
- data/hait_pairs.csv: the pairs from Hait et al. that we will use as ground truth
    Columns: meso_pdb, thermo_pdb, meso_seq, thermo_seq
    PDB ids and amino acid sequences of the pairs
"""
from yaml import safe_load as yaml_load
from yaml import dump as yaml_dump
import os
import sys
import shutil
import time
import tempfile
import requests
# set the dask config

import pandas as pd
import duckdb as ddb
import logging
import tqdm


import learn2therm.utils

if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
    LOGLEVEL = getattr(logging, LOGLEVEL)
else:
    LOGLEVEL = logging.INFO
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'
logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')

if __name__ == '__main__':
    # get hait data
    hait_ht1_pairs = pd.read_excel(
        './external/prot25866-sup-0001-datas1.xlsx',
        sheet_name='PAIRS',
        usecols="A:B",
        skiprows=8).rename(columns={'HT':'T'}).dropna()
    hait_t1_pairs = pd.read_excel(
        './external/prot25866-sup-0001-datas1.xlsx',
        sheet_name='PAIRS',
        usecols="D:E",
        skiprows=8).dropna().rename(columns={'T':'T', 'M.1':'M'})
    hait_ht2_pairs = pd.read_excel(
        './external/prot25866-sup-0001-datas1.xlsx',
        sheet_name='PAIRS',
        usecols="H:I",
        skiprows=8).rename(columns={'HT.1':'T', 'M.2':'M'}).dropna()
    hait_t2_pairs = pd.read_excel(
        './external/prot25866-sup-0001-datas1.xlsx',
        sheet_name='PAIRS',
        usecols="K:L",
        skiprows=8).rename(columns={'T.1':'T', 'M.3':'M'}).dropna()
    hait_pairs = pd.concat([hait_ht1_pairs, hait_t1_pairs, hait_ht2_pairs, hait_t2_pairs], join='outer', axis=0, ignore_index=True)
    logger.info(f"Got HAIT data with {len(hait_pairs)} pairs")
    # query for AA sequences from UniProt
    base_q = 'https://rest.uniprot.org/uniprotkb/search?query='
    pdbs = list(hait_pairs['M'].values.reshape(-1))
    pdbs.extend(list(hait_pairs['T'].values.reshape(-1)))
    pdbs = list(set(pdbs))
    logger.info(f"{len(pdbs)} unique PDB ids in HAIT data")
    # programatic access in chunks
    n = 50
    pdb_to_seq = {}
    pdb_to_pid = {}
    complete = False
    count_done = 0
    timeout_counter = 0
    with tqdm.tqdm(total=len(pdbs)) as pbar:
        while not complete:
            not_done = [p for p in pdbs if p not in pdb_to_seq]
            if len(not_done) == 0:
                complete = True
                break
            chunk = not_done[:n]
            q = [f'(xref:pdb-{p})OR' for p in chunk]
            q = ''.join(q)[:-2]
            q = base_q+q
            print(q)
            r = requests.get(q)
            r.raise_for_status()
            results = r.json()
            results = results['results']
            # map pdb id to seq
            for result in results:
                seq = result['sequence']['value']
                xrefs = result['uniProtKBCrossReferences']
                id_ = result['primaryAccession']
                for xref in xrefs:
                    if xref['database'] == 'PDB':
                        if not xref['id'] in chunk:
                            continue
                        pdb_to_seq[xref['id']] = seq
                        pdb_to_pid[xref['id']] = id_
            if len(pdb_to_seq) == count_done:
                timeout_counter += 1
                if timeout_counter > 5:
                    logger.info(f"Timeout counter exceeded. Exiting.")
                    break
                logger.info(f"Timeout counter: {timeout_counter}")
            pbar.update(len(pdb_to_seq) - count_done)
            count_done = len(pdb_to_seq)
            logger.info(f"Number retrieved: {len(pdb_to_seq)}")
            time.sleep(3)
    # create dataframe with seqs
    # assert len(pdb_to_seq) == len(pdbs)
    # assert len(pdb_to_pid) == len(pdbs)
    hait_pairs.rename(columns={'M': 'meso_pdb', 'T': 'thermo_pdb'}, inplace=True)
    hait_pairs['meso_seq'] = hait_pairs['meso_pdb'].map(pdb_to_seq)
    hait_pairs['thermo_seq'] = hait_pairs['thermo_pdb'].map(pdb_to_seq)
    hait_pairs['meso_pid'] = hait_pairs['meso_pdb'].map(pdb_to_pid)
    hait_pairs['thermo_pid'] = hait_pairs['thermo_pdb'].map(pdb_to_pid)
    hait_pairs.dropna(inplace=True)
    hait_pairs.drop_duplicates(subset=['meso_pid', 'thermo_pid'], inplace=True)
    hait_pairs.to_csv('./data/validation/hait_pairs.csv')
    logger.info(f"Got HAIT data with {len(hait_pairs)} pairs after dropping NA and duplicates")