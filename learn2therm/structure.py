"""
Contains functions to download structures from PDB or AlphaFold2 and run FATCAT on them.

Run on Hait pairs
"""

from Bio.PDB import PDBList
import os
import requests
import pandas as pd
import subprocess
import time
import tempfile
from functools import partial
import logging
import multiprocessing as mp
logger = logging.getLogger(__name__)


def download_structures(df, pdb_column, u_column, pdb_dir):
    """
    Download structures from PDB or AlphaFold2
    :param df: dataframe containing PDB IDs and UniProt IDs
    :param pdb_column: name of the column containing PDB IDs
    :param u_column: name of the column containing UniProt IDs
    :param pdb_dir: path to the directory where the structures will be downloaded
    :return: None
    """
    start_time = time.time()  # Start measuring time
    pdbl = PDBList()
    if not os.path.exists(pdb_dir):
        os.makedirs(pdb_dir)
        
    for i, row in df.iterrows():
        pdb_id = row[pdb_column]
        uniprot_id = row[u_column]
        if not pd.isna(pdb_id):  # check for NaN value in PDB IDs column
            # if we have the file already, skip
            if os.path.exists(os.path.join(pdb_dir, f'{pdb_id}.pdb')):
                continue

            pdbl.retrieve_pdb_file(pdb_id, pdir=pdb_dir, file_format='pdb')
            file_path = os.path.join(pdb_dir, f'pdb{pdb_id.lower()}.ent')
            if os.path.exists(file_path):
                os.rename(os.path.join(file_path), os.path.join(pdb_dir, f'{pdb_id}.pdb'))
            else:
                pass
        elif isinstance(uniprot_id, str):  # download structure using UniProt ID
            filename = f'{pdb_dir}/{uniprot_id}.pdb'
            # skip if we have
            if os.path.exists(filename):
                continue

            url = f'https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb'
            response = requests.get(url)
            if response.ok:
                with open(filename, 'wb') as f:
                    f.write(response.content)
                logger.info(f"Downloaded file for {uniprot_id} to {filename}")
            else:
                logger.info(f"Failed to download file for {uniprot_id}: {response.status_code} - {response.reason}")
        else:
            logger.info(f"No PDB ID or UniProt ID available for index {i}")
    end_time = time.time()  # Stop measuring time
    execution_time = end_time - start_time
    logger.info(f"Execution time: {execution_time} seconds")
    pass

def _do_one_fatcat(row, pdb_dir, exec='FATCAT'):
        index, row = row
        if not pd.isna(row['meso_pdb']):
            p1 = row['meso_pdb']
        else:
            p1 = row['meso_pid']
        
        if not pd.isna(row['thermo_pdb']):
            p2 = row['thermo_pdb']
        else:
            p2 = row['thermo_pid']
        
        # Check if the structure files exist in the 'checking' folder
        p1_file = f'{p1}.pdb'
        p2_file = f'{p2}.pdb'
        if not os.path.exists(os.path.join(pdb_dir, p1_file)) or not os.path.exists(os.path.join(pdb_dir, p2_file)):
            # Append the index of the row to the list of rows to be dropped
            return (index, pd.NA)

        # Set the FATCAT command and its arguments
        cmd = [exec, '-p1', p1_file, '-p2', p2_file, '-i', pdb_dir, '-q']
        
        # Run the FATCAT command and capture the output
        logger.info ("Running FATCAT for files %s and %s", p1_file, p2_file)
        process = subprocess.run(cmd, capture_output=True, text=True)
        stdout, stderr = process.stdout, process.stderr
        logger.debug(stdout)
        logger.info("Ran FATCAT for files %s and %s", p1_file, p2_file)

        # Find the line containing the p-value
        p_value_line = next(line for line in stdout.split('\n') if line.startswith("P-value"))

        # Extract the p-value
        p_value = float(p_value_line.split()[1])
        assert p_value >= 0 and p_value <= 1, f"Invalid p-value {p_value} for files {p1_file} and {p2_file}"
        logger.info("P-value for files %s and %s is %s", p1_file, p2_file, p_value)
        return (index, p_value)

def run_fatcat(df, pdb_dir, exec='FATCAT'):
    """
    Run FATCAT on the structures in the dataframe
    :param df: dataframe containing PDB IDs and UniProt IDs
    :param pdb_dir: path to the directory where the structures are downloaded
    :param exec: str, executable or location of fatcat executable
    :return: dataframe containing the p-values
    """
    p_values = {}  # List to store the extracted p-values
    rows_to_drop = []  # List to store the indices of rows to be dropped
        
    pool = mp.Pool(mp.cpu_count())
    outs = pool.map(partial(_do_one_fatcat, pdb_dir=pdb_dir, exec=exec), df.iterrows())
    p_values = dict(outs)

    # Drop the rows with missing structure files from the dataframe
    p_values = pd.Series(p_values)
    
    df['p_value'] = p_values
    return df