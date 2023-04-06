"""Ingest Uniprot

Parse raw uniprot files into parquet with minimal information
ready for downstream processing.
"""
import duckdb
import numpy as np
import pandas as pd
from yaml import dump as yaml_dump
from yaml import safe_load as yaml_load

import learn2therm.utils
import learn2therm.io

from codecarbon import OfflineEmissionsTracker

import logging
import os
import shutil


try:
    EMAIL = os.environ['ENV_EMAIL']
except KeyError:
    raise KeyError('Must set environmental variables `ENV_EMAIL`')
if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
    LOGLEVEL = getattr(logging, LOGLEVEL)
else:
    LOGLEVEL = logging.INFO
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'
# get the logger in subprocesses
logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')

def get_db_refs_from_xml_record(record):
    """Extract NCBI taxid from record of uniprot"""
    id_ = None
    alphafold = None
    for db_ref in record.dbxrefs:
        if db_ref.startswith("NCBI Taxonomy"):
            id_ = int(db_ref.split(':')[1])
        elif db_ref.startswith("AlphaFoldDB"):
            alphafold = str(db_ref.split(':')[1])
    return id_, alphafold

def uniprot_to_parquet_chunking(source_directory: str, endpoint_directory: str, ncbi_id_filter: list, max_filesize: int=100000):
    """Iteratres though downloaded uniprot files and produces fixed size parquet.
    
    Final files only contain sequence and NCBI id.

    Parameters
    ----------
    source_directory : str
        All .xml.gz files within directory are considered
    endpoint_directory: str
        Creates files of pattern 'uniprot_chunk_X.parquet'
    ncbi_id_filter : list of int
        ignore proteins that did not come from these organisms
    """
    source_files = [f for f in os.listdir(source_directory) if f.endswith('.xml.gz')]

    # start data structure
    data = []
    total_files = 0
    total_count = 0

    # each distinct file downloaded from uniprot has many sequences
    for filename in source_files:
        logger.info(f"Opening proteins in {filename}")
        records = learn2therm.io.seq_io_gnuzipped(source_directory+'/'+filename, filetype='uniprot-xml')
        for record in records:
            # check if we have this taxa, if so record, otherwise skip
            ncbi_id, alphafold_id = get_db_refs_from_xml_record(record)
            if ncbi_id is None or ncbi_id not in ncbi_id_filter:
                continue
            else:
                data.append((
                    record.id,
                    ncbi_id,
                    alphafold_id,
                    str(record.seq)
                ))

            # check if it is time to save and reset
            if len(data) >= max_filesize:
                df = pd.DataFrame(data=data, columns=['pid', 'taxid', 'alphafold_id', 'protein_seq'])
                df.to_parquet(endpoint_directory+'/'+f"uniprot_chunk_{total_files}.parquet")

                logger.info(f"File number {total_files+1} completed and saved with size {len(data)}.")
                # reset
                total_count += len(df)
                total_files += 1
                data = []
            else:
                pass
    # finish up
    df = pd.DataFrame(data=data, columns=['pid', 'taxid', 'alphafold_id', 'protein_seq'])
    df.to_parquet(endpoint_directory+'/'+f"uniprot_chunk_{total_files}.parquet")
    logger.info(f"File number {total_files+1} completed and saved with size {len(data)}.")
    total_count += len(df)
    return total_count

if __name__ == "__main__":
    # DVC tracked parameters
    with open("./params.yaml", "r") as stream:
        params = yaml_load(stream)['parse_proteins']
    logger.info(f"Loaded parameters: {params}")

    if not os.path.exists('./data/proteins'):
        os.mkdir('./data/proteins')
    else:
        shutil.rmtree('./data/proteins')
        os.mkdir('./data/proteins')

    tracker = OfflineEmissionsTracker(
        project_name=f"s0.2",
        output_dir='./data/',
        country_iso_code='USA',
        region='Washington'
    ) 
    tracker.start()

    # get the ncbi ids we have taxa data for
    ncbi_id_filter = list(pd.read_parquet('./data/taxa.parquet', columns=['taxid'])['taxid'])
    logger.info(f"Only considering proteins from taxa with ids {ncbi_id_filter}")

    # extract information downloaded into only needed information
    uniprot_to_parquet_chunking(
        source_directory='./data/uniprot',
        endpoint_directory='./data/proteins',
        ncbi_id_filter=ncbi_id_filter,
        max_filesize=params['max_prot_per_file'])
    
    logger.info(f"Finished extracting data from uniprot.")
    
    # get some metrics from the files using duckdb
    con = duckdb.connect()
    total_proteins = con.execute("SELECT COUNT(*) FROM './data/proteins/*.parquet'").fetchone()[0]

    # get some metadata about number of proteins per taxa
    protein_per_taxa_counts = con.execute("SELECT taxid, COUNT(*) FROM './data/proteins/*.parquet' GROUP BY taxid").df()
    protein_per_taxa_counts.to_csv('./data/metrics/s0.2_protein_per_data_distr.csv')

    # save metrics
    metrics = {}
    carbon = tracker.stop()
    metrics['n_proteins'] = int(total_proteins)
    metrics['s0.2_carbon'] = float(carbon)
    with open('./data/metrics/s0.2_metrics.yaml', "w") as stream:
        yaml_dump(metrics, stream)



