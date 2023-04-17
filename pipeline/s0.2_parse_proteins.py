"""Ingest Uniprot

Download raw uniprot data files, not tracked by dvc.
"""
import duckdb
import numpy as np
import pandas as pd
from yaml import dump as yaml_dump
from yaml import safe_load as yaml_load

import learn2therm.utils
import learn2therm.io

import codecarbon

import datetime
import logging
import os
import shutil
import tempfile

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
    pdb = None
    for db_ref in record.dbxrefs:
        if db_ref.startswith("NCBI Taxonomy"):
            id_ = int(db_ref.split(':')[1])
        elif db_ref.startswith("AlphaFoldDB"):
            alphafold = str(db_ref.split(':')[1])
        elif db_ref.startswith("PDB:"):
            pdb = str(db_ref.split(':')[1])
    return id_, alphafold, pdb

def uniprot_to_parquet_chunking(source_directory: str, endpoint_directory: str, ncbi_id_filter: list, max_filesize: int=100000, one_file: bool = False):
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
    if one_file:
        source_files = [f for f in source_files if 'sprot_archaea' in f]
    
    # start data structure
    data = []
    total_files = 0
    scanned_count = 0
    total_count = 0
    taken_count = 0

    # each distinct file downloaded from uniprot has many sequences
    for filename in source_files:
        logger.info(f"Opening proteins in {filename}")
        records = learn2therm.io.seq_io_gnuzipped(source_directory+'/'+filename, filetype='uniprot-xml')
        for record in records:
            logger.debug(f"Opened record {record}")
            # check if we have this taxa, if so record, otherwise skip
            ncbi_id, alphafold_id, pdb_id = get_db_refs_from_xml_record(record)
            total_count += 1
            scanned_count += 1
            # if we cannot match it to a taxa, we do not want it
            if ncbi_id is None or ncbi_id not in ncbi_id_filter:
                continue
            else:
                data.append((
                    record.id,
                    ncbi_id,
                    pdb_id,
                    alphafold_id,
                    str(record.seq)
                ))

            # check if it is time to save and reset
            if len(data) >= max_filesize:
                df = pd.DataFrame(
                    data=data,
                    columns=['pid', 'taxid', 'pdb_id', 'alphafold_id', 'protein_seq'])
                df.to_parquet(endpoint_directory+'/'+f"uniprot_chunk_{total_files}.parquet")

                logger.info(f"File number {total_files+1} complete. Total proteins scanned {total_count}.Found {len(data)} from {scanned_count} in chunk. Return ratio {len(data)/scanned_count}.")
                # reset
                taken_count += len(df)
                total_files += 1
                scanned_count = 0
                data = []
            else:
                pass
        logger.info(f"Completed parsing {filename}")
    # finish up
    df = pd.DataFrame(data=data, columns=['pid', 'taxid', 'pdb_id', 'alphafold_id', 'protein_seq'])
    df.to_parquet(endpoint_directory+'/'+f"uniprot_chunk_{total_files}.parquet")
    logger.info(f"File number {total_files+1} complete. Total proteins scanned {total_count}.Found {len(data)} from {scanned_count} in chunk. Return ratio {len(data)/scanned_count}.")
    taken_count += len(df)
    return taken_count

if __name__ == "__main__":
    # DVC tracked parameters
    with open("./params.yaml", "r") as stream:
        params = yaml_load(stream)['parse_proteins']
    logger.info(f"Loaded parameters: {params}")

    if not os.path.exists('./data/proteins'):
        os.mkdir('./data/proteins')
        
    data_tracker = codecarbon.OfflineEmissionsTracker( 
        project_name="data_prep_classifier",
        output_dir="./data/",
        country_iso_code="USA",
        region="washington"
    )
    data_tracker.start()

    # get the ncbi ids we have taxa data for
    ncbi_id_filter = list(pd.read_parquet('./data/taxa.parquet', columns=['taxid'])['taxid'])
    logger.info(f"Only considering proteins from taxa with ids {ncbi_id_filter}")

    # extract information downloaded into only needed information
    total_found = uniprot_to_parquet_chunking(
        source_directory='./data/uniprot',
        endpoint_directory='./data/proteins',
        ncbi_id_filter=ncbi_id_filter,
        max_filesize=params['max_prot_per_file'],
        one_file=params['dev_only_one_uniprot_file']
    )
    
    logger.info(f"Finished extracting data from uniprot, found {total_found/1000000.0}m")
    
    # get some metrics from the files using duckdb
    con = duckdb.connect()

    # get some metadata about number of proteins per taxa
    protein_per_taxa_counts = con.execute("SELECT taxid, COUNT(*) FROM './data/proteins/*.parquet' GROUP BY taxid").df()
    protein_per_taxa_counts.to_csv('./data/metrics/s0.2_protein_per_data_distr.csv')
    # how many have structures
    total_with_structures = con.execute("SELECT COUNT(*) FROM './data/proteins/*.parquet' WHERE pdb_id NOT NULL OR alphafold_id NOT NULL").fetchone()[0]

    # save metrics
    co2 = data_tracker.stop()
    metrics = {'s0.2_co2': float(co2)}
    metrics['n_proteins'] = int(total_found)
    metrics['percent_prot_w_struc'] = float(total_with_structures/total_found)
    metrics['protein_pulled_date'] = str(datetime.datetime.now().strftime("%m/%d/%Y"))
    with open('./data/metrics/s0.2_metrics.yaml', "w") as stream:
        yaml_dump(metrics, stream)



