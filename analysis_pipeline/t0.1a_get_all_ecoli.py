"""Extract all information from Uniprot for all E. coli proteins.

The persuit is to find a way to remove gene redundancy. We do not
want a dataset biased towards organisms with redudant genes like model
organisms. We want to remove redundancy as much as possible.

Loop through uniprot XML files, and make a new one with only E. coli

Inputs
------
- data/uniprot/uniprot_*.xml.gz : uniprot files

Outputs
-------
- data/uniprot/uniprot_ecoli.xml : uniprot file with only E. coli
"""

import duckdb
import numpy as np
import pandas as pd
from yaml import dump as yaml_dump
from yaml import safe_load as yaml_load

import learn2therm.utils
import learn2therm.io
from Bio import SeqIO

import datetime
from codecarbon import OfflineEmissionsTracker
import logging
import os
import shutil

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
    for db_ref in record.dbxrefs:
        if db_ref.startswith("NCBI Taxonomy"):
            id_ = int(db_ref.split(':')[1])
            break
    return id_

def uniprot_to_parquet_chunking(source_directory: str):
    """Iteratres though downloaded uniprot files and produces fixed size parquet.
    
    Final files only contain sequence and NCBI id.

    Parameters
    ----------
    source_directory : str
        All bacteria.xml.gz files within directory are considered
    """
    source_files = [f for f in os.listdir(source_directory) if f.endswith('bacteria.xml.gz')]
    
    # start data structure
    data = []
    scanned_count = 0
    taken_count = 0

    # each distinct file downloaded from uniprot has many sequences
    for filename in source_files:
        logger.info(f"Opening proteins in {filename}")
        records = learn2therm.io.seq_io_gnuzipped(source_directory+'/'+filename, filetype='uniprot-xml')
        for record in records:
            logger.debug(f"Opened record {record}")
            # check if we have this taxa, if so record, otherwise skip
            ncbi_id = get_db_refs_from_xml_record(record)
            scanned_count += 1
            # if we cannot match it to a taxa, we do not want it
            if ncbi_id is None or ncbi_id != 562:
                continue
            else:
                data.append(record)
                taken_count += 1
            if taken_count == 100:
                break
        logger.info(f"Completed parsing {filename}")
    logger.info(f"Finished parsing all files, scanned {scanned_count} proteins, of which {taken_count} were E. coli")
    return data

if __name__ == "__main__":
    # DVC tracked parameters
    with open("./params.yaml", "r") as stream:
        params = yaml_load(stream)['parse_proteins']
    logger.info(f"Loaded parameters: {params}")

    tracker = OfflineEmissionsTracker(
        project_name=f"t0.1a",
        output_dir='./data/',
        country_iso_code='USA',
        region='Washington'
    ) 
    tracker.start()

    # extract information downloaded into only needed information
    records = uniprot_to_parquet_chunking(
        source_directory='./data/uniprot',
    )
    
    # save records to xml with biopython
    SeqIO.write(data, "./data/uniprot/uniprot_ecoli.xml", "uniprot-xml")
    logger.info(f"Saved E. coli proteins to ./data/uniprot/uniprot_ecoli.xml")

    co2 = tracker.stop()
    metrics = {'t0.1a_co2': float(co2)}

    with open('./data/metrics/t0.1a_metrics.yaml', "w") as stream:
        yaml_dump(metrics, stream)



