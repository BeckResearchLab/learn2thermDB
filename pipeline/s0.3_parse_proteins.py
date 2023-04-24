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
    """Extract NCBI taxid and proteomes from record of uniprot"""
    id_ = None
    alphafold = None
    pdb = None
    proteomes = []
    for db_ref in record.dbxrefs:
        if db_ref.startswith("NCBI Taxonomy"):
            id_ = int(db_ref.split(':')[1])
        elif db_ref.startswith("AlphaFoldDB"):
            alphafold = str(db_ref.split(':')[1])
        elif db_ref.startswith("PDB:"):
            pdb = str(db_ref.split(':')[1])
        elif db_ref.startswith("Proteomes:"):
            proteomes.append(str(db_ref.split(':')[1]))
    return id_, alphafold, pdb, proteomes

def get_best_proteome_from_group(dataframe):
    """Find best proteome among a group of proteomes with the same species and or strain
    
    Select best proteome with priority:
    1. qualifier as 'Reference and representative proteome'
    2. qualifier as reference 'Reference proteome'
    3. qualifier as reference 'Representative proteome', most proteins
    4. remove 'Redundant proteome', and 'Excluded', whichever is left and has the most proteins is returned
    """
    if len(dataframe) == 1:
        return dataframe
    else:
        # get any 'Reference and representative proteome'
        ref_rep = dataframe[dataframe['qualifier'] == 'Reference and representative proteome']
        if len(ref_rep) == 1:
            return ref_rep
        elif len(ref_rep) > 1:
            return pd.DataFrame(ref_rep.iloc[np.argmax(ref_rep['num_proteins'])]).T 
        
        # now check for Reference proteome
        ref = dataframe[dataframe['qualifier'] == 'Reference proteome']
        if len(ref) == 1:
            return ref
        elif len(ref) > 1:
            return pd.DataFrame(ref.iloc[np.argmax(ref['num_proteins'])]).T
        
        # now check for Representative proteome
        rep = dataframe[dataframe['qualifier'] == 'Representative proteome']
        # get the one with the most proteins
        if len(rep) >= 1:
            return pd.DataFrame(rep.iloc[np.argmax(rep['num_proteins'])]).T
        
        # otherwise, remove 'Redundant proteome', and 'Excluded', whichever is left and has the most proteins is returned
        non_redun = dataframe[dataframe['qualifier'] != 'Redundant proteome']
        non_redun = non_redun[non_redun['qualifier'] != 'Excluded']
        if len(non_redun) == 1:
            return non_redun
        elif len(non_redun) > 1:
            return pd.DataFrame(non_redun.iloc[np.argmax(non_redun['num_proteins'])]).T

def get_best_strain_species_map(proteome_metadata):
    """Get a map of strain to species taxid from proteome metadata"""
    strain_species_map = dict(zip(proteome_metadata['strain_taxid'].values, proteome_metadata['species_taxid'].values))

    species_groups = proteome_metadata.groupby('species_taxid')
    species_aggregated_dfs = []
    for i, group in species_groups:
        best_proteome = get_best_proteome_from_group(group)
        species_aggregated_dfs.append(best_proteome)
    species_aggregated_df = pd.concat(species_aggregated_dfs)
    return strain_species_map, species_aggregated_df


def uniprot_to_parquet_chunking(
    source_directory: str,
    endpoint_directory: str,
    ncbi_id_filter: list,
    proteome_metadata: pd.DataFrame,
    max_filesize: int=100000,
    one_file: bool = False
):
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
    # get proteome metadata
    strain_species_map, species_aggregated_df = get_best_strain_species_map(proteome_metadata)
    species_aggregated_df.set_index('species_taxid', inplace=True)
    logger.info(f'Strain species map: {strain_species_map}')
    logger.info(f'Best proteomes for species: {species_aggregated_df}')

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
            ncbi_id, alphafold_id, pdb_id, proteomes = get_db_refs_from_xml_record(record)
            total_count += 1
            scanned_count += 1

            # convert strain taxid to species such that we can match to proteomes
            if ncbi_id in strain_species_map:
                ncbi_id = strain_species_map[ncbi_id]

            # now check for reasons to skip this protein
            # if we cannot match it to a taxa, we do not want it
            if ncbi_id is None or ncbi_id not in ncbi_id_filter:
                continue

            # if there is no reference proteome, we want to keep the protein
            if ncbi_id not in species_aggregated_df.index:
                proteome = None
                pass
            else:
                # there is a reference proteome
                reference_proteome = species_aggregated_df.loc[ncbi_id]
                reference_proteome_size = reference_proteome['num_proteins']
                reference_proteome = reference_proteome['pid']

                # if the reference proteome has enough proteins, we only want data from it
                # for this organism
                if reference_proteome_size > 2000:
                    if proteomes is not None and reference_proteome in proteomes:
                        proteome = reference_proteome
                        pass
                    else:
                        continue
                # small reference proteome, we want proteins from it and ones 
                # without any proteomes, but not proteins from other proteomes
                else:
                    if proteomes is None or reference_proteome in proteomes:
                        proteome = reference_proteome
                        pass
                    else:
                        continue

            # filters passed, we want this protein
            data.append((
                record.id,
                ncbi_id,
                pdb_id,
                alphafold_id,
                proteome,
                str(record.seq)
            ))

            # check if it is time to save and reset
            if len(data) >= max_filesize:
                df = pd.DataFrame(
                    data=data,
                    columns=['pid', 'taxid', 'pdb_id', 'alphafold_id', 'proteome', 'protein_seq'])
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
    df = pd.DataFrame(data=data, columns=['pid', 'taxid', 'pdb_id', 'alphafold_id', 'proteome', 'protein_seq'])
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

    tracker = OfflineEmissionsTracker(
        project_name=f"s0.2",
        output_dir='./data/',
        country_iso_code='USA',
        region='Washington'
    ) 
    tracker.start()

    # Load proteome metadata
    proteome_metadata = pd.read_csv('./data/uniprot/proteome_metadata.csv')
    

    # get the ncbi ids we have taxa data for
    ncbi_id_filter = list(pd.read_parquet('./data/taxa.parquet', columns=['taxid'])['taxid'])
    logger.info(f"Only considering proteins from taxa with ids {ncbi_id_filter}")

    # extract information downloaded into only needed information
    total_found = uniprot_to_parquet_chunking(
        source_directory='./data/uniprot',
        endpoint_directory='./data/proteins',
        ncbi_id_filter=ncbi_id_filter,
        proteome_metadata=proteome_metadata,
        max_filesize=params['max_prot_per_file'],
        one_file=params['dev_only_one_uniprot_file']
    )
    
    logger.info(f"Finished extracting data from uniprot, found {total_found/1000000.0}m")
    
    # get some metrics from the files using duckdb
    con = duckdb.connect()

    # get some metadata about number of proteins per taxa
    protein_per_taxa_counts = con.execute("SELECT taxid, COUNT(*) FROM './data/proteins/*.parquet' GROUP BY taxid").df()
    protein_per_taxa_counts.to_csv('./data/metrics/s0.3_protein_per_data_distr.csv')
    # how many have structures
    total_with_structures = con.execute("SELECT COUNT(*) FROM './data/proteins/*.parquet' WHERE pdb_id NOT NULL OR alphafold_id NOT NULL").fetchone()[0]

    # save metrics

    co2 = tracker.stop()
    metrics = {'s0.3_co2': float(co2)}
    metrics['n_proteins'] = int(total_found)
    metrics['percent_prot_w_struc'] = float(total_with_structures/total_found)

    with open('./data/metrics/s0.3_metrics.yaml', "w") as stream:
        yaml_dump(metrics, stream)



