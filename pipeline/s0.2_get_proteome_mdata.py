"""Download uniprot proteomes and parse them to find representative ones.

To minimize redundant proteins down the pipeline, proteins for organisms with
representative proteomes should only come from those proteomes.

This script downloads uniprot proteomes, then parses then to find representative
proteomes for each organism that exhibits proteomes.

Later, if a protein is not part of any proteome or is part of the representative
proteome for its organism, it will be kept.

Outputs
-------
- ./data/uniprot/proteome_metadata.csv
"""
import numpy as np
import pandas as pd
from yaml import dump as yaml_dump

import learn2therm.utils
import learn2therm.io

import logging
import tempfile
import os
import re
import gzip
import json
import requests
from requests.adapters import HTTPAdapter, Retry


if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
    LOGLEVEL = getattr(logging, LOGLEVEL)
else:
    LOGLEVEL = logging.INFO
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'
# get the logger in subprocesses
logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')

# functions for uniprot api
def get_next_link(headers):
        if "Link" in headers:
            match = re_next_link.match(headers["Link"])
            if match:
                return match.group(1)

def get_batch(batch_url):
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(response.headers)

# functions for parsing uniprot proteomes
def do_one_proteome(proteome):
    try:
        out = {}
        out['pid'] = proteome['id']
        try:
            lowest_lineage = proteome['taxonLineage'][-1]
            if lowest_lineage['rank'] == 'species':
                out['species_taxid'] = lowest_lineage['taxonId']
            else:
                raise KeyError()
        except KeyError:
            out['species_taxid'] = proteome['taxonomy']['taxonId']
        
        out['strain_taxid'] = proteome['taxonomy']['taxonId']
        out['qualifier'] = proteome['proteomeType']
        out['completeness'] = proteome['proteomeCompletenessReport']['cpdReport']['status']
        out['num_proteins'] = sum([c['proteinCount'] for c in proteome['components']])
    except:
        raise
    return out

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

if __name__ == '__main__':
    # downlaod uniprot proteomes
    # following tutorial at https://www.uniprot.org/help/api_queries
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))

    url = 'https://rest.uniprot.org/proteomes/search?compressed=false&format=json&query=%28%2A%29&size=500'
    proteomes = []
    for batch, total in get_batch(url):
        content = json.loads(batch.text)['results']
        proteomes.extend(content)
    logger.info("Completed downlaoding proteomes...")

    # parse into dataframe
    proteome_data = pd.DataFrame([do_one_proteome(p) for p in proteomes])
    logger.info(f"Parsed json to dataframe, shape: {proteome_data.shape}")

    # first aggregate on strains
    strain_groups = proteome_data.groupby('strain_taxid')
    strain_aggregated_dfs = []
    for i, group in strain_groups:
        best_proteome = get_best_proteome_from_group(group)
        strain_aggregated_dfs.append(best_proteome)
    strain_aggregated_df = pd.concat(strain_aggregated_dfs)
    logger.info(f"Aggregated over strains, remaining proteomes: {len(strain_aggregated_df)}")

    # now aggregate on species
    species_groups = strain_aggregated_df.groupby('species_taxid')
    species_aggregated_dfs = []
    for i, group in species_groups:
        best_proteome = get_best_proteome_from_group(group)
        species_aggregated_dfs.append(best_proteome)
    species_aggregated_df = pd.concat(species_aggregated_dfs)
    logger.info(f"Aggregated over species, remaining proteomes: {len(species_aggregated_df)}")
    species_aggregated_df.reset_index(drop=True).to_csv('./data/uniprot/proteome_metadata.csv')