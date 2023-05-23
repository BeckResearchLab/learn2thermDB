"""Compare lowest level of taxonomic match to 16s alignment score.

Inputs
------
- data/taxa_pairs/alignment : contains alignment scores for taxa pairs and their taxids
- data/taxa.parquet

Outputs
-------
- data/plots/t1.1_16s_vs_taxonomy.png : distribution of %id seperated by 
    lowest level of taxonomic match
"""
from yaml import safe_load as yaml_load
from yaml import dump as yaml_dump
import os
import shutil
import time
import tempfile
# set the dask config

import pandas as pd
import duckdb as ddb
import logging
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('whitegrid')
sns.set_context('paper')

import learn2therm.utils

if __name__ == '__main__':
    with tempfile.TemporaryDirectory(dir='./tmp', ) as tmpdir:
        conn = ddb.connect(tmpdir+'/proteins.db', read_only=False)
        conn.execute("CREATE TABLE taxa AS SELECT * FROM read_parquet('./data/taxa.parquet')")
        conn.commit()

        # now create table of taxa pairs
        conn.execute("CREATE TABLE pairs AS SELECT subject_id, query_id, scaled_local_symmetric_percent_id FROM read_parquet('./data/taxa_pairs/alignment/*.parquet')")
        conn.execute("CREATE INDEX meso_index ON pairs (subject_id)")
        conn.execute("CREATE INDEX thermo_index ON pairs (query_id)")
        conn.commit()

        # join the taxonomy over the scores
        data = conn.execute("""
            SELECT
                pairs.scaled_local_symmetric_percent_id AS percent_id,
                taxa_meso.superkingdom AS meso_superkingdom,
                taxa_meso.phylum AS meso_phylum,
                taxa_meso.class AS meso_class,
                taxa_meso.order AS meso_order,
                taxa_meso.family AS meso_family,
                taxa_meso.genus AS meso_genus,
                taxa_thermo.superkingdom AS thermo_superkingdom,
                taxa_thermo.phylum AS thermo_phylum,
                taxa_thermo.class AS thermo_class,
                taxa_thermo.order AS thermo_order,
                taxa_thermo.family AS thermo_family,
                taxa_thermo.genus AS thermo_genus
            FROM pairs 
            INNER JOIN taxa AS taxa_meso ON (pairs.subject_id = taxa_meso.taxid)
            INNER JOIN taxa AS taxa_thermo ON (pairs.query_id = taxa_thermo.taxid)
        """).df() 
    
    # determine highest rank of taxonomic match
    def do_one(row):
        lowest_match = 'None'
        for level in ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus']:
            if row[f'meso_{level}'] == row[f'thermo_{level}']:
                lowest_match = level
            else:
                break
        return lowest_match
    data['match'] = data.apply(do_one, axis=1)

    # plot results
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.boxplot(data=data, x='percent_id', y='match', order =['superkingdom', 'phylum', 'class', 'order', 'family', 'genus'], ax=ax)
    ax.set_xlabel('Percent id')
    ax.set_ylabel('Lowest match group')
    plt.savefig('./data/plots/t1.1_16s_vs_taxonomy.png', dpi=300, bbox_inches='tight')
        


