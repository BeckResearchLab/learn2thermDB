"""Assess the total search space for protein pairs after the taxa labeling
decision and runs a resource test.

Load the taxa information and the number of proteins for each taxa. Then
compute the pairwise search space possible as the difference for the
meso vs thermo space chosen.

Also run a resource test to estimate the cost of searching for all possible
protein pairs in pure alignment runtime not including worker overhead.
Note that for the test, random proteins are aligned instead of meso vs thermo 
proteins, which may not be representative of the actual search space.

Inputs
------
- data/taxa_thermophile_labels.parquet : contains taxaid and label for each taxa
- data/metrics/s0.3_protein_per_data_distr.csv : contains the number of proteins
    for each taxa
- data/proteins

The script looks at the quantity of proteins per mesophile and thermophile,
and extracts that to the search space available before taxa pair filering

Outputs
-------
- data/plots/protein_per_taxa_hist.png : histogram of the number of proteins per taxa
  colored by the taxa label as mesophile or thermophile
- data/plots/search_space_resources.png : plot of the resources used to compute
  protein pairs, extrapolated to the total search space
- data/metrics/t1.0_protein_search_space.yaml : contains search space size
"""
import os
import logging

import pandas as pd
import duckdb as ddb
from codecarbon import OfflineEmissionsTracker
import seaborn as sns
import matplotlib.pyplot as plt
from yaml import safe_load as yaml_load
from yaml import dump as yaml_dump
import tracemalloc
import tempfile

from learn2therm.blast import BlastAlignmentHandler, DiamondAlignmentHandler
import learn2therm.utils

if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
    LOGLEVEL = getattr(logging, LOGLEVEL)
else:
    LOGLEVEL = logging.INFO
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'

NUM_SAMPLE = 10000

def load_data(protein_counts):
    taxa = pd.read_parquet('./data/taxa_thermophile_labels.parquet')

    # outer join the two dataframes to get all taxa on taxid
    taxa = taxa.merge(protein_counts, on='taxid', how='outer', )
    # fill in the missing values for the proteins with 0
    taxa['n_proteins'] = taxa['n_proteins'].fillna(0)

    return taxa

def compute_search_space(taxa):
    mesophiles = taxa[taxa['thermophile_label'] == False]
    thermophiles = taxa[taxa['thermophile_label'] == True]

    protein_counts_meso = mesophiles['n_proteins'].values.reshape(-1,1)
    protein_counts_thermo = thermophiles['n_proteins'].values.reshape(1,-1)
    
    total_pairs = (protein_counts_meso * protein_counts_thermo).sum()

    return total_pairs

def main():
    logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, 'w')

    # now run a resource test
    # first get the params for running blast and diamond
    with open("./params.yaml", "r") as stream:
        params = yaml_load(stream)['get_protein_blast_scores']
    max_protein_length = params['max_protein_length']
    blast_params = params['method_blast_params']
    diamond_params = params['method_diamond_params']


    # get a random set of proteins
    with tempfile.TemporaryDirectory(dir='./') as tmpdir:
        conn = ddb.connect(tmpdir+'/proteins.db', read_only=False)
        conn.execute("CREATE TABLE proteins AS SELECT * FROM read_parquet('./data/proteins/*.parquet')")
        conn.commit()
        logger.info("Created temporary database with proteins table")

        # get some test proteins to run resource test alignment on
        # considering the max protein length
        query_proteins = conn.execute(
            f"""SELECT pid, protein_seq, LENGTH(protein_seq) AS len FROM proteins 
            WHERE len<={max_protein_length}
            ORDER BY RANDOM()
            LIMIT {NUM_SAMPLE}
            """).df()
        subject_proteins = conn.execute(
            f"""SELECT pid, protein_seq, LENGTH(protein_seq) AS len FROM proteins 
            WHERE len<={max_protein_length}
            ORDER BY RANDOM()
            LIMIT {NUM_SAMPLE}
            """).df()
        logger.info("Extracted random sample of proteins")

        # get some metadata about taxa protein join
        # protein total counts per organism was tracked in s0.3
        # but lets recompute that data considering a max protien length
        protein_counts = conn.execute(
            f"""SELECT taxid, COUNT(*) AS n_proteins
            FROM proteins
            WHERE LENGTH(protein_seq)<={max_protein_length}
            GROUP BY taxid""").df()


    # taxa metadata considering max protein length
    taxa = load_data(protein_counts)
    total_pairs = compute_search_space(taxa)
    logger.info(f"{total_pairs} possible pairings after sequence length filtering")
    scaling_factor = total_pairs / NUM_SAMPLE
    metrics = {'chosen_binary_possible_pairs_max_len': int(total_pairs)}

    sns.set(style='whitegrid')
    plt.figure(figsize=(10, 6))
    sns.histplot(data=taxa, x='n_proteins', hue='thermophile_label', bins=10)

    plt.xlabel('proteins per organism')
    plt.ylabel('Count')
    plt.savefig('./data/plots/protein_per_taxa_hist.png', dpi=300)
    plt.show()

    logger.info("Computed total search space: %s", total_pairs)


    # rename dfs to meet expected format
    query_proteins = query_proteins.rename(columns={'protein_seq': 'sequence'})
    subject_proteins = subject_proteins.rename(columns={'protein_seq': 'sequence'})

    # run blast
    # start carbon and ram tracking
    c_tracker = OfflineEmissionsTracker(
        project_name="t1.0_blast",
        output_dir='./data/',
        country_iso_code='USA',
        region='Washington'
    )
    c_tracker.start()
    tracemalloc.start()

    handler = BlastAlignmentHandler(
        query_proteins,
        subject_proteins,
        alignment_params=blast_params
    )
    _, blast_metadata = handler.run()
    blast_carbon = c_tracker.stop() * scaling_factor
    blast_ram = tracemalloc.get_traced_memory()[1] / 1e9
    blast_time = blast_metadata['execution_time']/60 * scaling_factor
    blast_hits = blast_metadata['hits'] * scaling_factor
    tracemalloc.stop()
    logger.info("Ran blast")

    # run diamond
    # start carbon and ram tracking
    c_tracker = OfflineEmissionsTracker(
        project_name="t1.0_diamond",
        output_dir='./data/',
        country_iso_code='USA',
        region='Washington'
    )
    c_tracker.start()
    tracemalloc.start()

    handler = DiamondAlignmentHandler(
        query_proteins,
        subject_proteins,
        alignment_params=diamond_params
    )
    _, diamond_metadata = handler.run()
    diamond_carbon = c_tracker.stop() * scaling_factor
    diamond_ram = tracemalloc.get_traced_memory()[1] / 1e9
    diamond_time = diamond_metadata['execution_time']/60 * scaling_factor
    diamond_hits = diamond_metadata['hits'] * scaling_factor
    logger.info("Ran diamond")
    tracemalloc.stop()

    # make plot with results
    x = ['BLASTP', 'DIAMOND']
    times = [blast_time, diamond_time]
    hits = [blast_hits, diamond_hits]
    ram = [blast_ram, diamond_ram]
    carbon = [blast_carbon, diamond_carbon]

    fig, ax = plt.subplots(1, 4, figsize=(20, 6))
    sns.barplot(x=x, y=times, ax=ax[0], label='total execution time')
    sns.barplot(x=x, y=hits, ax=ax[1], label='total hits')
    sns.barplot(x=x, y=ram, ax=ax[2], label='ram used')
    sns.barplot(x=x, y=carbon, ax=ax[3], label='carbon used')

    ax[0].set_ylabel('expected time (hr)')
    ax[1].set_ylabel('expected hits')
    ax[2].set_ylabel('expected min ram per processor (GB)')
    ax[3].set_ylabel('expected carbon (kg)')

    plt.savefig('./data/plots/search_space_resource_test.png', dpi=300)

    # save metrics to yaml
    with open('./data/metrics/t1.0_chosen_protein_search_space.yaml', 'w') as f:
        yaml_dump(metrics, f)

if __name__ == '__main__':
    main()
