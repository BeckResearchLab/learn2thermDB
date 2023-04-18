"""Assess the total search space for protein pairs after the taxa labeling
decision.

Load the taxa information and the number of proteins for each taxa. Then
compute the pairwise search space possible as the difference for the
meso vs thermo space chosen

Inputs
------
- data/taxa_thermophile_labels.parquet : contains taxaid and label for each taxa
- data/metrics/s0.2_protein_per_data_distr.csv : contains the number of proteins
    for each taxa

The script looks at the quantity of proteins per mesophile and thermophile,
and extracts that to the search space available before taxa pair filering

Outputs
-------
- data/plots/protein_per_taxa_hist.png : histogram of the number of proteins per taxa
  colored by the taxa label as mesophile or thermophile
- data/metrics/t1.0_protein_search_space.yaml : contains search space size
"""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import yaml

def load_data():
    taxa = pd.read_parquet('./data/taxa_thermophile_labels.parquet')
    protein_counts = pd.read_csv('./data/metrics/s0.2_protein_per_data_distr.csv', index_col=0)
    protein_counts = protein_counts.rename(columns={'count_star()': 'n_proteins'})

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
    taxa = load_data()
    total_pairs = compute_search_space(taxa)

    sns.set(style='whitegrid')
    plt.figure(figsize=(10, 6))
    sns.histplot(data=taxa, x='n_proteins', hue='thermophile_label', bins=10)

    plt.xlabel('proteins per organism')
    plt.ylabel('Count')
    plt.savefig('./data/plots/protein_per_taxa_hist.png', dpi=300)
    plt.show()

    # save metrics to yaml
    metrics = {'total_possible_pairs_no_filter': int(total_pairs)}
    with open('./data/metrics/t1.0_chosen_protein_search_space.yaml', 'w') as f:
        yaml.dump(metrics, f)

if __name__ == '__main__':
    main()
