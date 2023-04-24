"""Extract the total possible search space for the protein sequence pairs.

Load the taxa information and the number of proteins for each taxa. Then
compute the pairwise search space possible as the difference in taxa OGT
is increased.

Inputs
------
- data/taxa.parquet : contains taxaid and OGT for each taxa
- data/metrics/s0.3_protein_per_data_distr.csv : contains the number of proteins
    for each taxa

The script loops through binary OGT thresholds from 30.0 to 50.0 C and
windows OGT boundary widths of 5.0, 10.0, and 15.0 C. For each combination
of threshold and window, the script computes the number of protein pairs
that could be searched. The script creates a plot of the search space. 

Outputs
-------
- data/plots/total_search_space.png : plot of the total search space
- data/metrics/t0.2_total_search_space.csv : table of the total search space
"""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def load_data():
    taxa = pd.read_parquet('./data/taxa.parquet')
    protein_counts = pd.read_csv('./data/metrics/s0.3_protein_per_data_distr.csv')
    protein_counts = protein_counts.rename(columns={'count_star()': 'n_proteins'})

    return taxa, protein_counts

def compute_search_space(taxa, protein_counts, ogt_threshold, window_width):
    min_ogt = ogt_threshold - window_width
    max_ogt = ogt_threshold + window_width

    taxa_below_window = taxa[taxa['temperature'] < min_ogt]
    taxa_above_window = taxa[taxa['temperature'] > max_ogt]

    protein_counts_below = protein_counts[protein_counts['taxid'].isin(taxa_below_window['taxid'])]
    protein_counts_above = protein_counts[protein_counts['taxid'].isin(taxa_above_window['taxid'])]
    protein_counts_below = protein_counts_below['n_proteins'].values.reshape(-1,1)
    protein_counts_above = protein_counts_above['n_proteins'].values.reshape(1,-1)
    
    total_pairs = (protein_counts_below * protein_counts_above).sum()
    n_meso_total = protein_counts_below.sum()
    n_thermo_total = protein_counts_above.sum()
    n_meso_mean = protein_counts_below.mean()
    n_thermo_mean = protein_counts_above.mean()

    return n_meso_total, n_thermo_total, n_meso_mean, n_thermo_mean, total_pairs

def main():
    taxa, protein_counts = load_data()

    thresholds = [30.0, 35.0, 40.0, 45.0, 50.0]
    window_widths = [0.0, 5.0, 10.0, 15.0]

    search_space_data = []

    for ogt_threshold in thresholds:
        for window_width in window_widths:
            n_meso_total, n_thermo_total, n_meso_mean, n_thermo_mean, total_pairs = compute_search_space(taxa, protein_counts, ogt_threshold, window_width)
            search_space_data.append([ogt_threshold, window_width, n_meso_total, n_thermo_total, n_meso_mean, n_thermo_mean, total_pairs])

    search_space_df = pd.DataFrame(search_space_data, columns=['threshold', 'window_width', 'n_meso_total', 'n_thermo_total', 'n_meso_mean', 'n_thermo_mean', 'total_pairs'])
    search_space_df.to_csv('data/metrics/t0.2_total_search_space.csv', index=False)

    sns.set(style='whitegrid')
    fig, ax = plt.subplots(5,1,figsize=(10, 30))

    sns.lineplot(data=search_space_df, x='threshold', y='n_meso_total', hue='window_width', marker='o', ax=ax[0])
    ax[0].set_xlabel('OGT Threshold (°C)')
    ax[0].set_ylabel('Number of mesophilic proteins')

    sns.lineplot(data=search_space_df, x='threshold', y='n_thermo_total', hue='window_width', marker='o', ax=ax[1])
    ax[1].set_xlabel('OGT Threshold (°C)')
    ax[1].set_ylabel('Number of thermophilic proteins')

    sns.lineplot(data=search_space_df, x='threshold', y='total_pairs', hue='window_width', marker='o', ax=ax[2])
    ax[2].set_xlabel('OGT Threshold (°C)')
    ax[2].set_ylabel('Total possible pairs')

    sns.lineplot(data=search_space_df, x='threshold', y='n_meso_mean', hue='window_width', marker='o', ax=ax[3])
    ax[3].set_xlabel('OGT Threshold (°C)')
    ax[3].set_ylabel('Meso proteins per organism')

    sns.lineplot(data=search_space_df, x='threshold', y='n_thermo_mean', hue='window_width', marker='o', ax=ax[4])
    ax[4].set_xlabel('OGT Threshold (°C)')
    ax[4].set_ylabel('Thermo proteins per organism')

    plt.legend(title='Window Width (°C)', loc='upper left')
    plt.savefig('data/plots/total_search_space.png', dpi=300)
    plt.show()

    # now create a plot of the cumulative sum of proteins as number of proteins
    # per taxa increases
    protein_counts_sorted = protein_counts.sort_values('n_proteins', ascending=True)['n_proteins']
    protein_counts_cumsum = protein_counts_sorted.cumsum()

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(protein_counts_sorted, protein_counts_cumsum)
    ax.set_xlabel('Number of Proteins in organism')
    ax.set_ylabel('Cumulative Sum of Proteins')
    plt.savefig('data/plots/cumsum_proteins.png', dpi=300)

if __name__ == '__main__':
    main()
