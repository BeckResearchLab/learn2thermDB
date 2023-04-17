"""Extract the total possible search space for the protein sequence pairs.

Load the taxa information and the number of proteins for each taxa. Then
compute the pairwise search space possible as the difference in taxa OGT
is increased.

Inputs
------
- data/taxa.parquet : contains taxaid and OGT for each taxa
- data/metrics/s0.2_protein_per_data_distr.csv : contains the number of proteins
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
    protein_counts = pd.read_csv('./data/metrics/s0.2_protein_per_data_distr.csv')
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

    return total_pairs

def main():
    taxa, protein_counts = load_data()

    thresholds = [30.0, 35.0, 40.0, 45.0, 50.0]
    window_widths = [0.0, 5.0, 10.0, 15.0]

    search_space_data = []

    for ogt_threshold in thresholds:
        for window_width in window_widths:
            total_pairs = compute_search_space(taxa, protein_counts, ogt_threshold, window_width)
            search_space_data.append([ogt_threshold, window_width, total_pairs])

    search_space_df = pd.DataFrame(search_space_data, columns=['threshold', 'window_width', 'total_pairs'])
    search_space_df.to_csv('data/metrics/t0.2_total_search_space.csv', index=False)

    sns.set(style='whitegrid')
    plt.figure(figsize=(10, 6))
    sns.lineplot(data=search_space_df, x='threshold', y='total_pairs', hue='window_width', marker='o')

    plt.xlabel('OGT Threshold (°C)')
    plt.ylabel('Total Possible Protein Pairs')
    plt.legend(title='Window Width (°C)', loc='upper left')
    plt.savefig('data/plots/total_search_space.png', dpi=300)
    plt.show()

if __name__ == '__main__':
    main()
