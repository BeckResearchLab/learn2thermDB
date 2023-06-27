# DVC tracked files schema

```mermaid
flowchart TD
        node1["../../data/database.ddb"]
        node2["../../data/metrics/s0.3_metrics.yaml"]
        node3["../../data/metrics/s0.3_protein_per_data_distr.csv"]
        node4["../../data/metrics/s1.0_metrics.yaml"]
        node5["../../data/metrics/s1.1_metrics.yaml"]
        node6["../../data/metrics/s1.2_metrics.yaml"]
        node7["../../data/metrics/s1.3_metrics.yaml"]
        node8["../../data/metrics/t0.2_total_search_space.csv"]
        node9["../../data/metrics/t1.0_chosen_protein_search_space.yaml"]
        node10["../../data/metrics/t1.2_metrics.yaml"]
        node11["../../data/plots/cumsum_proteins.png"]
        node12["../../data/plots/ogt_hist.png"]
        node13["../../data/plots/protein_per_taxa_hist.png"]
        node14["../../data/plots/search_space_resource_test.png"]
        node15["../../data/plots/total_search_space.png"]
        node16["../../data/protein_pairs"]
        node17["../../data/proteins"]
        node18["../../data/taxa.parquet"]
        node19["../../data/taxa_pairs/alignment"]
        node20["../../data/taxa_pairs/pair_labels"]
        node21["../../data/taxa_thermophile_labels.parquet"]
        node22["../../data/uniprot/proteome_metadata.csv"]
        node23["../../data/uniprot/uniprot_pulled_timestamp"]
        node24["../../data/validation/hait_aligned_scores.csv"]
        node25["../../data/validation/hait_alignment"]
        node26["../../data/validation/hait_pairs.csv"]
        node27["../../data/validation/hmmer/Pfam-A.hmm"]
        node28["../../data/validation/hmmer/compare_jaccard_hist.png"]
        node29["../../data/validation/hmmer/hait_jaccard.png"]
        node30["../../data/validation/hmmer/hait_n_domains.png"]
        node31["../../data/validation/hmmer/hait_scores.csv"]
        node32["../../data/validation/hmmer/hmmer_labels"]
        node33["../../data/validation/hmmer/hmmer_outputs"]
        node34["../../data/validation/hmmer/s2.6_metrics.yaml"]
        node35["../../data/validation/hmmer/s2.7_metrics.yaml"]
        node36["../../data/validation/hmmer/s2.8_metrics.yaml"]
        node37["../../data/validation/hmmer/s2.9_metrics.yaml"]
        node38["../../data/validation/structure/hait_fatcat.csv"]
        node39["../../data/validation/structure/l2t_sample_fatcat.csv"]
        node40["../../data/validation/structure/sample_l2t_data.csv"]
        node41["../../data/validation/tm/metrics.yaml"]
        node42["../../data/validation/tm/ogt_vs_tm.csv"]
        node43["../../data/validation/tm/ogt_vs_tm_check.png"]
        node1-->node25
        node1-->node28
        node1-->node33
        node1-->node35
        node1-->node37
        node1-->node40
        node1-->node41
        node1-->node42
        node1-->node43
        node3-->node5
        node3-->node8
        node3-->node9
        node3-->node11
        node3-->node13
        node3-->node14
        node3-->node15
        node3-->node19
        node16-->node1
        node16-->node32
        node16-->node36
        node17-->node1
        node17-->node7
        node17-->node9
        node17-->node10
        node17-->node13
        node17-->node14
        node17-->node16
        node18-->node1
        node18-->node2
        node18-->node3
        node18-->node4
        node18-->node5
        node18-->node8
        node18-->node11
        node18-->node12
        node18-->node15
        node18-->node17
        node18-->node19
        node18-->node21
        node19-->node1
        node19-->node6
        node19-->node7
        node19-->node10
        node19-->node16
        node19-->node20
        node20-->node1
        node20-->node7
        node20-->node10
        node20-->node16
        node21-->node1
        node21-->node5
        node21-->node9
        node21-->node13
        node21-->node14
        node21-->node19
        node22-->node2
        node22-->node3
        node22-->node17
        node23-->node2
        node23-->node3
        node23-->node17
        node24-->node25
        node26-->node24
        node26-->node29
        node26-->node30
        node26-->node31
        node26-->node34
        node26-->node38
        node27-->node29
        node27-->node30
        node27-->node31
        node27-->node33
        node27-->node34
        node27-->node35
        node32-->node28
        node32-->node37
        node33-->node32
        node33-->node36
        node40-->node39
        node44["../../data/metrics/s0.0_metrics.yaml"]
        node45["../../data/metrics/s0.1_metrics.yaml"]
        node46["../../data/validation/hmmer/s2.5_metrics.yaml"]
```
## Metrics and plots:
- `data/metrics/s0.0_metrics.yaml`: Number of organisms with both OGT and 16s, date of taxa download
- `data/plots/ogt_hist.png`: histogram of taxa OGTs
- `data/metrics/s0.1_metrics.yaml`: Date of uniprot download
- `data/metrics/s0.3_metrics.yaml`: Number of proteins, fraction with PDB or AlphaFold, carbon cost of script
- `data/metrics/s0.3_protein_per_data_distr.csv`, `data/plots/protein_per_taxa_hist.png`: table and associated plot mapping number of uniprot proteins to species taxid
  - "taxid": ncbi species id
  - "count_star()": number of proteins per organism
- `data/metrics/t0.2_total_search_space.csv`, `data/plots/total_search_space.png`: data and associated plot for the total possible number of protein pairs to search for as thermophile/mesophile temperature and window changes
- `data/plots/cumsum_proteins.png`: plot of total number of proteins as proteins per organism increases
- `data/metrics/s1.0_metrics.yaml`: Number of thermophiles and mesophiles for chosen OGT binary threshold
- `data/metrics/t1.0_chosen_protein_search_space.yaml`: total possible protein pairs among all thermophilic vs all mesophilic proteins
- `data/plots/search_space_resource_test.png`: compares BLAST and DIAMOND in therms of time, hits, carbon, and RAM if we were to search for protein pairs among all thermophilic vs all mesophilic proteins
- `data/metrics/s1.1_metrics.yaml`: Fraction of total pairwise taxa alignment that were hit, and carbon cost of 16s alignment
- `data/metrics/s1.2_metrics.yaml`: number of taxa pairs identified from 16s alignment and fraction if the total possible that is
- `data/metrics/t1.2_metrics.yaml`:emissions, total pairs, time requirement, and fraction of total possible that we would get back from protein alignment based on a resource test
- `data/metrics/s1.3_metrics.yaml`: emissions, total pairs, time requirement, and fraction of total possible that we got back from protein alignment
- `data/plots/total_search_space.png`: Explores how various aspects like search space, and amount of data remain as delta OGT increases
- `data/validation/hmmer/s2.5_metrics.yaml`: date of pulling Pfam
- `data/validation/hmmer/s2.6_metrics.yaml`: results of running hmmer on Hait, how many were annotaed by Pfam, mean jaccard score, etc.
- `data/validation/hmmer/s2.7_metrics.yaml`: number of learn2therm proteins labeled with Pfam compared to total number
- `data/validation/hmmer/s2.8_metrics.yaml`: Fraction of learn2therm proteins labeled with Pfam compared to total number
- `data/validation/hmmer/s2.9_metrics.yaml`: Summary of camprison of hmmer labeling between learn2therm and Hait
- `data/validation/hmmer/hait_jaccard.png`: Distribution of jaccard scores for Hait
- `data/validation/hmmer/compare_jaccard_hist.png`: Distribution of jaccard scores for Hait and learn2therm
- `data/validation/hmmer/hait_n_domains.png`: Distribution of number of domains identified for Hait
- `data/validation/tm/metrics.yaml`: Counts of Tm data founf, correlation between OGT and TM
- `data/validation/tm/ogt_vs_tm_check.png`: Parity plot of OGT vs TM

## Data files:
- `data/taxa.parquet`: Table of organisms with taxid, 16s sequence and OGT
  - "taxid": ncbi taxonomy id of species
  - "16s_seq": nucleotide sequence of 16s rRNA
  - "16s_len": length of 16s sequence
  - "temperature": OGT label
  - "superkingdom", "phylum", "class", "order", "family", "genus": taxonomy information
- `data/uniprot/uniprot_pulled_timestamp`: In lieu of DVC tracking the uniprot zipped files (>1000 GB), simply tracks when uniprot was pulled
- `data/uniprot/proteome_metadata.csv`: Data table of proteomes in uniprot, ID, number of proteins, quality label (eg. representative, redundant, excluded, etc.). See uniprot docs for some fields: https://www.uniprot.org/help/proteome
  - "pid": uniprot proteome id
  - "species_taxid": ncbi taxid of species level for proteome
  - "strain_taxid": ncbi taxid of strian level for proteome
  - "qualifier": eg. 'redundant' or 'representative'
  - "completeness": statistical measure of how complete the proteome is compared to other proteomes of the strain
  - "num_proteins": count of proteins in proteome
- `data/proteins`: Contains many parquet files of proteins in chunks.
  - "pid": UPKB protein id
  - "taxid": ncbi species id of protein
  - "pdb_id": PDB databse identifier if present
  - "alphafold_id": AlphaFoldDB 
  - "proteome": UP proteome id
  - "protein_seq": amino acid sequences of proteins
- `data/taxa_thermophile_labels.parquet`: bolean labels of if thermophile or mesophile
  - "taxid": NCBI species id
  - "thermophile_label": boolean, True if thermophile
- `data/taxa_pairs/alignment`: Contains parquet files of taxa 16s alignment in chunks.
  - "query_id": ncbi taxid of thermophile
  - "subject_id": ncbi taxid of mesophile
  - "X": where X is some metric extracted from the alignments
- `data/taxa_pairs/pair_labels`: Contains parquet files of organism pair labels, indexes map 1:1 to `data/taxa_pairs/alignment`
  - "is_pair": boolean if metric thresholds were met and organism pair should be considered a pair for protein alignment
- `data/protein_pairs`: Contains many parquet files of the form `align_taxa_XX-YY.parquet` where XX is thermophile taxid and YY is mesophile taxid
  - "thermo_pid": UPKB protein ID for the thermophilic protein
  - "meso_pid": UPKB protein ID for the mesophilic protein
  - "X": where X is metrics of the protein alignment
  - "thermo_taxid": NCBI taxid of the organism of the thermophilic protein
  - "meso_taxid": NCBI taxid of the organism of the mesophilic protein
- `data/database.ddb`: Relational database of taxa, proteins, taxa pairs and protein pairs.
- `data/validation/hait_pairs.csv`: Protein pairs (PDB id) from Hait et al.
  - "thermo_pdb": PDB id of thermophilic protein
  - "meso_pdb": PDB id of mesophilic protein
  - "thermo_seq": amino acid sequence of thermophilic protein
  - "meso_seq": amino acid sequence of mesophilic protein
- `data/validation/hait_aligned_scores.csv`: Local alignment scores of protein pairs from Hait et al.
  - "thermo_pdb": PDB id of thermophilic protein
  - "meso_pdb": PDB id of mesophilic protein
  - "X": where X is some metric extracted from the alignments
- `data/validation/hmmer/Pfam-A.full.hmm`: HMM database of Pfam-A
- `data/validation/hmmer/hmmer_output/*.parquet`: Pfam labels of all proteins occurring in pairs.
  - "pid": UPKB protein ID
  - "accession": Pfam labels separated by `;`
- `data/validation/hmmer/hmmer_labels/*.parquet`: Scores of protein pairs
  - "thermo_pid": UPKB protein ID for the thermophilic protein
  - "meso_pid": UPKB protein ID for the mesophilic protein
  - "score": jaccard score of accession labels between proteins
- `data/validation/hmmer/hait_scores.csv`: Scores of protein pairs from Hait et al.
  - "thermo_pdb": PDB id of thermophilic protein
  - "meso_pdb": PDB id of mesophilic protein
  - "score": jaccard score of accession labels between proteins
- `data/validation/structure/hait_fatcat.csv`: Structural alignment scores of protein pairs from Hait et al.
  - "thermo_pdb": PDB id of thermophilic protein
  - "meso_pdb": PDB id of mesophilic protein
  - "p-value": FATCAT score of structural alignment between proteins
- `data/validation/structure/sample_l2t_data.csv`: Sample of protein pairs from our data, uniformly over some pairing metrics X
  - "thermo_pid": UPKB protein ID for the thermophilic protein
  - "meso_pid": UPKB protein ID for the mesophilic protein
  - "X": where X is some metric extracted from the alignments
  - "meso_sequence": amino acid sequence of mesophilic protein
  - "thermo_sequence": amino acid sequence of thermophilic protein
- `data/validation/structure/sample_l2t_fatcat.csv`: Structural alignment scores of protein pairs from our data, uniformly over some pairing metrics X
  - "thermo_pid": UPKB protein ID for the thermophilic protein
  - "meso_pid": UPKB protein ID for the mesophilic protein
  - "X": where X is some metric extracted from the alignments
  - "p-value": FATCAT score of structural alignment between proteins
- `data/validation/tm/ogt_vs_tm.csv`: OGT vs TM of proteins within the dataset
  - "pid": UPKB protein ID
  - "temperature": OGT of organism
  - "Tm": TM of protein
  - "from": 3rd party Tm dataset source
