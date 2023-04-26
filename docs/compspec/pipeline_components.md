# Script components

Scripts are executed to track the project from raw data to trained model. This docuemnt details the DVC tracked execution
Each DVC stage is associated with a script

1. `s0.0_get_raw_data_taxa.py`  
    Pull most recent NCBI 16s r RNA sequences, and OGT records from Enqvist 
    _Params_: `min_16s_len`, `max_16s_len` number of nucleotides required to keep and organism
    _Outputs_: `data/taxa.parquet`, columns include OGT, 16s sequence, taxid
    _Metrics_: `n_taxa` total number of labeled organisms, `taxa_pulled_date` when the data was retireved
2. `s0.1_get_raw_data_proteins.py`
    Retrieve single cell uniprot.  
    Uses FTP to download very large uniprot files    
    _Inputs_: `data/taxa.parquet`  
    _Outputs_: `data/taxa/uniprot/uniprot_pulled_timestamp` indicates when files
    were pulled
    _Metrics_: `taxa_pulled_date` when data was retrieved
    _Untracked Outputs_: The script produces `*.xml.gz` files that are untracked because they take so long to download. DVC ignores them, subsequent calls to the script skip downloading files already present. 
3. `s0.2_get_proteome_mdata.py`  
    Get metadata for UniProt proteomes.
    Selects one "best" proteome per organism. 
    _Outputs_: `data/uniprot/proteome_metadata.csv`  
4. `s0.3_parse_proteins.py`  
    Extract minimal protein data and store in efficient file format.  
    _Params_: `max_prot_per_file` size of parquet files
    _Inputs_: `data/taxa/uniprot/uniprot_pulled_timestamp`, `data/taxa.parquet` 
    _Outputs_: `data/proteins`, contains proteins in chunked files of the form `*.parquet`. Columns include protein sequence, database identifiers, and associated taxa IDs
    _Metrics_: `n_proteins` total protein count, `percent_prot_w_struc` fraction of proteins with PDB or alphafold id
5. `s1.0_label_taxa.py`
    Assign booleans for taxa as thermophile
    _Params_: `ogt_threshold` binary thermophile threshold
    _Inputs_: `data/taxa.parquet`
    _Outputs_: `data/taxa_thermophile_labels.parquet`
    _Metrics_: `n_meso`, `n_thermo`
6. `s1.1_get_16s_blast_scores.py`  
    Compute pairwise BLAST pairings of meso vs therma 16s rRNA sequences.   
    Uses `labels.csv`, grabs the 16s rRNA sequencies for each meso and therma. BLASTS all possible pairs.  
    _Params_: `16s_blast_parameters` (there are a number), `16s_blast_metrics`
    _Inputs_: `data/taxa*.parquet`  
    _Outputs_: `data/taxa_pairwise_16s_blast.parquet` table of id and BLAST scores. 
    _Metrics_: `percent_full_pairwise_16s_blast` fraction of the full pairwise table of thermo vs mesophiles that we have hits for. 
7. `s1.2_label_all_pairs.py`  
    Create a list of taxa pairs that meet a minumum 16s rRNA BLAST score.  
    _Params_: `blast_metric_thresholds` defines thresholds on 16s blast metrics to consider a pair  
    _Inputs_: `data/taxa_pairwise_16s_blast.parquet`  
    _Outputs_: `data/taxa_pair_labels.parquet` contains meso-thermo pairs (indexes for taxa) and booleans for pairing
    _Metrics_: `n_taxa_pairs` number of pairs that passed thresholds. Only this number will we blastp, `true_pair_ratio` fraction of pairs with metrics eg n_taxa_pairs/(n_therm\*n_meso) that will be searched
    for protein paris

# BELOW NEEDS UPDATING

8. `s1.4_get_protein_blast_scores.py`  
    For each meso-thermo taxa pair, compute pairwise BLAST scores of all protein sequences.  
    _Params_: `data_processing.protein_blast_parameters`  
    _Inputs_: `data/taxa_pairs/pair_labels.csv`, `data/taxa/proteins/`  
    _Outputs_: `data/taxa_pairs/protein_pairs_blast.csv` contains field pair_index which maps 1:1 to the index in pair_labels.csv, then protein indexes for meso and thermo protein, and their blast scores  
    _Metrics_: `blastp_kg_CO2_per_pair`, `percent_full_pairwise_blastp`: fraction of total possible protein pairs within labeled taxa pairs that we found hits, `blastp_time_per_pair`: minutes per pair to blast on average. Multiply these two per pair metrics by num_pairs_conservative to get total emissions and cpu time
9. `s1.5_spin_up_rdb.py`  
    Uses all files created so far to create an RDB of taxa, proteins, taxa pairs, and protein pairs.  
    The only constrictive choices at this point were meso-thermo OGT threshold and minumum 16s score to consider a pair.  
    _Inputs_:  
    - `data/taxa/taxa_info_and_ogt.csv`, `data/taxa/labels.csv` merge to create table of taxa
    - `data/taxa/proteins/` use index to create proteins table relational to taxa id
    - `data/taxa_pairs/pair_labels.csv` use to create table of taxa pairs, index of taxa are in columns
    - `data/taxa_pairs/protein_pairs_blast.csv` use to create table of protein pairs. Have to get tax ids from pair_labels.csv and then protein ids in the csv should match to proteins table

   _Outputs_: `data/learn2thermDB.sql`

## Notes

-[X] Humood will likely take the lead on blasting
