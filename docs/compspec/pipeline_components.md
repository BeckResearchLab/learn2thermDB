# Script components
Scripts are executed to track the project from raw data to trained model. This docuemnt details the DVC tracked execution
Each DVC stage is associated with a script

1. `s0.0_get_ncbi_refseq.sh`  
    Retrieve NCBI refseq genome files.  
    _Outputs_: `data/refseq` containing gbff files for all bacteria and archaea  
    _Metrics_: `data_extraction.number_taxa`  
2. `s0.1_get_bacdive_ogt.py`
    Retrieve optimal growth temperatures from bacdive.  
    Loops the refseq files and uses `taxid` to search bacdive.  
    _Params_: `data_extraction.n_jobs`  
    _Inputs_: `data/refseq`  
    _Outputs_: `data/taxa/taxa_info_and_ogt.csv` which maps taxid to ogt or ogt label, also contains auxilarry taxa information from NCBI for taxa  
    _Metrics_: `n_taxa, n_have_bacdive, n_have_ogt, n_have_growth_temp`   
3. `s1.0_label_taxa.py`  
    Label taxa as mesophile or thermophile.
    _Params_: `data_processing.ogt_threshold`,  `data_processing.ogt_determination_method`, define what do do with all of the data that we have temps for but not OGT  
    _Inputs_: `data/taxa/taxa_info_and_ogt.csv`   
    _Outputs_: `data/taxa/labels.csv`  
    _Metrics_: `n_meso`, `n_thermo`  
4. `s1.1_get_protein_sequences.py`  
    Extract sequences and 16s rRNA for each taxa.  
    _Inputs_: `data/taxa/taxa_info_and_ogt.csv`, `data/refseq`  
    _Outputs_: `data/taxa/16s_rRNA.csv`, `data/taxa/proteins/taxa_index_X.csv`, contains sequences and sequence labels  
    _Metrics_: `n_total_sequences`, `n_taxa_with_16srRNA`
5. `s1.2_get_16s_blast_scores.py`  
    Compute pairwise BLAST pairings of meso vs therma 16s rRNA sequences.   
    Uses `labels.csv`, grabs the 16s rRNA sequencies for each meso and therma. BLASTS all possible pairs.  
    _Params_: `data_processing.16s_blast_parameters` (there are a number), `data_processing.16s_blast_metrics`
    _Inputs_: `data/taxa/16s_rRNA.csv`, `data/taxa/labels.csv`  
    _Outputs_: `data/taxa_pairs/pairwise_16s_blast.csv` table of id and BLAST scores. Meso and thermo as rows and columns.  
    _Metrics_: `percent_full_pairwise_16s_blast` fraction of the full pairwise table of thermo vs mesophiles that we have hits for. Sources of less than 1 are: taxa without 16s, blast not returning a hit
6. `s1.3_label_all_pairs.py`  
    Create a list of taxa pairs that meet a minumum 16s rRNA BLAST score.  
    _Params_: `data_processing.thresholds` defines thresholds on 16s blast metrics to consider a pair  
    _Inputs_: `data/taxa_pairs/pairwise_16s_blast.csv`  
    _Outputs_: `data/taxa_pairs/pair_labels.csv` contains meso-thermo pairs (indexes for taxa) that met threshold and tracks blast score.  
    _Metrics_: `num_pairs_conservative` number of pairs that passed thresholds. Only this number will we blastp, `true_pair_ratio` fraction of pairs with metrics (eg n_therm\*n_meso\*percent_full_pairwise_16s_blast) that we kept based on thresholding  
7. `s1.4_get_protein_blast_scores.py`  
    For each meso-thermo taxa pair, compute pairwise BLAST scores of all protein sequences.  
    _Params_: `data_processing.protein_blast_parameters`  
    _Inputs_: `data/taxa_pairs/pair_labels.csv`, `data/taxa/proteins/`  
    _Outputs_: `data/taxa_pairs/protein_pairs_blast.csv` contains field pair_index which maps 1:1 to the index in pair_labels.csv, then protein indexes for meso and thermo protein, and their blast scores  
    _Metrics_: `blastp_kg_CO2_per_pair`, `percent_full_pairwise_blastp`: fraction of total possible protein pairs within labeled taxa pairs that we found hits, `blastp_time_per_pair`: minutes per pair to blast on average. Multiply these two per pair metrics by num_pairs_conservative to get total emissions and cpu time
8. `s1.5_spin_up_rdb.py`  
    Uses all files created so far to create an RDB of taxa, proteins, taxa pairs, and protein pairs.  
    The only constrictive choices at this point were meso-thermo OGT threshold and minumum 16s score to consider a pair.  
    _Inputs_:  
    - `data/taxa/taxa_info_and_ogt.csv`, `data/taxa/labels.csv` merge to create table of taxa
    - `data/taxa/proteins/` use index to create proteins table relational to taxa id
    - `data/taxa_pairs/pair_labels.csv` use to create table of taxa pairs, index of taxa are in columns
    - `data/taxa_pairs/protein_pairs_blast.csv` use to create table of protein pairs. Have to get tax ids from pair_labels.csv and then protein ids in the csv should match to proteins table

    _Outputs_: `data/learn2thermDB.sql`

## Notes
- Humood will likely take the lead on blasting