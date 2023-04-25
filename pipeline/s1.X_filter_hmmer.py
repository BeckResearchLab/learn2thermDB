"""Script that labels putative BLAST protein pairs using domains from HMMER.

HMMER is a tool for searching sequence databases for homologs of protein.
HMM databases include Pfam, BFD.

This script HMMER to search all proteins against the HMM database. It then
labels putative BLAST protein pairs as homologous if they share some
amount of similarity on identified domains.

Eg. if a protein is labeled with 4 domains, and three are shared with
the other protein which has 5 domains, then the two proteins have 
Jaccard similarity of 3/9 = 0.33. Some threshold is set for this
to be validated as a protein pair.

Inputs to script:
-----------------
1 For DB version 1.0 using refseq and bacdive: proteins in ./data/taxa/proteins
  in the form of csv files of prien sequences and associated taxa ids.
1 For DB version 1.1 using uniprot: proteins in ./data/proteins in the form
  parquet files
2 parquet files at ./data/taxa_pairs/protein_alignment of the form "taxa_pair_XX-YY.parquet"
  Where XX and YY are the taxa ids of the two taxa in the pair. Contents are protein IDs and
  BLAST metics.

Outputs of script:
------------------
1 parquet files at ./data/taxa_pairs/hmmer_val of the form "taxa_pair_XX-YY.parquet"
  that are 1:1 mirrors of the blast files, with domains IDs identified for each protein,
  and jaccard score

Expected steps, approximately:
------------------------------
1. Load the protein sequences and associated taxonomic information from the input CSV files
  for DB version 1.0, or from the input Parquet files for DB version 1.1.
2. Use HMMER to search all proteins against the HMM database, which includes Pfam and BFD.
3. Identify the domains for each protein using HMMER output, and calculate the Jaccard similarity
  between the domains of each pair of proteins.
    In accomplishing this step, you may need to save the HMMER domains to file, and then load and 
    compute the Jaccard similarity in a separate step.
    Otherwise, do both at runtime 
4. Save the results, including protein IDs in parquet files at 
  ./data/taxa_pairs/hmmer_val of the form "taxa_pair_XX-YY.parquet", which are 1:1
  mirrors of the blast files. Each parquet file should contain the domains IDs identified for each protein
  and the calculated Jaccard score.

Make sure to handle errors and exceptions appropriately and document any assumptions or limitations of the script.
"""

import os
import pandas as pd
from pyhmmer import easel, plan7
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from joblib import Parallel, delayed

def load_protein_data(db_version, input_dir):
    protein_data = []

    if db_version == "1.0":
        for file in os.listdir(input_dir):
            if file.endswith(".csv"):
                filepath = os.path.join(input_dir, file)
                df = pd.read_csv(filepath)
                protein_data.append(df)

    elif db_version == "1.1":
        for file in os.listdir(input_dir):
            if file.endswith(".parquet"):
                filepath = os.path.join(input_dir, file)
                df = pd.read_parquet(filepath)
                protein_data.append(df)

    return pd.concat(protein_data)

def run_hmmer_on_all_proteins(proteins, hmmer_db, chunk_size):
    hmmer_results = {}

    for i in range(0, len(proteins), chunk_size):
        chunk = proteins[i:i + chunk_size]
        seqs = [SeqRecord(Seq(seq), id=pid) for pid, seq in zip(chunk["protein_id"], chunk["sequence"])]

        with easel.SequenceFile.create(b"tmp.fasta", "fasta") as seq_file:
            seq_file.write(seqs)

        with plan7.Builder(hmmer_db) as builder:
            hmmer_output = builder.search(seq_file.name)

        for hit in hmmer_output:
            domains = [dom.hmm_name.decode("utf-8") for dom in hit.domains]
            hmmer_results[hit.query_name.decode("utf-8")] = domains

    return hmmer_results

def jaccard_similarity(set1, set2):
    intersection = set1.intersection(set2)
    union = set1.union(set2)
    return len(intersection) / len(union)

def process_taxa_pair_with_precomputed_results(taxa_pair, input_dir, output_dir, hmmer_results):
    input_path = os.path.join(input_dir, taxa_pair)
    df = pd.read_parquet(input_path)

    df["domains_1"] = df["protein_id_1"].apply(lambda x: set(hmmer_results.get(x, [])))
    df["domains_2"] = df["protein_id_2"].apply(lambda x: set(hmmer_results.get(x, [])))
    df["jaccard_similarity"] = df.apply(lambda row: jaccard_similarity(row["domains_1"], row["domains_2"]), axis=1)

    output_path = os.path.join(output_dir, taxa_pair)
    df.to_parquet(output_path)

def validate_protein_pairs(db_version, protein_input_dir, taxa_pair_input_dir, output_dir, hmmer_db, n_jobs, chunk_size):
    # Load protein data
    proteins = load_protein_data(db_version, protein_input_dir)

    # Run HMMER on all proteins and store the results in a dictionary (protein_id: domain_list)
    hmmer_results = run_hmmer_on_all_proteins(proteins, hmmer_db, chunk_size)

    # Process taxa pairs using precomputed HMMER results
    taxa_pairs = os.listdir(taxa_pair_input_dir)
    Parallel(n_jobs=n_jobs)(delayed(process_taxa_pair_with_precomputed_results)(
        taxa_pair, taxa_pair_input_dir, output_dir, hmmer_results) for taxa_pair in taxa_pairs)

if __name__ == "__main__":
    # start the logger

    # load any DVC paramters, we can work to parametrize the function later, for now you can hardcode
    # parameters as global variables if you like.

    # Set up the necessary arguments and paths for your specific use case
    db_version = "1.0"  # Use "1.0" for refseq and bacdive, "1.1" for uniprot
    protein_input_dir = "./data/taxa/proteins"  # Change this path to your protein input directory
    taxa_pair_input_dir = "./data/taxa_pairs/protein_alignment"  # Change this path to your taxa pairs input directory
    output_dir = "./data/taxa_pairs/hmmer_val"  # Change this path to your desired output directory
    hmmer_db = "path/to/your/hmmer/database"  # Change this to the path of your HMMER database (Pfam, BFD, etc.)
    n_jobs = -1  # Use all available CPUs for parallel processing
    chunk_size = 100  # Adjust the chunk size based on your requirements and system resources

    # Call validate_protein_pairs() with the specified arguments
    validate_protein_pairs(db_version, protein_input_dir, taxa_pair_input_dir, output_dir, hmmer_db, n_jobs, chunk_size)
