"""Ingest raw taxa data over HTTP and FTP.

Sources
-------
Optimal Growth Temperature:
    Engqvist, M. K. M. Correlating Enzyme Annotations with a Large Set of Microbial Growth 
        Temperatures Reveals Metabolic Adaptations to Growth at Diverse Temperatures. 
        BMC Microbiology 2018, 18 (1), 177. https://doi.org/10.1186/s12866-018-1320-7.

    Data hosted on Zenodo: https://zenodo.org/record/1175609#.ZCxX4uzMIqs

16s rRNA sequences:
    Via NCBI Nucleotide Entrez
    Projects 33175[BioProject] OR 33317[BioProject]

"""
import time
import datetime
import pandas as pd
import numpy as np
from yaml import dump as yaml_dump
from yaml import safe_load as yaml_load

from Bio import Entrez, SeqIO, SeqRecord

import learn2therm.utils

import logging
import os

try:
    EMAIL = os.environ['ENV_EMAIL']
except KeyError:
    raise KeyError('Must set environmental variables `ENV_EMAIL`')
try:
    NCBI_API_KEY = os.environ['NCBI_API_KEY']
except KeyError:
    raise KeyError('Must set environmental variables `NCBI_API_KEY`')

if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
    LOGLEVEL = getattr(logging, LOGLEVEL)
else:
    LOGLEVEL = logging.INFO
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'

def nucleotide_entrez_iterator(ids: list, chunksize: int=5000, batch_size: int=200):
    """Creates an iterator of Entrez posts in chunks.
    
    Parameters
    ----------
    chunksize : int
        number of ids to post at a time
    batch_size : int
        number of ids to retrieve at a time
    """
    n_chunks = int(len(ids)/chunksize)+1
    id_chunks = np.array_split(ids, n_chunks)

    for chunk in id_chunks:
        request = Entrez.epost(db="nucleotide", id=','.join(list(chunk)))
        result = Entrez.read(request)
        webEnv = result["WebEnv"]
        queryKey = result["QueryKey"]
        
        # split further into batches for retrieval
        n_batches = int(len(chunk)/batch_size)+1
        id_batches = np.array_split(chunk, n_batches)
        check_time = time.time()

        for i, batch in enumerate(id_batches):
            if i > 0 and (time.time() - check_time < 1.0):
                time.sleep(1)
            check_time = time.time()
            handle = Entrez.efetch(
                db="nucleotide",
                id=batch,
                rettype="gb",
                retmode="text",
                retmax=batch_size+1, 
                queryKey=queryKey,
                webEnv=webEnv)

            records = SeqIO.parse(handle, 'genbank')
            for record in records:
                yield record


def get_16s_taxid(record: SeqRecord):
    """Read one record from seqio and retrieve the desired information as a dictionary."""
    sequence = str(record.seq)
    source_feature = record.features[0]
    if not source_feature.type == "source":
        raise ValueError("Expected feature from record to be 'source'")
    taxid = None
    for ref in dict(source_feature.qualifiers)['db_xref']:
        if ref.startswith('taxon'):
            taxid = int(ref.split(':')[-1])
    if taxid == None:
        raise ValueError(f"Could not interpret taxid from {dict(source_feature.qualifiers)['db_xref']}")
    return (taxid, sequence)

if __name__ == "__main__":
    # get the logger in subprocesses
    logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')
    # DVC tracked parameters
    with open("./params.yaml", "r") as stream:
        params = yaml_load(stream)['get_raw_data_taxa']
    logger.info(f"Loaded parameters: {params}")

    # start with taxa information
    # ogt data
    ogt_data=pd.read_csv('https://zenodo.org/record/1175609/files/temperature_data.tsv', sep='\t')
    # average over multiple experiments for the same taxa id
    logger.info(f"{len(ogt_data)} OGT found.")
    ogt_data = ogt_data.groupby('taxid').mean().reset_index()
    logger.info(f"OGT for {len(ogt_data)} organisms after averaging.")

    # 16s information via entrez
    Entrez.email = EMAIL
    Entrez.tool = "learn2therm"
    Entrez.api_key = NCBI_API_KEY
    # retrieve ids for the project
    count_16s = Entrez.read(Entrez.esearch(db="nucleotide", term="33175[BioProject] OR 33317[BioProject]"))['Count']
    ids_16s = Entrez.read(Entrez.esearch(db="nucleotide", term="33175[BioProject] OR 33317[BioProject]", retmax=count_16s))['IdList']
    logger.info(f"Found {len(ids_16s)} 16s sequences from NCBI")

    # post the request for the sequences
    records = nucleotide_entrez_iterator(ids_16s, chunksize=5000, batch_size=200)

    # parse the records for taxid and 16s
    data_16s = []
    for record in records:
        data_16s.append(get_16s_taxid(record))
    data_16s = pd.DataFrame(data=data_16s, columns=['taxid', '16s_seq'])
    logger.info(f"Found {len(data_16s)} 16s.")
    # drop outside of 16s sequence length
    data_16s['16s_len'] = data_16s['16s_seq'].apply(len)
    data_16s = data_16s[data_16s['16s_len'] > params['min_16s_len']]
    data_16s = data_16s[data_16s['16s_len'] < params['max_16s_len']]
    logger.info(f"{len(data_16s)} remaining 16s after dropping those shorter or longer than {params['min_16s_len']}, {params['max_16s_len']}")
    data_16s = data_16s.sort_values(by='16s_len', ascending=False)
    data_16s = data_16s.drop_duplicates(subset='taxid', keep='first')
    logger.info(f"{len(data_16s)} 16s sequences for nonduplicate taxa")

    # join the tables
    taxa_df = data_16s.merge(ogt_data, on='taxid', how='inner')
    logger.info(f"{len(taxa_df)} total taxa with overlapping data")
    taxa_df.to_parquet('./data/taxa.parquet')

    # save metrics
    metrics = {}
    metrics['n_taxa'] = len(taxa_df)
    metrics['taxa_pulled_date'] = str(datetime.datetime.now().strftime("%m/%d/%Y"))
    with open('./data/metrics/s0.0_metrics.yaml', "w") as stream:
        yaml_dump(metrics, stream)



