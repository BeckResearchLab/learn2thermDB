"""Compare ESM distance of observed protein pairs to random pairs.

Actual pairs should be statistically significantly different than random pairs.
"""
import os
import json
from yaml import safe_load as yaml_load
from yaml import dump as yaml_dump
import pathlib
import tempfile
import requests
import time

import learn2therm.blast
import learn2therm.utils
import numpy as np
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt
import seaborn as sns
import torch
import torch.nn
from esm import Alphabet, FastaBatchedDataset, pretrained, MSATransformer
import duckdb as ddb
sns.set_style('whitegrid')
sns.set_context('paper')

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from typing import List
import logging

if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
    LOGLEVEL = getattr(logging, LOGLEVEL)
else:
    LOGLEVEL = logging.INFO
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'
# get the logger in subprocesses
logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')

LOAD_FILE = True
HAIT_URL = 'https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fprot.25866&file=prot25866-sup-0001-DataS1.xlsx'

def process_fasta_file(
    model_location: str,
    fasta_file: pathlib.Path,
    toks_per_batch: int = 4096,
    repr_layers: List[int] = [-1],
    include: List[str] = ["mean", "per_tok", "bos", "contacts"],
    truncation_seq_length: int = 300,
    nogpu: bool = False,
):
    """Run ESM on a fasta file and save the results.
    
    Adapted from https://github.com/facebookresearch/esm/blob/main/scripts/extract.py
    """
    model, alphabet = pretrained.load_model_and_alphabet(model_location)
    model.eval()
    if isinstance(model, MSATransformer):
        raise ValueError(
            "This script currently does not handle models with MSA input (MSA Transformer)."
        )
    if torch.cuda.is_available() and not nogpu:
        print('CUDA')
        model = model.cuda()
        logger.info("Transferred model to GPU")
    else:
        print('NOCUDA')
    dataset = FastaBatchedDataset.from_file(fasta_file)
    batches = dataset.get_batch_indices(toks_per_batch, extra_toks_per_seq=1)
    data_loader = torch.utils.data.DataLoader(
        dataset, collate_fn=alphabet.get_batch_converter(truncation_seq_length), batch_sampler=batches
    )
    logger.info(f"Read {fasta_file} with {len(dataset)} sequences")

    return_contacts = "contacts" in include

    assert all(-(model.num_layers + 1) <= i <= model.num_layers for i in repr_layers)
    repr_layers = [(i + model.num_layers + 1) % (model.num_layers + 1) for i in repr_layers]

    results = []

    with torch.no_grad():
        for batch_idx, (labels, strs, toks) in enumerate(data_loader):
            logger.info(
                f"Processing {batch_idx + 1} of {len(batches)} batches ({toks.size(0)} sequences)"
            )
            if torch.cuda.is_available() and not nogpu:
                toks = toks.to(device="cuda", non_blocking=True)

            out = model(toks, repr_layers=repr_layers, return_contacts=return_contacts)

            logits = out["logits"].to(device="cpu")
            representations = {
                layer: t.to(device="cpu") for layer, t in out["representations"].items()
            }
            if return_contacts:
                contacts = out["contacts"].to(device="cpu")

            for i, label in enumerate(labels):
                result = {"label": label}
                truncate_len = min(truncation_seq_length, len(strs[i]))
                if "per_tok" in include:
                    result["representations"] = {
                        layer: t[i, 1 : truncate_len + 1].clone()
                        for layer, t in representations.items()
                    }
                if "mean" in include:
                    result["mean_representations"] = {
                        layer: t[i, 1 : truncate_len + 1].mean(0).clone()
                        for layer, t in representations.items()
                    }
                if "bos" in include:
                    result["bos_representations"] = {
                        layer: t[i, 0].clone() for layer,t in representations.items()
                    }
                if return_contacts:
                    result["contacts"] = contacts[i, : truncate_len, : truncate_len].clone()
                results.append(result)
    return results

def parse_results(results):
    """Parse results from ESM into a dataframe.
    
    Original structure is list of dict of {label, mean_representations{layer_number: tensor}}}
    """
    new_dicts = []
    for result in results:
        new_dict = {}
        new_dict['pid'] = result['label']
        new_dict['tensor'] = list(result['mean_representations'].values())[-1]
        new_dicts.append(new_dict)
    return pd.DataFrame(data=new_dicts)

def run_esm_on_df(dataframe, model='esm2_t36_3B_UR50D'):
    """Create fasta files and run ESM"""

    # need iter to make blast files
    # if there are duplicate sequences, only keep one, they will be merged later
    records = []
    ids_ = []
    for i, row in dataframe.iterrows():
        if row['pid'] in ids_:
            continue
        else:
            record = SeqRecord(Seq(row['protein_seq']), id=row['pid'], description='')
            records.append(record)
            ids_.append(row['pid'])

    
    # write fasta file
    fasta_file = tempfile.NamedTemporaryFile(mode='w', dir='./tmp/', delete=False)
    SeqIO.write(records, fasta_file, 'fasta')
    fasta_file.close()
    fp = pathlib.Path(fasta_file.name)
    logger.info(f"Created fasta file {fasta_file}")

    # run ESM
    result = process_fasta_file(
        model_location=model,
        fasta_file=fp,
        toks_per_batch=4096,
        repr_layers=[-1],
        include=["mean"],
        truncation_seq_length=300,
        nogpu=False,
    )
    result = parse_results(result)

    result = result.merge(dataframe, on='pid', how='inner')
    assert len(result) == len(dataframe)

    logger.info(f"Ran ESM on {fp}")
    return result

cos = torch.nn.CosineSimilarity(dim=1, eps=1e-6)

if __name__ == '__main__':

    # load params
    with open("./params.yaml", "r") as stream:
        params = yaml_load(stream)

    temps = np.linspace(params['label_taxa']['ogt_threshold'], 75.0, 3)
    logger.info(f"Got temperatures to use for random thermo proteins: {temps}")

    if LOAD_FILE:
        random_distances = torch.load('./tmp/random_distances.pt')
        actual_distances = torch.load('./tmp/actual_distances.pt')
        hait_distances = torch.load('./tmp/hait_distances.pt')
        logger.info('Loaded precomputed distances')
    else:
        # get hait data
        hait_pairs = pd.read_csv('./data/validation/hait_pairs.csv')

        # get ESM distances
        df_ = hait_pairs[['meso_seq', 'meso_pid']].rename(columns={'meso_seq':'protein_seq', 'meso_pid':'pid'})
        hait_meso_esm_result = run_esm_on_df(df_)
        df_ = hait_pairs[['thermo_seq', 'thermo_pid']].rename(columns={'thermo_seq':'protein_seq', 'thermo_pid':'pid'})
        hait_thermo_esm_result = run_esm_on_df(df_)
        # get distances
        meso_tensor = torch.vstack(tuple(hait_meso_esm_result['tensor'].values))
        thermo_tensor = torch.vstack(tuple(hait_thermo_esm_result['tensor'].values))
        assert len(meso_tensor) == len(thermo_tensor)
        hait_distances = cos(meso_tensor, thermo_tensor)
        logger.info(f"Mean, std of ESM distance for Hait pairs: {distances.mean()}, {distances.std()}")
        torch.save(distances, './tmp/hait_distances.pt')

        # connect to database and get protein pairs, as well as random pairs
        con = ddb.connect('./data/database.ddb', read_only=True)
        random_thermo_proteins_dict = {}
        for T in temps:
            random_thermo_proteins_dict[T] = con.execute(f"""
                SELECT pid, protein_seq, temperature FROM proteins
                INNER JOIN taxa ON proteins.taxid = taxa.taxid
                WHERE temperature>{T}
                ORDER BY RANDOM()
                LIMIT 100
            """).df()
        random_meso_proteins = con.execute(f"""
            SELECT pid, protein_seq, temperature FROM proteins
            INNER JOIN taxa ON proteins.taxid = taxa.taxid
            WHERE temperature<{params['label_taxa']['ogt_threshold']}
            ORDER BY RANDOM()
            LIMIT 100
        """).df()
        logger.info(f"Got  random mesophilic proteins and sets of random themo ones based on T")
        # now get actual pairs
        actual_protein_pairs = {}
        for T in temps:
            data = con.execute(f"""
                SELECT 
                    meso_proteins.pid AS meso_pid,
                    thermo_proteins.pid AS thermo_pid, 
                    meso_proteins.protein_seq AS meso_seq,
                    thermo_proteins.protein_seq AS thermo_seq,
                    meso_taxa.temperature AS meso_temp,
                    thermo_taxa.temperature AS thermo_temp
                FROM pairs
                INNER JOIN taxa AS thermo_taxa ON pairs.thermo_taxid = thermo_taxa.taxid
                INNER JOIN taxa AS meso_taxa ON pairs.meso_taxid = meso_taxa.taxid
                INNER JOIN proteins AS thermo_proteins ON pairs.thermo_pid = thermo_proteins.pid
                INNER JOIN proteins AS meso_proteins ON pairs.meso_pid = meso_proteins.pid
                WHERE thermo_taxa.temperature>={T}
                ORDER BY RANDOM()
                LIMIT 100
            """).df()
            actual_protein_pairs[T] = data
        con.close()
        logger.info(f"Got actual meso and thermo proteins based on T")

        # run esm and check outputs
        random_meso_results = run_esm_on_df(random_meso_proteins)
        random_distances = {}
        for T, random_thermo_proteins in random_thermo_proteins_dict.items():
            # run thermo proteins in esm
            random_thermo_results = run_esm_on_df(random_thermo_proteins)
            meso_tensor = torch.vstack(tuple(random_meso_results['tensor'].values))
            thermo_tensor = torch.vstack(tuple(random_thermo_results['tensor'].values))
            assert len(meso_tensor) == len(thermo_tensor)
            print(meso_tensor.shape)

            # compute distance between meso and thermo random proteins
            distances = cos(meso_tensor, thermo_tensor)
            random_distances[T] = distances
            logger.info(f"Mean, std of ESM distance for random pairs with thermohpile T>{T}: {distances.mean()}, {distances.std()}")
        torch.save(random_distances, './tmp/random_distances.pt')

        # repeat the process, but use actual protein pairs
        actual_distances = {}
        for T, data in actual_protein_pairs.items():
            # split the data into meso and thermo stuff
            meso_sequences = data[['meso_pid', 'meso_seq']].rename(columns={'meso_pid': 'pid', 'meso_seq': 'protein_seq'})
            thermo_sequences = data[['thermo_pid', 'thermo_seq']].rename(columns={'thermo_pid': 'pid', 'thermo_seq': 'protein_seq'})

            # run thermo proteins in esm
            meso_results = run_esm_on_df(meso_sequences)
            thermo_results = run_esm_on_df(thermo_sequences)
            meso_tensor = torch.vstack(tuple(meso_results['tensor'].values))
            thermo_tensor = torch.vstack(tuple(thermo_results['tensor'].values))
            assert len(meso_tensor) == len(thermo_tensor)

            # compute distance between meso and thermo random proteins
            distances = cos(meso_tensor, thermo_tensor)
            logger.info(f"Mean, std of ESM distance for actual pairs with thermohpile T>{T}: {distances.mean()}, {distances.std()}")
            actual_distances[T] = distances
        torch.save(actual_distances, './tmp/actual_distances.pt')

    # make plot of all distributions
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    for T, distances in random_distances.items():
        break
    sns.kdeplot(data=random_distances[T], ax=ax, label=f'Random pairs, T thermo>{T}')
    sns.kdeplot(data=actual_distances[T], ax=ax, label=f'Actual pairs, T thermo>{T}')
    sns.kdeplot(data=hait_distances, ax=ax, label=f'Hait et al. pairs')
    ax.set_xlabel('ESM distance')
    ax.set_ylabel('Observed probability')
    ax.legend()
    plt.savefig('./data/plots/esm_distance_hist.png', dpi=300, bbox_inches='tight')
    
    print(actual_distances.items)

    # do some t tests
    t_tests_to_hait = {}
    for T, random_distance in random_distances.items():
        for T2, actual_distance in actual_distances.items():
            if np.isclose(T, T2):
                break
        _, p_actual = scipy.stats.ttest_ind(hait_distances, actual_distance, alternative='less')
        _, p_random = scipy.stats.ttest_ind(hait_distances, random_distance, alternative='less')
        t_tests_to_hait[float(T)] = {'our_pairs_p_value': float(p_actual), 'random_pairs_p_value': float(p_random)}
    
    # save metrics
    with open('./data/metrics/s2.X_esm.yaml', 'w') as f:
        yaml_dump(t_tests_to_hait, f)



    
