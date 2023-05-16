"""Compare ESM distance of observed protein pairs to random pairs.

Actual pairs should be statistically significantly different than random pairs.
"""
import os
from yaml import safe_load as yaml_load
import pathlib
import tempfile
import torch
from esm import Alphabet, FastaBatchedDataset, pretrained, MSATransformer
import duckdb as ddb

import learn2therm.blast
import learn2therm.utils
import numpy as np
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt
import seaborn as sns
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
    records = []
    for i, row in dataframe.iterrows():
        record = SeqRecord(Seq(row['protein_seq']), id=row['pid'], description='')
        records.append(record)
    
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

if __name__ == '__main__':

    # load params
    with open("./params.yaml", "r") as stream:
        params = yaml_load(stream)

    temps = np.linspace(params['label_taxa']['ogt_threshold'], 75.0, 3)
    logger.info(f"Got temperatures to use for random thermo proteins: {temps}")

    # connect to database and get protein pairs, as well as random pairs
    con = ddb.connect('./tmp/database', read_only=True)
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
        ids = con.execute(f"""
            SELECT meso_pid, thermo_pid FROM pairs
            INNER JOIN taxa AS thermo_taxa ON pairs.thermo_taxid = thermo_taxa.taxid
            WHERE thermo_taxa.temperature>={T}
        """).df()
        meso_sequences = con.execute(f"""
            SELECT pid, protein_seq, temperature FROM proteins
            INNER JOIN taxa ON proteins.taxid = taxa.taxid
            WHERE pid IN {tuple(ids['meso_pid'])}
        """).df()
        thermo_sequences = con.execute(f"""
            SELECT pid, protein_seq, temperature FROM proteins
            INNER JOIN taxa ON proteins.taxid = taxa.taxid
            WHERE pid IN {tuple(ids['thermo_pid'])}
        """).df()
        actual_protein_pairs[T] = (meso_sequences, thermo_sequences)
    con.close()

    # run esm and check outputs
    random_meso_results = run_esm_on_df(random_meso_proteins)
    random_intervals = {}
    for T, random_thermo_proteins in random_thermo_proteins_dict.items():
        # run thermo proteins in esm
        random_thermo_results = run_esm_on_df(random_thermo_proteins)
        meso_tensor = torch.vstack(tuple(random_meso_results['tensor'].values))
        thermo_tensor = torch.vstack(tuple(random_thermo_results['tensor'].values))
        assert len(meso_tensor) == len(thermo_tensor)
        print(meso_tensor.shape)

        # compute distance between meso and thermo random proteins
        distances = (meso_tensor - thermo_tensor).pow(2).sum(axis=1).sqrt().numpy()
        logger.info(f"Mean, std of ESM distance for random pairs with thermohpile T>{T}: {distances.mean()}, {distances.std()}")

        # estimate CI on low distance
        res = scipy.stats.bootstrap(
            (distances,),
            statistic = lambda x: np.quantile(x, 0.05),
            confidence_level = 0.95,
        )
        interval = res.confidence_interval
        logger.info(f"Confidence interval on 5th percentile ESM distance: {interval}")
        random_intervals[T] = interval

    # repeat the process, but use actual protein pairs
    actual_intervals = {}
    for T, (meso_sequences, thermo_sequences) in actual_protein_pairs.items():
        # run thermo proteins in esm
        meso_results = run_esm_on_df(meso_sequences)
        thermo_results = run_esm_on_df(thermo_sequences)
        meso_tensor = torch.vstack(tuple(meso_results['tensor'].values))
        thermo_tensor = torch.vstack(tuple(thermo_results['tensor'].values))
        assert len(meso_tensor) == len(thermo_tensor)

        # compute distance between meso and thermo random proteins
        distances = (meso_tensor - thermo_tensor).pow(2).sum(axis=1).sqrt().numpy()
        logger.info(f"Mean, std of ESM distance for actual pairs with thermohpile T>{T}: {distances.mean()}, {distances.std()}")

        # estimate CI on low distance
        res = scipy.stats.bootstrap(
            (distances,),
            statistic = lambda x: np.quantile(x, 0.95),
            confidence_level = 0.95,
        )
        interval = res.confidence_interval
        logger.info(f"Confidence interval on 95th percentile ESM distance: {interval}")
        actual_intervals[T] = interval

    # make plot of distance interval vs T
    fig, ax = plt.subplots()
    random_low = [i[0] for i in random_intervals.values()]
    random_high = [i[1] for i in random_intervals.values()]
    actual_low = [i[0] for i in actual_intervals.values()]
    actual_high = [i[1] for i in actual_intervals.values()]
    ax.fill_between(temps, random_low, random_high, alpha=0.5, label="5pctile random protein pair distance")
    ax.fill_between(temps, actual_low, actual_high, alpha=0.5, label="95pctile actual protein pair distance")
    ax.set_xlabel('Thermophilic temperature (C)')
    ax.set_ylabel('ESM Distance')
    ax.set_title('5th percentile ESM distance CI for meso <40C')
    fig.savefig('./data/plots/esm_distance_ci.png', dpi=300, bbox_inches='tight'
    )



    
