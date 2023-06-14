"""Select a sample of l2t proteins that occur in pairs equal to the fraction of Atlas used,
and compute and save the ESM embeddings.
"""
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import torch
import numpy as np
import duckdb as ddb
import pandas as pd
from esm import Alphabet, FastaBatchedDataset, pretrained, MSATransformer
import pathlib
from typing import List
import tempfile
import os

import logging
logging.basicConfig(filename='l2t_esm.log',
                    filemode='w',
                    format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                    datefmt='%H:%M:%S',
                    level=logging.INFO)

logger = logging.getLogger('l2t_esm')

FRAC = 0.06010428815245909 # fraction of original dataset to use

def process_fasta_file(
    model_location: str,
    fasta_file: pathlib.Path,
    toks_per_batch: int = 18000,
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
        logger.info('CUDA FOUND!')
        model = model.cuda()
    else:
        logger.info('NOCUDA')
    dataset = FastaBatchedDataset.from_file(fasta_file)
    batches = dataset.get_batch_indices(toks_per_batch, extra_toks_per_seq=1)
    data_loader = torch.utils.data.DataLoader(
        dataset, collate_fn=alphabet.get_batch_converter(truncation_seq_length), batch_sampler=batches
    )
    logger.info(f"Read {fasta_file} with {len(dataset)} sequences")

    return_contacts = "contacts" in include

    assert all(-(model.num_layers + 1) <= i <= model.num_layers for i in repr_layers)
    repr_layers = [(i + model.num_layers + 1) % (model.num_layers + 1) for i in repr_layers]

    with torch.no_grad():
        tmp_results = []
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
                tmp_results.append(result)
            torch.save(tmp_results, f'./data/pt/l2t_esm_results_batch_{batch_idx}.pt')
            tmp_results = []
    return

def run_esm_on_df(dataframe, model='esm2_t36_3B_UR50D'):
    """Create fasta files and run ESM"""

    
    # write fasta file
    fasta_file_loc = './data/esm_fasta_input.fasta'
    if os.path.exists(fasta_file_loc):
        fasta_file = fasta_file_loc
    else:
        records = []
        ids_ = []
        for i, row in dataframe.iterrows():
            if row['pid'] in ids_:
                continue
            else:
                record = SeqRecord(Seq(row['protein_seq']), id=row['pid'], description='')
                records.append(record)
                ids_.append(row['pid'])
        fasta_file = open(fasta_file_loc, 'w')
        SeqIO.write(records, fasta_file, 'fasta')
        fasta_file.close()
        fasta_file = fasta_file.name
        del records
    fp = pathlib.Path(fasta_file)
    logger.info(f"Using fasta file {fasta_file}")

    del dataframe

    # run ESM
    result = process_fasta_file(
        model_location=model,
        fasta_file=fp,
        toks_per_batch=4096,
        repr_layers=[-1],
        include=["mean"],
        truncation_seq_length=250,
        nogpu=False,
    )

    logger.info(f"Ran ESM on {fp}")
    return result

if __name__ == '__main__':

    # get data from duck db
    con = ddb.connect('../../data/database.ddb', read_only=True)
    count = con.execute("SELECT COUNT(pid) FROM proteins WHERE proteins.pid IN (SELECT DISTINCT(pairs.meso_pid) FROM pairs) OR proteins.pid IN (SELECT DISTINCT(pairs.thermo_pid) FROM pairs)").fetchone()[0]
    n_to_keep = int(FRAC * count)

    data = con.execute(f"SELECT pid, protein_seq FROM proteins WHERE proteins.pid IN (SELECT DISTINCT(pairs.meso_pid) FROM pairs) OR proteins.pid IN (SELECT DISTINCT(pairs.thermo_pid) FROM pairs) ORDER BY RANDOM() LIMIT {n_to_keep}").df()
    logger.info(f'Got {len(data)} proteins')

    # run ESM
    run_esm_on_df(data, model='esm2_t36_3B_UR50D')
    
