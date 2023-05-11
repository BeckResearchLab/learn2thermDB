"""Compare ESM distance of observed protein pairs to random pairs.

Actual pairs should be statistically significantly different than random pairs.
"""
import os
import pathlib
import torch
from esm import Alphabet, FastaBatchedDataset, pretrained, MSATransformer

import learn2therm.blast
import learn2therm.utils

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
    output_dir: pathlib.Path,
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
        model = model.cuda()
        logger.info("Transferred model to GPU")

    dataset = FastaBatchedDataset.from_file(fasta_file)
    batches = dataset.get_batch_indices(toks_per_batch, extra_toks_per_seq=1)
    data_loader = torch.utils.data.DataLoader(
        dataset, collate_fn=alphabet.get_batch_converter(truncation_seq_length), batch_sampler=batches
    )
    logger.info(f"Read {fasta_file} with {len(dataset)} sequences")

    output_dir.mkdir(parents=True, exist_ok=True)
    return_contacts = "contacts" in include

    assert all(-(model.num_layers + 1) <= i <= model.num_layers for i in repr_layers)
    repr_layers = [(i + model.num_layers + 1) % (model.num_layers + 1) for i in repr_layers]

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
                output_file = output_dir / f"{label}.pt"
                output_file.parent.mkdir(parents=True, exist_ok=True)
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
    return result

if __name__ == '__main__':

    # create fasta files and run ESM on both, compare distance
