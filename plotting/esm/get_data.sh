#!/bin/bash
curl -O 'https://raw.githubusercontent.com/facebookresearch/esm/main/scripts/atlas/v2023_02/full/esm2_embeddings/tm_.90_1_plddt_.90_1.txt'
aria2c --dir=./data/atlas/ --input-file='./data/tm_.90_1_plddt_.90_1.txt'
