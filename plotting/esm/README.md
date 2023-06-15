# code for creating comparison of protein space between datasets


Data: "Atlas" is ESM protein embeddings with >90 TM and plddt, our learn2therm data, and hait data

## steps

1. create environment and install via `pip install -r requirements.txt`
2. Get atlas data: `./get_data.sh`
3. Convert data to fast loading numpy arrays `python load_atlas_into_npy.py`
4. Select a sample of Atlas to keep randomly, because TSNE is too expensive: `python select_atlas_subset.py`
5. Run esm on Hait proteins `python run_esm_on_hait.py` (makes sure gpu is available)
6. Select a sample of l2t data: `python make_l2t_fasta_for_esm.py`
7. Run esm on l2t proteins `bash run_esm_l2t.sh` (make sure gpu is available)
8. Run tsne `python run_tsne.py` (make sure many cpus are available, 32 was used here)
