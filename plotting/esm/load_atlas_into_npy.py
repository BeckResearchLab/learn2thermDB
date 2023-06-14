"""Run multicore tsne data on ESM2 Atlas proteins
and the proteins within the learn2therm dataset
"""
# from MulticoreTSNE import MulticoreTSNE as TSNE
import numpy as np
import os
import tqdm

# load the esm atlas data
ATLAS_DATA_DIR = './data/atlas/'
files = os.listdir(ATLAS_DATA_DIR)
files = [f for f in files if f.endswith('.npz')]

def do_one(f):
    with np.load(os.path.join(ATLAS_DATA_DIR, f)) as data_file:
        n_components = len(data_file.values())
        vecs = []
        with tqdm.tqdm(total=n_components) as pbar:
            for v in data_file.values():
                vecs.append(v.reshape(1,-1))
                pbar.update(1)
    return np.vstack(vecs)

for f in tqdm.tqdm(files):
    if f+'.npy' not in os.listdir('./data/atlas_process/'):
        array = do_one(f)
        np.save('./data/atlas_process/' + f, array)
    else:
        pass

size = 0
for f in os.listdir('./data/atlas_process/'):
    if f.endswith('.npy'):
        array = np.load('./data/atlas_process/' + f)
        size += len(array)
print(f"Embeddings of size {array.shape[1]} loaded for {size} proteins")

