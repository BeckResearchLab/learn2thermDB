"""Select a random set of embeddings to use for tsne"""

import numpy as np
import os
import tqdm

N_SAMPLE = 500000

files = os.listdir('./data/atlas_process/')
files = [f for f in files if f.endswith('.npy')]

sizes = []
for f in tqdm.tqdm(files):
    array = np.load('./data/atlas_process/' + f)
    sizes.append(len(array))
print(f"Embeddings of size {array.shape[1]} loaded for {sum(sizes)} proteins")
with open('log', 'a') as f:
    f.write(f"Number of Atlas proteins: {sum(sizes)}\n")
    f.write(f"Fraction of Atlas proteins used: {N_SAMPLE/sum(sizes)}\n")

# generate random indexes
indexes = np.random.choice(np.arange(sum(sizes)), size=N_SAMPLE, replace=False)
# sort indexes
indexes = np.sort(indexes)

# loop throuhg files, offset the indexes by the running file size, and select the embeddings
vecs = []
for i, f in enumerate(tqdm.tqdm(files)):
    array = np.load('./data/atlas_process/' + f)
    # get only the indexes that are within this file
    file_indexes = indexes[(indexes >= sum(sizes[:i])) & (indexes < sum(sizes[:i+1]))]
    # offset the indexes by the running file size
    file_indexes = file_indexes - sum(sizes[:i])
    # select the embeddings
    vecs.append(array[file_indexes,:])
# concatenate the embeddings
vecs = np.vstack(vecs)
# save the embeddings
np.save('./data/atlas_subset.npy', vecs)
    