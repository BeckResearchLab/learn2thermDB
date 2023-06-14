"""Run multithreaded t-SNE on a dataset."""
from MulticoreTSNE import MulticoreTSNE as TSNE
import numpy as np
import time
import matplotlib.pyplot as plt
import torch
import numpy as np
import pandas as pd
import os

import logging
logging.basicConfig(filename='tsne.log',
                    filemode='w',
                    format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                    datefmt='%H:%M:%S',
                    level=logging.INFO)

logger = logging.getLogger('l2t_esm')

FRAC = 0.06010428815245909 # fraction of original dataset to use, here only applies to Hait, L2t is already subsetted
SAMPLE = None

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

if __name__ == "__main__":

     # if we have already run tsne, load it
    if not os.path.exists('./data/tsne.npy'):

        # load the l2t embeddings
        if os.path.exists('./data/l2t_subset.npy'):
            l2t = np.load('./data/l2t_subset.npy')
        else:
            result_files = os.listdir('./data/pt')
            result_files = [f'./data/pt/{f}' for f in result_files if f.endswith('.pt')]
            results = []
            logger.info(f"Loading l2t ESM embeddings from {len(result_files)} files")
            for result_file in result_files:
                results.extend(torch.load(result_file))
            result = parse_results(results)['tensor']
            l2t = np.vstack(result.values)
            np.save('./data/l2t_subset.npy', l2t)
        logger.info(f"l2t Embeddings of shape {l2t.shape}")

        # load the hait embeddings
        if os.path.exists('./data/hait_subset.npy'):
            hait = np.load('./data/hait_subset.npy')
        else:
            result_files = os.listdir('./data/pt_hait')
            result_files = [f'./data/pt_hait/{f}' for f in result_files if f.endswith('.pt')]
            results = []
            logger.info(f"Loading Hait ESM embeddings from {len(result_files)} files")
            for result_file in result_files:
                results.extend(torch.load(result_file))
            result = parse_results(results)['tensor']
            hait = np.vstack(result.values)
            np.save('./data/hait_subset.npy', hait)
        logger.info(f"Hait Embeddings of shape {hait.shape}")
        hait = hait[np.random.choice(hait.shape[0], int(hait.shape[0] * FRAC), replace=False)]
        logger.info(f"Sampling fraction {FRAC} of hait embeddings to retain dataset size rations, new shape {hait.shape}")


        # load the atlas embeddings
        atlas = np.load('./data/atlas_subset.npy')
        logger.info(f"Atlas Embeddings of shape {atlas.shape}")
        
        # sampling
        if SAMPLE:
            l2t = l2t[np.random.choice(l2t.shape[0], int(l2t.shape[0] * SAMPLE), replace=False)]
            atlas = atlas[np.random.choice(atlas.shape[0], int(atlas.shape[0] * SAMPLE), replace=False)]
            logger.info(f"Sampling fraction {SAMPLE} of embeddings for development")

        # stick them otgether and retain the labels
        labels = np.array(['l2t'] * l2t.shape[0] + ['atlas'] * atlas.shape[0] + ['hait'] * hait.shape[0])
        x = np.vstack([l2t, atlas, hait])
        logger.info(f"Combined embeddings of shape {x.shape}")
        
   
        tsne = TSNE(
            2,
            n_jobs=40,
            n_iter_early_exag=500,
            n_iter=2000,
            verbose=1)
        t0 = time.time()
        print("running tsne...")
        p = tsne.fit_transform(x)
        t1 = time.time()
        print(f"tsne took {(t1-t0)/60} minutes")
        np.save('./data/tsne.npy', p)
        np.save('./data/tsne_labels.npy', labels)
    else:
        p = np.load('./data/tsne.npy')
        labels = np.load('./data/tsne_labels.npy')

    # get the data back out seperated by source
    l2t = p[labels == 'l2t']
    atlas = p[labels == 'atlas']
    hait = p[labels == 'hait']

    # plot the map
    plt.style.use('dark_background')
    fig, ax = plt.subplots(figsize=(10,10))
    ax.scatter(atlas[:,0], atlas[:,1], c='#0079FF', s=12, alpha=1, linewidths=1, edgecolors="#0062CF", label='ESM Atlas')
    ax.scatter(l2t[:,0], l2t[:,1], c='#ffff99', s=8, alpha=0.5, linewidths=1, edgecolors="#FFE03D", label='learn2therm')
    ax.scatter(hait[:,0], hait[:,1], c='#FF5416', s=15, alpha=1.0, linewidths=0.0, label='Hait et al.')
    # remove the spline and labels, only want internal plot
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    leg = plt.legend()
    for lh in leg.legendHandles: 
        lh.set_alpha(1)
    ax.set_xticks([])
    ax.set_yticks([])
    plt.savefig('./tsne.png', dpi=300)

