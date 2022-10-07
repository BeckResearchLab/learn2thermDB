"""Assigns mesophilic, thermophilic labels.

A method for dealing with taxa with data but not OGT should be passed.
"""
import ast
from locale import normalize
import logging
import os

import matplotlib.pyplot as plt
import pandas as pd
from yaml import dump as yaml_dump
from yaml import safe_load as yaml_load

import learn2therm.utils

# get environmental variables
if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
    LOGLEVEL = getattr(logging, LOGLEVEL)
else:
    LOGLEVEL = logging.INFO
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'

def string_temp_to_float(temp: str):
    """Get floating point Temperature from a string.
    
    Some temperatures were returned as ranges, so the center is taken.
    """
    Ts = [float(T) for T in temp.split('-')]
    if len(Ts) == 1:
        return Ts[0]
    elif len(Ts) == 2:
        return sum(Ts)/2.0
    else:
        raise ValueError(f"Cannot interpret temp input {temp}")

"""Below are options to be called on each raw OGT input."""
def true_ogt_only(ogt_input):
    """Only take Actualy OGT's from bacdive as OGTs"""
    if type(ogt_input) == str:
        # we got OGT from bacdive
        return string_temp_to_float(ogt_input)
    elif type(ogt_input) == int or type(ogt_input) == float:
        return float(ogt_input)
    elif ogt_input == None:
        # we did not get anything from bacdive
        return None
    elif type(ogt_input) == dict:
        # we got something from bacdive but not an OGT
        return None 
    else:
        raise ValueError(f"Cannot interpret OGT inout {ogt_input}")

if __name__ == '__main__':
    # connect to log file
    logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')

    # DVC tracked parameters
    with open("./pipeline/s1_data_processing_params.yaml", "r") as stream:
        params = yaml_load(stream)['label_taxa']
    logger.info(f"Loaded parameters: {params}")

    # load the filepaths and their identities
    taxa = pd.read_csv('./data/taxa/taxa_info_and_ogt.csv', index_col=0,  usecols=[0,7])['ogt_raw'].fillna(value='None')
    

    # convert strings in csv to python obj
    logger.info("Converting raw OGT records saved to file to python objects")
    def eval_handled(string_in):
        try:
            return ast.literal_eval(string_in)
        except:
            return string_in
    taxa = taxa.apply(eval_handled)

    # apply the metric
    logger.info("Computing OGT from raw records")
    if params['ogt_determination_method'] == 'true ogt only':
        ogts = taxa.apply(true_ogt_only)
    ogts.name = 'ogt'
    ogts = ogts.astype(float)
    
    # extract meso or thermo based on threshold
    logger.info("Labeling thermophiles")
    thermo_bools = ogts.map(lambda ogt: ogt > params['ogt_threshold'], na_action='ignore')
    thermo_bools.name = "thermophile_label"
    
    # save to file
    labels = pd.DataFrame([ogts, thermo_bools]).T
    labels['thermophile_label'] = labels['thermophile_label'].astype('boolean')
    print(labels.columns)
    labels.to_csv('./data/taxa/labels.csv')

    # save metrics
    metrics = {}
    metrics['n_meso'] = int((labels['thermophile_label'] == False).sum())
    metrics['n_therm'] = int((labels['thermophile_label'] == True).sum())
    with open('./data/metrics/s1.0_metrics.yaml', "w") as stream:
        yaml_dump(metrics, stream)
    
    # save plot
    fig, ax = plt.subplots(figsize=(5,5))
    ax.set_xlabel('OGT [C]')
    labels['ogt'].plot.hist(bins=15, ax=ax)
    ax.vlines(x=[params['ogt_threshold']], ymin=0, ymax=ax.get_ylim()[1], colors=['r'])
    plt.savefig('./data/plots/ogt_hist.png', bbox_inches='tight', dpi=250)
