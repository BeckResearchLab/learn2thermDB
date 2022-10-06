"""Use NCBI TAXid to search bacdive for OGT. and record the results"""
import logging
import os
import time

from joblib import delayed, Parallel
import pandas as pd
from yaml import safe_load as yaml_load
from yaml import dump as yaml_dump

import learn2therm.bacdive
import learn2therm.io
import learn2therm.utils

# get environmental variables
try:
    USERNAME = os.environ['ENV_EMAIL']
    PASSWORD = os.environ['BACDIVE_PASSWORD']
except KeyError:
    raise KeyError('Must set environmental variables `ENV_EMAIL` and `BACDIVE_PASSWORD`')
if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
    LOGLEVEL = getattr(logging, LOGLEVEL)
else:
    LOGLEVEL = logger.inFO
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'

def get_one_ogt_from_taxid(client: learn2therm.bacdive.BacdiveClient, taxid: int):
    """Attempts to retrieve an OGT from BACDive from TAXID, and handles processing if failed.
    Returns
    -------
    bacdive id, OGT or list of growth temperatures
    """
    logger.info(f"Start bacdive search for TAXID {taxid}")
    # make the query for bacdiveID from taxid
    if client.getIDByNCBITaxID(taxid) == 1: # returns 1 if found items, 0 if not
        bacdive_output = next(client.retrieve())
        bacdive_id = bacdive_output['General']['BacDive-ID']
        # we must parse for temperatures carefully, multiple experiments may be listed.
        # lets get them all and we can process the results into a single label downstream
        temp_labels = {'growth': [], 'max': None, 'min': None}
        logger.debug(f"Growth conditions for taxid {taxid}: {bacdive_output['Culture and growth conditions']}")
        if 'culture temp' in bacdive_output['Culture and growth conditions']:
            if type(bacdive_output['Culture and growth conditions']['culture temp']) == list:
                temp_dicts = bacdive_output['Culture and growth conditions']['culture temp']
            elif type(bacdive_output['Culture and growth conditions']['culture temp']) == dict:
                temp_dicts = [bacdive_output['Culture and growth conditions']['culture temp']]
            else:
                raise ValueError(f"Unknown type {type(bacdive_output['Culture and growth conditions']['culture temp'])} for 'culture temp' in 'Culture and growth conditions' of bacdive record {bacdive_id}")

            for temp_dic in temp_dicts:
                # some instances have too little information to use, just skip them
                if 'type' not in temp_dic:
                    continue
                # if we hit OGT just take it and run
                elif temp_dic['type'] == 'optimum':
                    temp_labels = str(temp_dic['temperature'])
                    return (bacdive_id, temp_labels)

                # only take measurements that resulted in positive growth : we don't care about a temperature
                # if it killed the bugs
                elif temp_dic['type'] == 'growth':
                    if 'growth' in temp_dic and temp_dic['growth'] not in ['yes', 'positive']:
                        continue # we want only positive growth
                    if 'temperature' in temp_dic:
                        temp_labels['growth'].append(temp_dic['temperature'])
                    elif 'range' in temp_dic:
                        temp_labels['growth'].append(temp_dic['range'])
                    else:
                        # if we won't have a tmperature or at least a label, why are we here
                        continue
                elif temp_dic['type'] in ['maximum', 'max']:
                    temp_labels['max'] = temp_dic['temperature']
                elif temp_dic['type'] in ['minimum', 'min']:
                    temp_labels['min'] = temp_dic['temperature']
                else:
                    continue
            # check if we actually added anything, if not replace with None
            if len(temp_labels['growth']) == 0 and temp_labels['min'] == None and temp_labels['max'] == None:
                temp_labels = None
            return (bacdive_id, temp_labels) # we found a bacdive and it had temperatures
        else:
            return (bacdive_id, None) # we found a bacdive but no temperatures
    else:
        return (None, None) # we didn't even find a bacdive

def process_one_gbff_zip_file(filepath: str):
    """Loads, reads, queries bacdive, and formats output for a single file.
    
    Parameters
    ----------
    filepath : str
        path to zipped gbff file
    
    Returns
    -------
    dict of information from the gbff file and bacdive
    """
    client = learn2therm.bacdive.BacdiveClient(USERNAME, PASSWORD)

    # get the logger in subprocesses
    logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL)

    logger.debug(f"Extracting file {filepath}")
    refseq_records = learn2therm.io.seq_io_gnuzipped(filepath, 'genbank')

    # check that each record has the same taxid fo sanity
    taxids = []
    for i, r in enumerate(refseq_records):
        if r.features[0].type != 'source':
            raise ValueError("Feature at position 0 of filepath record {i} in {filepath} is not a source feature, cannot check taxid")
        else:
            taxids.append(dict(r.features[0].qualifiers)['db_xref'][-1])
    if len(set(tuple(taxids))) > 1:
        raise ValueError(f"Inconsistant taxids accross records in {filepath}: {set(tuple(taxids))}")

    # we only need one record to get metadata
    record = refseq_records[0]
    source_feature = record.features[0]
    
    # sometimes there are multipble database references, we want the one that startwith "taxon" to get taxid
    taxid = None
    for ref in dict(source_feature.qualifiers)['db_xref']:
        if ref.startswith('taxon'):
            taxid = int(ref.split(':')[-1])
    if taxid == None:
        raise ValueError(f"Could not interpret taxid from {dict(source_feature.qualifiers)['db_xref']}")

    # save information from the refseq file
    output = {
        'taxid': taxid,
        'record_name': record.name,
        'filepath': filepath,
        'taxonomy': ' '.join(record.annotations['taxonomy']),
        'organism': record.annotations['organism'],
    }

    # query bacdive and save the results
    bacdive_id, ogt_raw = get_one_ogt_from_taxid(client, taxid)
    if bacdive_id == None:
        logger.info(f"Could not find BACDIVE record for taxid {taxid}")
    elif ogt_raw == None:
        logger.info(f"Could not valid temperature experiments for taxid {taxid}")
    elif type(ogt_raw) == str:
        logger.info(f"Found OGT for taxid {taxid}")
    else:
        logger.info(f"Found temperatures that led to growth for TAXID {taxid}: {ogt_raw['growth']} with min and max ({ogt_raw['min'], ogt_raw['max']})")

    # add info to the output
    output['bacdive_id'] = bacdive_id
    output['ogt_raw'] = ogt_raw
    return output

if __name__ == '__main__':
    # connect to log file
    logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')

    # DVC tracked parameters
    with open("./pipeline/s0_data_extraction_params.yaml", "r") as stream:
        params = yaml_load(stream)['get_bacdive_ogt']
    logger.info(f"Loaded parameters: {params}")

    # # connect to bacdive to retrieve data
    # client = learn2therm.bacdive.BacdiveClient(USERNAME, PASSWORD)
    # logger.info(f"Connected to BacDive")

    # get the filepaths to execute on
    file_list = []
    for path, subdirs, files in os.walk('./data/refseq'):
        for file in files:
            file_list.append(os.path.join(path, file))
    
    # select a sample if specified
    if params['n_sample']:
        file_list = file_list[:params['n_sample']]
        logger.info(f"Using only the first {params['n_sample']} files.")
    
    # execute in parallel
    outputs = Parallel(n_jobs=params['n_jobs'])(delayed(lambda file: process_one_gbff_zip_file(file))(file) for file in file_list)
    df = pd.DataFrame(outputs)
    
    # save data
    os.makedirs('./data/taxa', exist_ok=True)
    df.to_csv('./data/taxa/taxa_info_and_ogt.csv')

    # save metrics
    metrics = {}
    metrics['n_taxa'] = len(df)
    metrics['n_have_bacdive'] = int((~df['bacdive_id'].isna()).sum())
    metrics['n_have_ogt'] = int(df['ogt_raw'].apply(lambda x: type(x) == str).sum())
    metrics['n_have_growth_temp'] = int(df['ogt_raw'].apply(lambda x: type(x) == dict).sum())
    with open('./data/metrics/s0.1_metrics.yaml', "w") as stream:
        yaml_dump(metrics, stream)


