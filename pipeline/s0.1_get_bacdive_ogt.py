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


def get_one_ogt_from_taxid(client: learn2therm.bacdive.BacdiveClient, taxid: int):
    """Attempts to retrieve an OGT from BACDive from TAXID, and handles processing if failed.
    Returns
    -------
    bacdive id, OGT or list of growth temperatures
    """
    # make the query for bacdiveID from taxid
    if client.getIDByNCBITaxID(taxid) == 1: # returns 1 if found items, 0 if not
        bacdive_output = next(client.retrieve())
        bacdive_id = bacdive_output['General']['BacDive-ID']
        # we must parse for temperatures carefully, multiple experiments may be listed.
        # lets get them all and we can process the results into a single label downstream
        temp_labels = []
        logging.debug(f"Growth conditions for taxid {taxid}: {bacdive_output['Culture and growth conditions']}")
        if 'culture temp' in bacdive_output['Culture and growth conditions']:
            if type(bacdive_output['Culture and growth conditions']['culture temp']) == list:
                temp_dicts = bacdive_output['Culture and growth conditions']['culture temp']
            elif type(bacdive_output['Culture and growth conditions']['culture temp']) == dict:
                temp_dicts = [bacdive_output['Culture and growth conditions']['culture temp']]
            else:
                raise ValueError(f"Unknown type {type(bacdive_output['Culture and growth conditions']['culture temp'])} for 'culture temp' in 'Culture and growth conditions' of bacdive record {bacdive_id}")

            for temp_dic in temp_dicts:
                # if we hit OGT just take it and run
                if temp_dic['type'] == 'optimum':
                    temp_labels = int(temp_dic['temperature'])
                    break

                # only take measurements that resulted in positive growth : we don't care about a temperature
                # if it killed the bugs
                if temp_dic['growth'] == 'positive':
                    if 'temperature' in temp_dic:
                        temp_labels.append(temp_dic['temperature'])
                    elif 'range' in temp_dic:
                        temp_labels.append(temp_dic['range'])
                    else:
                        # if we won;t have a tmperature or at least a label, why are we here
                        continue
                else:
                    continue
            return (bacdive_id, temp_labels) # we found a bacdive and it had temperatures
        else:
            return (bacdive_id, None) # we found a bacdive but no temperatures
    else:
        return (None, None) # we didn't even find a bacdive

def process_one_gbff_zip_file(filepath: str, client: learn2therm.bacdive.BacdiveClient):
    """Loads, reads, queries bacdive, and formats output for a single file.
    
    Parameters
    ----------
    filepath : str
        path to zipped gbff file
    
    Returns
    -------
    dict of information from the gbff file and bacdive
    """
    refseq_records = learn2therm.io.seq_io_gnuzipped(filepath, 'genbank')

    # check that each record has the same taxid fo sanity
    taxids = []
    for i, r in enumerate(refseq_records):
        if r.features[0].type != 'source':
            raise ValueError("Feature at position 0 of filepath record {i} in {filepath} is not a source feature, cannot check taxid")
        else:
            taxids.append(dict(r.features[0].qualifiers)['db_xref'][0])
    if len(set(tuple(taxids))) > 1:
        raise ValueError(f"Inconsistant taxids accross records in {filepath}: {set(tuple(taxids))}")

    # we only need one record to get metadata
    record = refseq_records[0]
    source_feature = record.features[0]
    taxid = int(dict(source_feature.qualifiers)['db_xref'][0].split(':')[1])

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
        logging.info(f"Could not find BACDIVE record for taxid {taxid}")
    elif ogt_raw == None:
        logging.info(f"Could not valid temperature experiments for taxid {taxid}")
    elif type(ogt_raw) != list:
        logging.info(f"Found OGT for taxid {taxid}")
    else:
        logging.info(f"Found temperatures that led to growth for TAXID {taxid}")

    # add info to the output
    output['bacdive_id'] = bacdive_id
    output['ogt_raw'] = ogt_raw
    return output

if __name__ == '__main__':
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
        LOGLEVEL = logging.INFO

    # connect to log file
    logger = logging.basicConfig(filename=f'./logs/{os.path.basename(__file__)}.log', level=LOGLEVEL)

    # DVC tracked parameters
    with open("./pipeline/s0_data_extraction_params.yaml", "r") as stream:
        params = yaml_load(stream)['get_bacdive_ogt']
    logging.info(f"Loaded parameters: {params}")

    # connect to bacdive to retrieve data
    client = learn2therm.bacdive.BacdiveClient(USERNAME, PASSWORD)
    logging.info(f"Connected to BacDive")

    # get the filepaths to execute on
    file_list = []
    for path, subdirs, files in os.walk('./data/refseq'):
        for file in files:
            file_list.append(os.path.join(path, file))
    
    # select a sample if specified
    if params['n_sample']:
        file_list = file_list[:params['n_sample']]
        logging.info(f"Using only the first {params['n_sample']} files.")
    
    # execute in parallel
    outputs = Parallel(n_jobs=params['n_jobs'])(delayed(lambda file: process_one_gbff_zip_file(file, client))(file) for file in file_list)
    df = pd.DataFrame(outputs)
    
    # save data
    os.makedirs('./data/taxa', exist_ok=True)
    df.to_csv('./data/taxa/taxa_info_and_ogt.csv')

    # save metrics
    metrics = {}
    metrics['n_taxa'] = len(df)
    metrics['n_have_bacdive'] = int((~df['bacdive_id'].isna()).sum())
    metrics['n_have_ogt'] = int(df['ogt_raw'].apply(lambda x: type(x) == int).sum())
    metrics['n_have_growth_temp'] = int(df['ogt_raw'].apply(lambda x: type(x) == list).sum())
    with open('./data/metrics/s0.1_metrics.yaml', "w") as stream:
        yaml_dump(metrics, stream)


