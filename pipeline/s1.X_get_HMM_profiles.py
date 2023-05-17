"""Ingest raw PFAM HMMs

Sources
-------
[1] R. D. Finn et al., “The Pfam protein families database: towards a more sustainable future,” 
Nucleic Acids Research, vol. 44, no. D1, pp. D279–D285, Jan. 2016, doi: 10.1093/nar/gkv1344.

[2] T. Paysan-Lafosse et al., “InterPro in 2022,” 
Nucleic Acids Research, vol. 51, no. D1, pp. D418–D427, Jan. 2023, doi: 10.1093/nar/gkac993.


Notes
-----
TODO:
Edit __name__: "__main__" with Evan

https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/
"""
# system dependecies
import datetime
from ftplib import FTP
import logging
import os
from tqdm import tqdm

# library dependecies
from yaml import dump as yaml_dump
from yaml import safe_load as yaml_load

# local dependencies
import learn2therm.database
import learn2therm.utils

try:
    EMAIL = os.environ['ENV_EMAIL']
except KeyError:
    raise KeyError('Must set environmental variables `ENV_EMAIL`')
if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
    LOGLEVEL = getattr(logging, LOGLEVEL)
else:
    LOGLEVEL = logging.INFO
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'
# get the logger in subprocesses
logger = learn2therm.utils.start_logger_if_necessary(LOGNAME, LOGFILE, LOGLEVEL, filemode='w')

# set up ftp
FTP_ADDRESS = 'ftp.ebi.ac.uk'
FTP_DIR = '/pub/databases/Pfam/current_release/'

def ftp_get_file_progress_bar(filename, endpoint_dir):
    ftp = FTP(FTP_ADDRESS)
    ftp.login(user="anonymous", passwd=EMAIL)
    ftp.cwd(FTP_DIR) 

    file_size = ftp.size(filename)
    with open(endpoint_dir+f'{filename}', 'wb') as file:
        pbar = tqdm(range(file_size))
        def write_file(data):
            file.write(data)
            pbar.n += len(data)
            pbar.refresh()

        ftp.retrbinary(f"RETR {filename}", write_file, blocksize=262144)
    ftp.close()

if __name__ == "__main__":
    # DVC tracked parameters
    with open("./params.yaml", "r") as stream:
        params = yaml_load(stream)['get_raw_data_HMM']
    logger.info(f"Loaded parameters: {params}")

    if not os.path.exists('./data/HMM'):
        os.mkdir('./data/HMM')

    # download the raw data into temporary files
    dir = './data/HMM'
    addresses = ['Pfam-A.hmm.gz']
    if params['dev_only_one_uniprot_file']:
        addresses = addresses[:1]
        logger.info(f"Downloading only {addresses}")

    # download each file from uniprot
    for i, address in enumerate(addresses):
        if address in os.listdir(dir):
            logger.info(f"Address exists, skipping: {address}")
            continue
        ftp_get_file_progress_bar(address, endpoint_dir='./data/HMM/')
        logger.info(f"Completed download of {address}")

    # save metrics
    date_pulled = str(datetime.datetime.now().strftime("%m/%d/%Y"))
    with open('./data/HMM/HMM_pulled_timestamp', 'w') as file:
        file.write(date_pulled)
                      
    metrics = {}
    metrics['HMM_pulled_date'] = date_pulled
    with open('./data/metrics/s0.1_metrics.yaml', "w") as stream:
        yaml_dump(metrics, stream)
