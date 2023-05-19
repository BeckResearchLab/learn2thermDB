"""Ingest raw PFAM HMMs
TODO: how to unzip

Sources
-------
[1] R. D. Finn et al., “The Pfam protein families database: towards a more sustainable future,” 
Nucleic Acids Research, vol. 44, no. D1, pp. D279–D285, Jan. 2016, doi: 10.1093/nar/gkv1344.

[2] T. Paysan-Lafosse et al., “InterPro in 2022,” 
Nucleic Acids Research, vol. 51, no. D1, pp. D418–D427, Jan. 2023, doi: 10.1093/nar/gkac993.


Notes
-----
TODO
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

## get environmental variables
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

def download_ftp_file(server, remote_file, local_file):
    ftp = FTP(server)
    ftp.login(user="anonymous", passwd='')
    ftp.cwd('/')  # Set the current directory (if needed)

    with open(local_file, 'wb') as file:
        ftp.retrbinary('RETR ' + remote_file, file.write)

    ftp.quit()

if __name__ == "__main__":
    # DVC tracked parameters
    with open("./params.yaml", "r") as stream:
        params = yaml_load(stream)['get_raw_data_HMM']
    logger.info(f"Loaded parameters: {params}")

    if not os.path.exists('./data/HMM'):
        os.mkdir('./data/HMM')

    # download the raw data into temporary files 
    local_dir = './data/HMM/Pfam-A.hmm.gz'
    remote_file = '/pub/databases/Pfam/current_release/Pfam-A.hmm.gz'
    
    download_ftp_file(FTP_ADDRESS, remote_file, local_dir)
    logger.info("download from FTP")

    # save metrics
    date_pulled = str(datetime.datetime.now().strftime("%m/%d/%Y"))
    with open('./data/HMM/HMM_pulled_timestamp', 'w') as file:
        file.write(date_pulled)
                      
    metrics = {}
    metrics['HMM_pulled_date'] = date_pulled
    with open('./data/metrics/s0.1_metrics.yaml', "w") as stream:
        yaml_dump(metrics, stream)
