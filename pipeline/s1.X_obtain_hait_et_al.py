"""Download Hait data from the internet


Sources
-------
[1] [1] S. Hait, S. Mallik, S. Basu, and S. Kundu, “Finding the generalized molecular principles of protein thermal stability,” 
Proteins: Structure, Function, and Bioinformatics, vol. 88, no. 6, pp. 788–808, 2020, doi: 10.1002/prot.25866.




Notes
-----
https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fprot.25866&file=prot25866-sup-0001-DataS1.xlsx

TODO: check copyright stuff for this?
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
import requests
from urllib.parse import urlparse, parse_qs

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

# URL of the data file
url = 'https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fprot.25866&file=prot25866-sup-0001-DataS1.xlsx'

def download_file(url, local_file_path):
    """
    Download a file from a given URL.

    Parameters
    ----------
    url : str
        The URL of the file to download.
    local_file_path : str
        The local file path where the file will be saved.

    Returns
    -------
    local_file_path : str
        The local path where the file has been saved.

    Raises
    ------
    HTTPError
        If there was an HTTP error while making the request.
    """
    # Because Wiley is stuipd
    headers = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.3',}
    response = requests.get(url, stream=True, headers=headers)
    response.raise_for_status()
    with open(local_file_path, 'wb') as f:
        for chunk in response.iter_content(chunk_size=8192):
            if chunk:
                f.write(chunk)
    return local_file_path

def get_data(url, local_file_path):
    """
    Retrieve data file from a URL or a local file path. 
    
    If the local file doesn't exist, try to download it. If the download fails,
    fall back to the local file.

    Parameters
    ----------
    url : str
        The URL of the file to download.
    local_file_path : str
        The local file path where the file will be saved.

    Returns
    -------
    local_file_path : str
        The local path where the file has been saved.

    Raises
    ------
    Exception
        If no data file found at the local file path.
    """
    # If the local file doesn't exist, try to download it
    if not os.path.isfile(local_file_path):
        try:
            download_file(url, local_file_path)
        except requests.exceptions.HTTPError as e:
            logger.warning(f'Warning: failed to download file from {url} due to HTTP error: {e}. Falling back to local copy')

    # At this point. the local file should exist (either we just downloaded it, or it was already there)
    if os.path.isfile(local_file_path):
        return local_file_path
    else:
        raise Exception(f'Error: no data file found at {local_file_path}')

if __name__ == "__main__":
    # DVC tracked parameters
    with open("./params.yaml", "r") as stream:
        params = yaml_load(stream)['get_raw_data_HMM']
    logger.info(f"Loaded parameters: {params}")

    try:
        os.makedirs('./data/Hait', exist_ok=True)
    except OSError as e:
        logger.error(f'Error creating directory: {e}')

    # Parse the URL into components
    parsed_url = urlparse(url)

    # Parse the query string into a dictionary
    query_dict = parse_qs(parsed_url.query)

    # Get the filename from the 'file' query parameter
    filename = query_dict.get('file', [''])[0]

    logger.info("Parse URL to find file")

    # Path to save the data file locally
    local_file_path = os.path.join('./data/Hait', filename)

    get_data(url, local_file_path)
    logger.info("downloaded data from the internet")

    # save metrics
    date_pulled = str(datetime.datetime.now().strftime("%m/%d/%Y"))
    with open('./data/Hait/Hait_pulled_timestamp', 'w') as file:
        file.write(date_pulled)
                      
    metrics = {}
    metrics['Hait_pulled_date'] = date_pulled
    with open('./data/metrics/s1.X_metrics.yaml', "w") as stream:
        yaml_dump(metrics, stream)