"""Congregates csv files into duckdb file using learn2therm.databsae.L2TDatabase class

No parameters here, just operating on already created data.
"""

import learn2therm.database
import learn2therm.utils
import os
import logging

if 'LOGLEVEL' in os.environ:
    LOGLEVEL = os.environ['LOGLEVEL']
    LOGLEVEL = getattr(logging, LOGLEVEL)
else:
    LOGLEVEL = logging.INFO
LOGNAME = __file__
LOGFILE = f'./logs/{os.path.basename(__file__)}.log'

if __name__ == "__main__":
    logger = learn2therm.utils.start_logger_if_necessary("", LOGFILE, LOGLEVEL, filemode='w')

    db = learn2therm.database.L2TDatabase.from_files(files_path='./data/', db_path='./data/database')