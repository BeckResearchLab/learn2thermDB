"""Random utilies."""
import logging
import distributed

def start_logger_if_necessary(logger_name: str, log_file: str, log_level, filemode: str = 'a'):
    """Quickly configure and return a logger that respects parallel processes.
    
    Parameters
    ----------
    logger_name : str
        name of logger to start or retrieve
    log_file : str
        path to file to log to
    log_level
        log level to respect
    worker: str
        name of worker using this logger
    filemode : str
        mode to apply to log file eg "a" for append
    """
    logger = logging.getLogger(logger_name)
    if len(logger.handlers) == 0:
        logger.setLevel(log_level)
        fh = logging.FileHandler(log_file, mode=filemode)
        fh.setFormatter(logging.Formatter("%(filename)s - %(asctime)s %(levelname)-8s %(message)s"))
        logger.addHandler(fh)
    return logger

class WorkerLogger(logging.Logger):
    """Logger that will supply worker name when logging."""
    def _log(self, level, msg, args, exc_info=None, extra=None):
        if extra is None:
            worker = distributed.get_worker()
            extra = {'worker': worker.name}
        super(WorkerLogger, self)._log(level, msg, args, exc_info, extra)