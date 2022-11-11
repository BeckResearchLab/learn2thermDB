"""Random utilies."""
import logging
import distributed


class WorkerLogger(logging.Logger):
    """Logger that will supply worker name when logging."""
    def _log(self, level, msg, args, exc_info=None, extra=None):
        if extra is None:
            extra = {}
        worker = distributed.get_worker()
        worker_name = worker.name

        extra.update({'worker': worker_name})
        super(WorkerLogger, self)._log(level, msg, args, exc_info, extra)     
        
def start_logger_if_necessary(logger_name: str, log_file: str, log_level, filemode: str = 'a', worker: bool = False):
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
    logger.setLevel(log_level)
    fh = logging.FileHandler(log_file, mode=filemode)
    if worker:
        worker_name = distributed.get_worker().name
        fh.setFormatter(logging.Formatter("%(filename)s %(worker)s - %(asctime)s %(levelname)-8s %(message)s"))
        if len(logger.handlers) == 0:
            logger.addHandler(fh)
        else:
            logger.handlers[-1] = fh
        logger = logging.LoggerAdapter(logger, extra={'worker': worker_name})
    else:
        fh.setFormatter(logging.Formatter("%(filename)s - %(asctime)s %(levelname)-8s %(message)s"))
        if len(logger.handlers) == 0:
            logger.addHandler(fh)
        else:
            logger.handlers[-1] = fh
    return logger