import os
import sys
import time
import logging
import logging.handlers as handlers

formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

logger = logging.getLogger('main')
logger.setLevel(logging.INFO)
logger.propagate = False  # third-party libs (e.g. jpype/orekit) may configure the
                          # root logger; without this every line would be duplicated

ch = logging.StreamHandler(sys.stdout)  # log to the stdout
ch.setFormatter(formatter)
logger.addHandler(ch)


def init_file(log_path):
    """Attach the main.log file handler. Called by main.py once the output
    directory is known (--output-dir), so it cannot happen at import time.
    A previous log is removed first so every run starts fresh; sync clients
    (OneDrive) may hold a transient lock on it, hence the retry."""
    for attempt in range(5):
        try:
            if os.path.exists(log_path):
                os.remove(log_path)
            break
        except PermissionError:
            time.sleep(1.0)
    else:
        logger.warning(f'Could not remove previous {log_path} (file locked); appending')
    file_handler = handlers.RotatingFileHandler(log_path, maxBytes=5*1024*1024)
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
