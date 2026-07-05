import sys
import logging
import logging.handlers as handlers

# logging.info('started')  # Most simple

formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

logger = logging.getLogger('main')
logger.setLevel(logging.INFO)
logger.propagate = False  # third-party libs (e.g. jpype/orekit) may configure the
                          # root logger; without this every line would be duplicated

logHandler = handlers.RotatingFileHandler('../output/main.log', maxBytes=5*1024*1024)  # log to the main.log file
logHandler.setLevel(logging.INFO)
logHandler.setFormatter(formatter)
logger.addHandler(logHandler)

ch = logging.StreamHandler(sys.stdout)  # log to the stdout
ch.setFormatter(formatter)
logger.addHandler(ch)