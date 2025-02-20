import logging
import logging.handlers


def get_logger(name):
    logFormatter = logging.Formatter('[%(asctime)s] %(levelname)8s %(message)s ' +
                                     '(%(filename)s:%(lineno)s)', datefmt='%Y-%m-%d %H:%M:%S')
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)

    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    consoleHandler.setLevel(logging.INFO)
    logger.addHandler(consoleHandler)

    fileHandler = logging.handlers.RotatingFileHandler('nicerutil.log', mode='w')
    fileHandler.setLevel(logging.DEBUG)
    fileHandler.setFormatter(logFormatter)
    logger.addHandler(fileHandler)

    return logger
