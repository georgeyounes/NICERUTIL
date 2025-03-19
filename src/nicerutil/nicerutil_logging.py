import logging


def get_logger(name):
    logging.basicConfig(level=logging.INFO)

    logger = logging.getLogger(name)
    logFormatter = logging.Formatter('[%(asctime)s] %(levelname)8s %(message)s ' +
                                     '(%(filename)s:%(lineno)s)', datefmt='%Y-%m-%d %H:%M:%S')

    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    consoleHandler.setLevel(logging.INFO)

    fileHandler = logging.FileHandler('nicerutil.log')
    fileHandler.setFormatter(logFormatter)
    fileHandler.setLevel(logging.INFO)

    logger.addHandler(consoleHandler)
    logger.addHandler(fileHandler)

    logger.propagate = False

    return logger
