
import logging
import sys

g_logger=None

def init_logger():
    global g_logger
    g_logger=logging.getLogger('NumBAT log')

def report(msg):
    if g_logger is None:
        init_logger()
    g_logger.warning(msg)

def report_and_exit(msg):
    if g_logger is None:
        init_logger()
    g_logger.error('\nFatal error: \n %s', msg)
#g_logger.error('from logger')
    sys.exit(1)
