# -*- coding: utf-8 -*-
"""
Set of functions used to test the availability of different libraries required by Dispa-SET

@author: Sylvain Quoilin
"""
import platform


def check_logging():
    import logging
    import logging.config
    # select ColorStreamHandler based on platform
    
    _LOGCONFIG = {
         "version": 1,
         "disable_existing_loggers": False,
         'formatters': {
            'standard': {
                'format': '%(asctime)s [%(levelname)-8s] (%(funcName)s): %(message)s',
                'datefmt': '%y/%m/%d %H:%M:%S'
            },
            'notime': {
                'format': '[%(levelname)-8s] (%(funcName)s): %(message)s',
                'datefmt': '%y/%m/%d %H:%M:%S'
            },                
         },
         "handlers": {
             "console": {
                "class": "logging.StreamHandler",
                 "stream": "ext://sys.stderr",
    #             "stream": "sys.stdout",
                 "level": "INFO",
                 'formatter': 'notime',
             },
    
             "error_file": {
                 "class": "logging.FileHandler",
                 "level": "WARNING",
                 'formatter': 'standard',
                 'filename': 'warn.log',
                 'encoding': 'utf8'
    
             }
         },
    
         "root": {
             "level": "INFO",
             "handlers": ["console", "error_file"],
         }
    }
    
    try:
        logging.config.dictConfig(_LOGCONFIG)
            # Display a few messages:
        logging.warn('This is a warning')
        logging.info('This is information')
        logging.critical('This is critical')
        logging.error('This is an error')
        return True
    except Exception, e:
        print('Logging cannot be handled. The following error msg was generated:')
        print(e)
        return False