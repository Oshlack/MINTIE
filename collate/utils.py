import os
import sys
import logging
import pickle

PROGRAM_NAME = 'MAKE_SUPERTRANSCRIPT'
def exit_with_error(message, exit_status):
    '''
    Print an error message to stderr, prefixed by the program name and 'ERROR'.
    Then exit program with supplied exit status.
    Arguments:
        message: an error message as a string.
        exit_status: a positive integer representing the exit status of the
            program.
    '''
    logging.error(message)
    print("{} ERROR: {}, exiting".format(PROGRAM_NAME, message), file=sys.stderr)
    sys.exit(exit_status)

def init_logging(log_filename):
    '''
    If the log_filename is defined, then
    initialise the logging facility, and write log statement
    indicating the program has started, and also write out the
    command line from sys.argv
    Arguments:
        log_filename: either None, if logging is not required, or the
            string name of the log file to write to
    Result:
        None
    '''
    if log_filename is not None:
        logging.basicConfig(filename=log_filename,
                            level=logging.DEBUG,
                            filemode='w',
                            format='%(asctime)s %(levelname)s - %(message)s',
                            datefmt='%m-%d-%Y %H:%M:%S')
        logging.info('program started')
    logging.info('command line: %s', ' '.join(sys.argv))

def cached(cachefile):
    '''
    source: https://datascience.blog.wzb.eu/2016/08/12/a-tip-for-the-impatient-simple-caching-with-python-pickle-and-decorators/
    A function that creates a decorator which will use "cachefile"
    for caching the results of the decorated function "fn".
    '''
    def decorator(fn):  # define a decorator for a function "fn"
        def wrapped(*args, **kwargs):   # define a wrapper that will finally call "fn" with all arguments
          # if cache exists -> load it and return its content
          if os.path.exists(cachefile):
              with open(cachefile, 'rb') as cachehandle:
                logging.info("using cached result from '%s'" % cachefile)
                return pickle.load(cachehandle)

          # execute the function with all arguments passed
          res = fn(*args, **kwargs)

          # write to cache file
          with open(cachefile, 'wb') as cachehandle:
            logging.info("saving result to cache '%s'" % cachefile)
            pickle.dump(res, cachehandle)

          return res

        return wrapped

    return decorator   # return this "customized" decorator that uses "cachefile"
