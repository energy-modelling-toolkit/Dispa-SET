import logging
import sys


def shrink_to_64(x, N=64):
    """
    Function that reduces the length of the keys to be written to 64 (max admissible length for GAMS)

    :param x:   String or list of strings
    :param N:   Integer with the maximum string length (if different from 64)

    :returns:   Shrinked string or list of strings
    """

    def shrink_singlestring(key, N):
        if len(key) >= N:
            return key[:20] + ' ... ' + key[-20:]
        else:
            return key

    if type(x) == str:
        return shrink_singlestring(x, N)
    elif type(x) == list:
        return [shrink_singlestring(xx, N) for xx in x]
    else:
        logging.critical('Argument type not supported')
        sys.exit(1)


def clean_strings(x, exclude_digits=False, exclude_punctuation=False):
    """
    Function to convert strange unicode
    and remove characters punctuation

    :param x: any string or list of strings

    Usage::

        df['DisplayName'].apply(clean_strings)

    """
    if sys.version_info >= (3,0): # Skip this funcion if python version is >3. Have to test better TODO
        return x
    import unicodedata
    import string
    def clean_singlestring(x):
        if exclude_digits:  # modify the following depending on what you need to exclude
            exclude1 = set(string.punctuation)
            # exception to the exclusion:
            exclude1.remove('_')
            exclude1.remove('-')
            exclude1.remove('[')
            exclude1.remove(']')
        else:
            exclude1 = set([])
        if exclude_punctuation:
            exclude2 = set(string.digits)
        else:
            exclude2 = set([])
        exclude = exclude1 | exclude2

        # http://stackoverflow.com/questions/2365411/python-convert-unicode-to-ascii-without-errors
        x = str(x).decode('utf-8')  # to string byte and then unicode
        x = unicodedata.normalize('NFKD', x).encode('ascii', 'ignore')  # convert utf characters and to ascii

        # x = x.upper() #to UPPERCASE
        x = ''.join(ch for ch in x if ch not in exclude)  # remove numbers and punctuation
        return x
    if isinstance(x, str):
        return clean_singlestring(x)
    elif isinstance(x, list):
        return [clean_singlestring(xx) for xx in x]
    else:
        logging.error('Argument type not supported')
        sys.exit(1)

def force_str(x):
    """ Used to get a str object both in python 2 and 3 although they represent different objects (byte vs unicode)
    It is small hack for py2->3 compatibility of gams APIs which require a str object
    """
    if isinstance(x, str):
        return x
    elif isinstance(x,bytes):
        return str(x.decode('utf-8'))
    else:
        return x.encode()