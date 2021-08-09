
import pickle
from astropy.io import ascii


def dataSave(data, file_out, wm='wb', fdict={}):
    """
    Create output data file with: MP values / MCMC samples / mass values.
    """
    if wm == 'wb':
        with open(file_out, wm) as f:
            pickle.dump(data, f)
    elif wm == 'w':
        ascii.write(data, file_out, overwrite=True, formats=fdict)


def dataRead(clust_name, file_in, rm='rb'):
    """
    Read data file with: MP values / MCMC samples / mass values.
    """
    if rm == 'rb':
        file_in = file_in.replace('output', 'input')
        file_in = file_in.replace(clust_name + '/', '')
        with open(file_in, rm) as f:
            data = pickle.load(f)
    elif rm == 'r':
        data = ascii.read(file_in)

    return data
