
import numpy as np
import warnings


def main(cld, err_max):
    """
    Accept stars with photometric errors < e_max in mag and color.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        acpt_indx = np.flatnonzero(
            ((cld['em'] < err_max).all(0) &
                (cld['ec'] < err_max).all(0))).tolist()
        rjct_indx = np.setdiff1d(
            np.arange(len(cld['em'][0])), acpt_indx).tolist()

    return acpt_indx, rjct_indx
