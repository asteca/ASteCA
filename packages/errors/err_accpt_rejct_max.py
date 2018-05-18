
import numpy as np
import warnings


def main(cld, err_max):
    """
    Accept stars with photometric errors < e_max in *all* its mags and colors.

    All 'nan' values are kept, i.e.: evaluated to True.
    Source: https://stackoverflow.com/a/48584644/1391441

    This means that there will be a different number of stars with valid data
    in each photometric dimension.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        m_msk = np.logical_or(cld['em'] < err_max, np.isnan(cld['em']))
        c_msk = np.logical_or(cld['ec'] < err_max, np.isnan(cld['ec']))

    acpt_indx = np.flatnonzero((m_msk.all(0) & c_msk.all(0))).tolist()
    rjct_indx = np.setdiff1d(
        np.arange(len(cld['em'][0])), acpt_indx).tolist()

    return acpt_indx, rjct_indx
