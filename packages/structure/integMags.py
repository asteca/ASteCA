
import numpy as np
from scipy import spatial


def main(clp, x, y, mags, N_dists=50, **kwargs):
    """
    Calculate integrated magnitude versus distance to center.

    HARDCODED
    N_dists : number of distance steps.
    """
    xy = np.array([x, y]).T
    rdp_dists_cent = spatial.distance.cdist([clp['kde_cent']], xy)[0]

    dists = np.linspace(
        np.percentile(rdp_dists_cent, 1), np.percentile(rdp_dists_cent, 99),
        N_dists)

    integ_dists, integ_mags = [], []
    # d0 = 0.
    for d1 in dists:
        msk = (rdp_dists_cent < d1)  # (rdp_dists_cent >= d0) &
        m = mags[0][msk]
        # Remove nans
        m = m[~np.isnan(m)]
        if len(m) > 0:
            integ_dists.append(d1)
            integ_mags.append(-2.5 * np.log10(np.sum(10 ** (m / -2.5))))
        # d0 = d1

    clp['integ_dists'], clp['integ_mags'] = integ_dists, integ_mags
    return clp
