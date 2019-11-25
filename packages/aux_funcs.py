
from collections import defaultdict, Iterable
import numpy as np
from scipy import stats
from astropy.stats import sigma_clipped_stats


def flatten(l):
    """
    Source: https://stackoverflow.com/a/2158532/1391441
    """
    for el in l:
        if isinstance(el, Iterable) and not\
                isinstance(el, (str, bytes)):
            for sub in flatten(el):
                yield sub
        else:
            yield el


def list_duplicates(seq):
    """
    Find and report duplicates in list.

    Source: https://stackoverflow.com/a/5419576/1391441
    """
    tally = defaultdict(list)
    for i, item in enumerate(seq):
        tally[item].append(i)
    dups = ((key, map(str, locs)) for key, locs in tally.items()
            if len(locs) > 1)
    return dups


def reject_outliers(data, m=4.):
    """
    Outlier rejection.
    """
    mean, median, std = sigma_clipped_stats(data)
    msk = (data > mean - m * std) & (data < mean + m * std)
    return data[msk]


def kde1D(data):
    """
    1D KDE.
    """
    # Cluster region KDE curve.
    xmin, xmax = np.min(data), np.max(data)
    # Define KDE limits.
    x_rang = .1 * (xmax - xmin)
    kde_x = np.mgrid[xmin - x_rang:xmax + x_rang:1000j]
    kernel_cl = stats.gaussian_kde(data)
    # KDE for plotting.
    kde = np.reshape(kernel_cl(kde_x).T, kde_x.shape)

    return kde_x, kde


def circFrac(cent, rad, x0, x1, y0, y1, N_tot=250000):
    """
    Use Monte Carlo to estimate the fraction of the area of a circle centered
    in (cx, cy) with a radius of 'rad', that is located within the frame given
    by the limits 'x0, x1, y0, y1'.
    """
    cx, cy = cent
    # Source: https://stackoverflow.com/a/50746409/1391441
    r = rad * np.sqrt(np.random.uniform(0., 1., N_tot))
    theta = np.random.uniform(0., 1., N_tot) * 2 * np.pi
    xr = cx + r * np.cos(theta)
    yr = cy + r * np.sin(theta)

    # Points within the circle that are within the frame.
    msk_xy = (xr > x0) & (xr < x1) & (yr > y0) & (yr < y1)

    # The area is the points within circle and frame over the points within
    # circle.
    return msk_xy.sum() / N_tot
