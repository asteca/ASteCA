
from collections import defaultdict
from collections.abc import Iterable
import numpy as np
from scipy import stats
from astropy.stats import sigma_clipped_stats


def flatten(lst):
    """
    Source: https://stackoverflow.com/a/2158532/1391441
    """
    for el in lst:
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


def kde1D(data, xmin=None, xmax=None, bw=None, xr_fr=.1):
    """
    1D KDE.
    """
    if xmin is None and xmax is None:
        xmin, xmax = np.min(data), np.max(data)
    # Define KDE limits.
    x_rang = xr_fr * (xmax - xmin)
    kde_x = np.mgrid[xmin - x_rang:xmax + x_rang:1000j]
    try:
        kernel_cl = stats.gaussian_kde(data, bw_method=bw)
        # KDE for plotting.
        kde = np.reshape(kernel_cl(kde_x).T, kde_x.shape)
    except np.linalg.LinAlgError:
        kde = np.array([])

    return kde_x, kde


def circFrac(cent, rad, x0, x1, y0, y1, rand_01_MC, cos_t, sin_t):
    """
    Use Monte Carlo to estimate the fraction of the area of a circle centered
    in (cx, cy) with a radius of 'rad', that is located within the frame given
    by the limits 'x0, x1, y0, y1'.
    """
    N_tot = len(rand_01_MC)

    cx, cy = cent
    # Source: https://stackoverflow.com/a/50746409/1391441
    # r = rad * np.sqrt(np.random.uniform(0., 1., N_tot))
    # theta = np.random.uniform(0., 1., N_tot) * 2 * np.pi
    rand_01_MC = rand_01_MC * rad
    xr = cx + rand_01_MC * cos_t
    yr = cy + rand_01_MC * sin_t

    # Points within the circle that are within the frame.
    msk_xy = (xr > x0) & (xr < x1) & (yr > y0) & (yr < y1)

    # The area is the points within circle and frame over the points within
    # circle.
    return msk_xy.sum() / N_tot


def ellipFrac(
        cent, a, theta, ell, x0, x1, y0, y1, rand_01_MC, cos_t, sin_t):
    """
    Use Monte Carlo to estimate the fraction of the area of an ellipse centered
    in (cx, cy) and rotated an angle theta, with a semi-major axis 'a' and
    ellipticity 'ell', that is located within the frame given by the limits
    'x0, x1, y0, y1'.
    """
    N_tot = len(rand_01_MC)

    cx, cy = cent
    b = a * (1 - ell)
    cos_th = np.cos(theta)
    sin_th = np.sin(theta)

    # r = rad * np.sqrt(np.random.uniform(0., 1., N_tot))
    # theta = np.random.uniform(0., 1., N_tot) * 2 * np.pi
    rand_01_MC_a = rand_01_MC * a
    rand_01_MC_b = rand_01_MC * b

    # https://math.stackexchange.com/questions/2645689/what-is-the-parametric-equation-of-a-rotated-ellipse-given-the-angle-of-rotatio

    xr = cx + rand_01_MC_a * cos_t * cos_th - rand_01_MC_b * sin_t * sin_th
    yr = cy + rand_01_MC_a * cos_t * sin_th + rand_01_MC_b * sin_t * cos_th

    # Points within the ellipse that are within the frame.
    msk_xy = (xr > x0) & (xr < x1) & (yr > y0) & (yr < y1)

    return msk_xy.sum() / N_tot


def monteCarloPars(N_MC=10000):
    """
    Parameters for the Monte Carlo estimation of the fraction of circle/ellipse
    area within a rectangle. Used by the circFrac(), ellipFrac() functions.

    N_MC: number of Monte Carlo points to use when estimating circular areas.
    """

    # Obtain these values here so they are not estimated every time the MC
    # in circFrac() is called.
    rand_01_MC = np.sqrt(np.random.uniform(0., 1., N_MC))
    theta = np.random.uniform(0., 1., N_MC) * 2 * np.pi
    cos_t, sin_t = np.cos(theta), np.sin(theta)

    return rand_01_MC, cos_t, sin_t
