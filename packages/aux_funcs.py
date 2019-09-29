
import numpy as np
from scipy import stats
from astropy.stats import sigma_clipped_stats


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
