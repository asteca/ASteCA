
import numpy as np
from scipy import stats


def reject_outliers(data, m=4.):
    """
    Outlier rejection.
    Source: https://stackoverflow.com/a/16562028/1391441
    """
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d / (mdev if mdev else 1.)
    return data[s < m]


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
