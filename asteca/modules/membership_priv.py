import numpy as np


def get_dims_norm(
    N_clust_min, dims_norm, lon, lat, pmRA, pmDE, plx, xy_c, vpd_c, plx_c, msk
):
    """
    Normalize dimensions using twice the median of the selected probable
    members.
    """
    data_5d = np.array([lon, lat, pmRA, pmDE, plx]).T
    cents_5d = np.array([xy_c + vpd_c + [plx_c]])

    if msk is None or len(msk) < N_clust_min:
        return data_5d, cents_5d

    data_mvd = data_5d - cents_5d

    if dims_norm is None:
        dims_norm = 2 * np.nanmedian(abs(data_mvd[msk]), 0)
        data_norm = data_mvd / dims_norm
        # from sklearn.preprocessing import RobustScaler
        # data_norm = RobustScaler().fit(data_5d).transform(data_5d)
    else:
        data_norm = data_mvd / dims_norm

    cents_norm = np.array([[0.0, 0.0, 0.0, 0.0, 0.0]])

    return data_norm, cents_norm
