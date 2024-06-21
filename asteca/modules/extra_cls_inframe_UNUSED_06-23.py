import numpy as np
from . import cluster_priv as cp

# :param centers_ex: List of dictionaries, one dictionary for each object that
#     shares the frame with the cluster and that should be ignored. The
#     dictionaries must have at most three keys, 'xy', 'pms', 'plx', each with a
#     list containing the center value(s) in those dimensions. Defaults to
#     ``None``. Examples:
#     one object:
#     ``[{'xy': [105.39, 0.9], 'plx': [1.3]}]``,
#     two objects:
#     ``[{'xy': [105.39, 0.9], 'pms': [3.5, -0.7]},
#     {'xy': [0.82, -4.5], 'pms': [3.5, -0.7], 'plx': [3.5]}]``
# :type centers_ex: list[dict] | None


def get_extra_cls(centers_ex):
    """ """

    # Prepare dictionary of parameters for extra clusters in frame (if any)
    extra_cls_dict, dims_msk = prep_extra_cl_dict(centers_ex)


def prep_extra_cl_dict(centers_ex):
    """
    The parameter 'centers_ex' must be a list of dictionaries, one
    dictionary for each extra cluster in the frame. The dictionaries
    must have at most three keys, 'xy', 'pms', 'plx', each with a list
    containing the center value(s) in those dimensions.

    Examples:

    centers_ex = [{'xy': [105.39, 0.9], 'plx': [1.3]}]

    centers_ex = [
        {'xy': [105.39, 0.9], 'pms': [3.5, -0.7]},
        {'xy': [0.82, -4.5], 'pms': [3.5, -0.7], 'plx': [3.5]}
    ]
    """
    extra_cls_dict = {"run_flag": False}

    if centers_ex is None:
        return extra_cls_dict, None

    extra_cls_dict["run_flag"] = True

    # Set the distribution of dimensions
    dims_msk = {
        "xy": np.array([0, 1]),
        "pms": np.array([2, 3]),
        "plx": np.array([4]),
        "xy_pms": np.array([0, 1, 2, 3]),
        "xy_plx": np.array([0, 1, 4]),
        "pms_plx": np.array([2, 3, 4]),
        "xy_pms_plx": np.array([0, 1, 2, 3, 4]),
    }

    # Read centers of the extra clusters in frame
    dims = ["xy", "pms", "plx"]
    dim_keys, cents = [], []
    for extra_cl in centers_ex:
        dims_ex = extra_cl.keys()
        dims_id, centers = "", []
        for d in dims:
            if d in dims_ex:
                center = extra_cl[d]
                if center is not None:
                    dims_id += "_" + d
                    centers += center
        # Store dimensions ids and centers for each extra cluster
        dim_keys.append(dims_id[1:])
        cents.append(centers)

    extra_cls_dict["dim_keys"] = dim_keys
    extra_cls_dict["centers"] = cents

    return extra_cls_dict, dims_msk


def filter_cls_in_frame(
    lon,
    lat,
    pmRA,
    pmDE,
    plx,
    xy_c,
    vpd_c,
    plx_c,
    idx_survived,
    extra_cls_dict,
    dims_msk,
    N_clust_min,
    dims_norm,
):
    """
    Filter extra clusters in frame (if any)
    """
    # If there are no extra clusters to remove, skip
    if extra_cls_dict["run_flag"] is False:
        return idx_survived

    # Add centers to the end of these lists to normalize their values too
    dims_plus_cents = [list(_) for _ in (lon, lat, pmRA, pmDE, plx)]
    for cent in extra_cls_dict["centers"]:
        # For each dimension, some won't be present
        for i in range(5):
            try:
                dims_plus_cents[i].append(cent[i])
            except IndexError:
                dims_plus_cents[i].append(np.nan)
    lon2, lat2, pmRA2, pmDE2, plx2 = dims_plus_cents

    # Data normalization
    msk = np.full(len(lon2), True)
    data, cents = get_dims_norm(
        N_clust_min, dims_norm, lon2, lat2, pmRA2, pmDE2, plx2, xy_c, vpd_c, plx_c, msk
    )
    data = data.T

    # Extract normalized centers for extra clusters
    new_cents = []
    for i in range(len(extra_cls_dict["centers"]), 0, -1):
        new_cents.append([_ for _ in data[:, -i] if not np.isnan(_)])
    # Remove center values from the 'data' array
    data = data[:, : -len(new_cents)]

    # Find the distances to the center for all the combinations of data
    # dimensions in the extra clusters in frame
    dists = {}
    for dims in list(set(extra_cls_dict["dim_keys"])):
        msk = dims_msk[dims]
        dists[dims] = cp.get_Nd_dists(cents[:, msk], data[msk, :].T, True)

    # Identify stars closer to the cluster's center than the centers of
    # extra clusters
    msk_d = np.full(len(idx_survived), True)
    for i, cent_ex in enumerate(new_cents):
        dims_ex = extra_cls_dict["dim_keys"][i]
        msk = dims_msk[dims_ex]
        # Distance to this extra cluster's center
        dists_ex = cp.get_Nd_dists(np.array([cent_ex]), data[msk, :].T, True)

        # If the distance to the selected cluster center is smaller than
        # the distance to this extra cluster's center, keep the star
        msk_d = msk_d & (dists[dims_ex] <= dists_ex)

    # Never return less than N_clust_min stars
    if msk_d.sum() < N_clust_min:
        return idx_survived

    return idx_survived[msk_d]
