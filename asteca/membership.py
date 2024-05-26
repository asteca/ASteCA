import numpy as np
from .modules import cluster_priv as cp
from .modules.fastmp import fastMP


def fastmp(
    cluster,
    xy_c=None,
    vpd_c=None,
    plx_c=None,
    fixed_centers=False,
    N_cluster=None,
    N_clust_min=25,
    N_clust_max=5000,
    centers_ex=None,
    N_resample=1000,
):
    """Assign membership probabilities.

    Estimate the probability of being a true cluster member for all observed
    stars using the fastMP algorithm. The algorithm was described in detail in
    the `article <https://academic.oup.com/mnras/article/526/3/4107/7276628>`__
    were we introduced the
    `Unified Cluster Catalogue (UCC) <https://ucc.ar/>`__.

    :param cluster: :py:class:`asteca.cluster` object with the loaded data for the
        observed cluster
    :type cluster: :py:class:`asteca.cluster`
    :param xy_c: Estimated value for the (RA, DEC) center, defaults to ``None``
    :type xy_c: tuple(float, float), optional
    :param vpd_c: Estimated value for the (pmRA, pmDE) center, defaults to ``None``
    :type vpd_c: tuple(float, float), optional
    :param plx_c: Estimated value for the plx center, defaults to ``None``
    :type plx_c: float, optional
    :param fixed_centers: If ``True`` any center estimate (xy_c, vpd_c, plx_c)
        given will be kept fixed throughout the process, defaults to ``False``
    :type fixed_centers: bool
    :param N_cluster: Estimated number of members, defaults to ``None``
    :type N_cluster: int, optional
    :param N_clust_min: Minimum number of cluster members, defaults to ``25``
    :type N_clust_min: int, optional
    :param N_clust_max: Maximum number of cluster members, defaults to ``5000``
    :type N_clust_max: int, optional
    :param centers_ex: List of dictionaries, one dictionary for each object that
        shares the frame with the cluster and that should be ignored. The
        dictionaries must have at most three keys, 'xy', 'pms', 'plx', each with a
        list containing the center value(s) in those dimensions. Defaults to
        ``None``. Examples:
        one object:
        ``[{'xy': [105.39, 0.9], 'plx': [1.3]}]``,
        two objects:
        ``[{'xy': [105.39, 0.9], 'pms': [3.5, -0.7]},
        {'xy': [0.82, -4.5], 'pms': [3.5, -0.7], 'plx': [3.5]}]``
    :type centers_ex: list(dict), optional
    :param N_resample: Maximum number of resamples, defaults to ``1000``
    :type N_resample: int, optional

    :return: Membership probabilities for all stars in the frame
    :rtype: np.array
    """
    if any(
        [
            _ is None
            for _ in (
                cluster.pmra_v,
                cluster.pmde_v,
                cluster.plx_v,
                cluster.e_pmra_v,
                cluster.e_pmde_v,
                cluster.e_plx_v,
            )
        ]
    ):
        raise ValueError(
            "fastMP requires (ra, dec, pmra, pmde, plx) data and\n"
            + "their uncertainties to be present in the 'cluster' object"
        )
    for k, N in {
        "N_cluster": N_cluster,
        "N_clust_min": N_clust_min,
        "N_clust_max": N_clust_max,
        "N_resample": N_resample,
    }.items():
        if N is not None and not isinstance(N, int):
            raise ValueError(f"{k}={N}, must be either 'None' or an 'int'")
    if fixed_centers is True and xy_c is None and vpd_c is None and plx_c is None:
        raise ValueError("fixed_centers=True but no center values were given")

    # Convert (RA, DEC) to (lon, lat)
    glon, glat = cp.radec2lonlat(cluster.ra_v.values, cluster.dec_v.values)
    if xy_c is not None:
        xy_c = cp.radec2lonlat(xy_c[0], xy_c[1])

    # Generate input data array for fastMP
    X = np.array(
        [
            glon,
            glat,
            cluster.pmra_v,
            cluster.pmde_v,
            cluster.plx_v,
            cluster.e_pmra_v,
            cluster.e_pmde_v,
            cluster.e_plx_v,
        ]
    )

    print("\nRunning fastMP...")
    print(f"fixed_centers : {fixed_centers}")
    print(f"N_cluster     : {N_cluster}")
    print(f"N_clust_min   : {N_clust_min}")
    print(f"N_clust_max   : {N_clust_max}")
    centers_ex_flag = False
    if centers_ex is not None:
        centers_ex_flag = True
    print(f"centers_ex    : {centers_ex_flag}")
    print(f"N_resample    : {N_resample}")
    probs = fastMP(
        X,
        xy_c,
        vpd_c,
        plx_c,
        fixed_centers,
        N_cluster,
        N_clust_min,
        N_clust_max,
        centers_ex,
        N_resample,
    )

    return probs
