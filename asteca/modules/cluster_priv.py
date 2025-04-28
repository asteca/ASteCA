import warnings

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.modeling.models import KingProjectedAnalytic1D
from scipy import spatial, stats
from scipy.optimize import least_squares


def radec2lonlat(ra: float | np.ndarray, dec: float | np.ndarray) -> np.ndarray:
    """Convert from right ascension and declination to galactic longitude and latitude.

    :param ra: Right ascension.
    :type ra: float | np.ndarray
    :param dec: Declination.
    :type dec: float | np.ndarray

    :return: Galactic longitude and latitude.
    :rtype: np.ndarray
    """
    gc = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)  # pyright: ignore
    lb = gc.transform_to("galactic")
    return np.array([lb.l.value, lb.b.value])  # pyright: ignore


def lonlat2radec(lon: float | np.ndarray, lat: float | np.ndarray) -> np.ndarray:
    """Convert from galactic longitude and latitude to right ascension and declination.

    :param lon: Galactic longitude.
    :type lon: float | np.ndarray
    :param lat: Galactic latitude.
    :type lat: float | np.ndarray

    :return: Right ascension and declination.
    :rtype: np.ndarray
    """
    gc = SkyCoord(l=lon * u.degree, b=lat * u.degree, frame="galactic")  # pyright: ignore
    ra, dec = gc.fk5.ra.value, gc.fk5.dec.value  # pyright: ignore
    return np.array([ra, dec])


def reject_nans(arr_data: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Remove nans in arr_data

    :param arr_data: Array of data.
    :type arr_data: np.ndarray

    :return: Indexes of non-nan data and the non-nan arr_data.
    :rtype: tuple[np.ndarray, np.ndarray]
    """
    msk_all = []
    # Process each dimension separately
    for arr in arr_data:
        # Identify non-nan arr_data
        msk = ~np.isnan(arr)
        # Keep this non-nan arr_data
        msk_all.append(msk.data)
    # Combine into a single mask
    msk_accpt = np.logical_and.reduce(msk_all)

    # Indexes that survived
    idx_clean = np.arange(arr_data.shape[1])[msk_accpt]

    return idx_clean, arr_data.T[msk_accpt].T


def get_knn_5D_center(
    lon: np.ndarray,
    lat: np.ndarray,
    pmRA: np.ndarray,
    pmDE: np.ndarray,
    plx: np.ndarray,
    xy_c: tuple[float, float] | None,
    vpd_c: tuple[float, float] | None,
    plx_c: float | None,
    N_clust_min: int,
    N_clust_max: int,
) -> tuple[float, float, float, float, float]:
    """Estimate the 5-dimensional center of a cluster.

    1. Keep only 'N_cent' stars if xy_c or plx_c are given
    2. (Re)Estimate the center in PMs (the value can be given as input)
    3. Obtain the 'N_cent' stars closest to the available center values
    4. Estimate the 5-dimensional final center using kNN

    :param lon: Galactic Longitude.
    :type lon: np.ndarray
    :param lat: Galactic Latitude.
    :type lat: np.ndarray
    :param pmRA: Proper motion in Right Ascension.
    :type pmRA: np.ndarray
    :param pmDE: Proper motion in Declination.
    :type pmDE: np.ndarray
    :param plx: Parallax.
    :type plx: np.ndarray
    :param xy_c: Center coordinates in (lon, lat).
    :type xy_c: tuple[float, float] | None
    :param vpd_c: Center coordinates in proper motions (pmRA, pmDE).
    :type vpd_c: tuple[float, float] | None
    :param plx_c: Center coordinate in parallax.
    :type plx_c: float | None
    :param N_clust_min: Minimum number of stars in the cluster.
    :type N_clust_min: int
    :param N_clust_max: Maximum number of stars in the cluster.
    :type N_clust_max: int

    :return: Center coordinates in (lon, lat, pmRA, pmDE, plx).
    :rtype: tuple[float, float, float, float, float]
    """
    N_tot = len(lon)
    # N_clust_min < N_cent < 250
    N_cent = max(N_clust_min, min(250, int(0.1 * N_tot)))

    # Get filtered stars close to given xy+Plx centers (if given) to use
    # in the PMs center estimation.
    if xy_c is None and plx_c is None:
        # If neither xy_c nor plx_c were given and the number of stars in the frame
        # is larger than N_clust_max, select twice the N_clust_max stars closest to the
        # center of the XY frame (i.e.: this assumes that the cluster is centered in
        # XY). This prevents very large frames to deviate from the actual cluster's
        # center because most stars in the VPD are distributed around another center
        # value
        if N_tot > N_clust_max:
            data = np.array([lon, lat]).T
            cent = np.array([np.median(data, 0)])
            idx = get_Nd_dists(cent, data)[: 2 * N_clust_max]
            pmRA_i, pmDE_i = pmRA[idx], pmDE[idx]
        else:
            pmRA_i, pmDE_i = np.array(pmRA), np.array(pmDE)
    else:
        pmRA_i, pmDE_i = filter_pms_stars(
            xy_c, plx_c, lon, lat, pmRA, pmDE, plx, N_cent
        )

    # (Re)estimate VPD center
    vpd_c_i = get_pms_center(vpd_c, N_clust_min, pmRA_i, pmDE_i)

    # Get N_cent stars closest to vpd_c and given xy and/or plx centers
    lon_i, lat_i, pmRA_i, pmDE_i, plx_i = get_stars_close_center(
        lon, lat, pmRA, pmDE, plx, xy_c, vpd_c_i, plx_c, N_cent
    )

    # kNN center
    x_c, y_c, pmra_c, pmde_c, plx_c = get_kNN_center(
        N_clust_min, np.array([lon_i, lat_i, pmRA_i, pmDE_i, plx_i]).T
    )

    return x_c, y_c, pmra_c, pmde_c, plx_c


def get_Nd_dists(cents: np.ndarray, data: np.ndarray) -> np.ndarray:
    """Obtain indexes and distances of stars to the given center

    :param cents: Center coordinates.
    :type cents: np.ndarray
    :param data: Array of data.
    :type data: np.ndarray

    :return: Indexes or distances of stars to the given center.
    :rtype: np.ndarray
    """
    # Distances to center
    dist_Nd = spatial.distance.cdist(data, cents).T[0]

    # Indexes that sort the distances
    d_idxs = dist_Nd.argsort()

    return d_idxs


def filter_pms_stars(
    xy_c: tuple[float, float] | None,
    plx_c: float | None,
    lon: np.ndarray,
    lat: np.ndarray,
    pmRA: np.ndarray,
    pmDE: np.ndarray,
    plx: np.ndarray,
    N_cent: int,
) -> tuple[np.ndarray, np.ndarray]:
    """If either xy_c or plx_c values are given, select the 'N_cent' stars
    closest to this 1D/2D/3D center, and return their proper motions.

    :param xy_c: Center coordinates in (lon, lat).
    :type xy_c: tuple[float, float] | None
    :param plx_c: Center coordinate in parallax.
    :type plx_c: float | None
    :param lon: Galactic Longitude.
    :type lon: np.ndarray
    :param lat: Galactic Latitude.
    :type lat: np.ndarray
    :param pmRA: Proper motion in Right Ascension.
    :type pmRA: np.ndarray
    :param pmDE: Proper motion in Declination.
    :type pmDE: np.ndarray
    :param plx: Parallax.
    :type plx: np.ndarray
    :param N_cent: Number of stars to select.
    :type N_cent: int

    :raises ValueError: If xy_c and plx_c are both None.

    :return: Proper motions of the selected stars.
    :rtype: tuple[np.ndarray, np.ndarray]
    """
    # Create arrays with required shape
    if xy_c is None and plx_c is not None:
        cent = np.array([[plx_c]])
        data = np.array([plx]).T
    elif xy_c is not None and plx_c is None:
        cent = np.array([xy_c])
        data = np.array([lon, lat]).T
    elif xy_c is not None and plx_c is not None:
        cent = np.array([list(xy_c) + [plx_c]])
        data = np.array([lon, lat, plx]).T
    else:
        raise ValueError("Either xy_c or plx_c must be given")

    # Closest stars to the selected center
    idx = get_Nd_dists(cent, data)[:N_cent]
    pmRA_i, pmDE_i = pmRA[idx], pmDE[idx]

    return pmRA_i, pmDE_i


def get_pms_center(
    vpd_c: tuple[float, float] | None,
    N_clust_min: int,
    pmRA: np.ndarray,
    pmDE: np.ndarray,
    N_bins: int = 50,
    zoom_f: int = 4,
    N_zoom: int = 10,
) -> tuple[float, float]:
    """Estimate the center in proper motion space.

    :param vpd_c: Center coordinates in proper motions (pmRA, pmDE).
    :type vpd_c: tuple[float, float] | None
    :param N_clust_min: Minimum number of stars in the cluster.
    :type N_clust_min: int
    :param pmRA: Proper motion in Right Ascension.
    :type pmRA: np.ndarray
    :param pmDE: Proper motion in Declination.
    :type pmDE: np.ndarray
    :param N_bins: Number of bins for the 2D histogram, defaults to 50
    :type N_bins: int
    :param zoom_f: Zoom factor for the iterative center estimation, defaults to 4
    :type zoom_f: int
    :param N_zoom: Number of zoom iterations, defaults to 10
    :type N_zoom: int

    :raises ValueError: If the PMs center could not be estimated

    :return: Center coordinates in proper motions (pmRA, pmDE).
    :rtype: tuple[float, float]
    """
    vpd = np.array([pmRA, pmDE]).T

    # Center in PMs space
    cxym = None
    for _ in range(N_zoom):
        N_stars = vpd.shape[0]
        if N_stars < N_clust_min:
            break

        # Find center coordinates as max density
        x, y = vpd.T
        H, edgx, edgy = np.histogram2d(x, y, bins=N_bins)
        flat_idx = H.argmax()
        cbx, cby = np.unravel_index(flat_idx, H.shape)
        cx = (edgx[cbx + 1] + edgx[cbx]) / 2.0
        cy = (edgy[cby + 1] + edgy[cby]) / 2.0
        # Store the auto center values
        cxym = (cx, cy)

        # If a manual center was set, use these values to zoom in
        if vpd_c is not None:
            cx, cy = vpd_c

        # Zoom in
        rx, ry = edgx[1] - edgx[0], edgy[1] - edgy[0]
        msk = (
            (x < (cx + zoom_f * rx))
            & (x > (cx - zoom_f * rx))
            & (y < (cy + zoom_f * ry))
            & (y > (cy - zoom_f * ry))
        )
        vpd = vpd[msk]

    if cxym is not None:
        cx, cy = cxym
    else:
        if vpd_c is None:
            raise ValueError("Could not estimate the PMs center value")
        cx, cy = vpd_c
        warnings.warn("Could not estimate a better PMs center value")

    return cx, cy


def get_stars_close_center(
    lon: np.ndarray,
    lat: np.ndarray,
    pmRA: np.ndarray,
    pmDE: np.ndarray,
    plx: np.ndarray,
    xy_c: tuple[float, float] | None,
    vpd_c_i: tuple[float, float],
    plx_c: float | None,
    N_cent: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Select the 'N_cent' stars closest to the given center.

    Distances to centers using the vpd_c_i and other available data

    :param lon: Galactic Longitude.
    :type lon: np.ndarray
    :param lat: Galactic Latitude.
    :type lat: np.ndarray
    :param pmRA: Proper motion in Right Ascension.
    :type pmRA: np.ndarray
    :param pmDE: Proper motion in Declination.
    :type pmDE: np.ndarray
    :param plx: Parallax.
    :type plx: np.ndarray
    :param xy_c: Center coordinates in (lon, lat).
    :type xy_c: tuple[float, float] | None
    :param vpd_c_i: Center coordinates in proper motions (pmRA, pmDE).
    :type vpd_c_i: tuple[float, float]
    :param plx_c: Center coordinate in parallax.
    :type plx_c: float | None
    :param N_cent: Number of stars to select.
    :type N_cent: int

    :raises ValueError: If xy_c and plx_c are both None.

    :return: Coordinates and proper motions of the selected stars.
    :rtype: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]
    """
    if xy_c is None and plx_c is None:
        cent = np.array([vpd_c_i])
        data = np.array([pmRA, pmDE]).T
    elif xy_c is None and plx_c is not None:
        cent = np.array([list(vpd_c_i) + [plx_c]])
        data = np.array([pmRA, pmDE, plx]).T
    elif xy_c is not None and plx_c is None:
        cent = np.array([list(xy_c) + list(vpd_c_i)])
        data = np.array([lon, lat, pmRA, pmDE]).T
    elif xy_c is not None and plx_c is not None:
        cent = np.array([list(xy_c) + list(vpd_c_i) + [plx_c]])
        data = np.array([lon, lat, pmRA, pmDE, plx]).T
    else:
        raise ValueError("Either xy_c or plx_c must be given")

    # Closest stars to the selected center
    idx = get_Nd_dists(cent, data)[:N_cent]

    return lon[idx], lat[idx], pmRA[idx], pmDE[idx], plx[idx]


def get_kNN_center(
    N_clust_min: int, data: np.ndarray
) -> tuple[float, float, float, float, float]:
    """Estimate 5D center with kNN.

    :param N_clust_min: Minimum number of stars in the cluster.
    :type N_clust_min: int
    :param data: Array of data.
    :type data: np.ndarray

    :return: Center coordinates in (lon, lat, pmRA, pmDE, plx).
    :rtype: tuple[float, float, float, float, float]
    """
    # Better results are obtained not using the parallax data?
    data_noplx = data[:, :4]  # <-- HARDCODED, TODO

    tree = spatial.KDTree(data_noplx)
    inx = tree.query(data_noplx, k=N_clust_min + 1)
    NN_dist = np.max(inx[0], axis=1)
    # Convert to densities
    dens = 1.0 / NN_dist
    # Sort by largest density
    idxs = np.argsort(-dens)

    # # Use the single star with the largest density
    # cent = np.array([data[idxs[0]]])[0]

    # Use median of 'N_clust_min' stars with largest densities
    x_c, y_c, pmra_c, pmde_c, plx_c = np.median(data[idxs[:N_clust_min]], 0)

    return x_c, y_c, pmra_c, pmde_c, plx_c


def get_2D_center(
    x: np.ndarray, y: np.ndarray, N_max: int = 10_000
) -> tuple[float, float]:
    """Estimate the 2-dimensional center of a cluster, using only its coordinates.

    Find the KDE maximum value pointing to the center coordinates.

    :param x: X coordinates.
    :type x: np.ndarray
    :param y: Y coordinates.
    :type y: np.ndarray
    :param N_max: Maximum number of stars to use for performance reasons, defaults
     to 10_000
    :type N_max: int

    :return: Center coordinates in (x, y).
    :rtype: tuple[float, float]
    """
    values = np.vstack([x, y])
    # Use maximum number for performance
    if values.shape[-1] > N_max:
        # idx = np.random.choice(values.shape[-1], N_max, replace=False)
        xc, yc = np.median([x, y], 1)
        warnings.warn(
            f"\nUsing closest {N_max} stars to center of frame ({xc:.2f}, {yc:.2f})"
            + " to avoid performance issues."
        )
        dist = np.sqrt((x - xc) ** 2 + (y - yc) ** 2)
        idx = np.argsort(dist)[:N_max]
        values = values[:, idx]
        x, y = values

    # Approximate center values
    x_cent_pix, y_cent_pix = get_KDE_cent(values, gd=50)

    # Restrict the KDE to a smaller area to improve performance
    if values.shape[1] > 500:
        rad_x = np.percentile(x, 60) - np.percentile(x, 40)
        rad_y = np.percentile(y, 60) - np.percentile(y, 40)
        # Generate zoom around approx center value to speed things up.
        xmin, xmax = x_cent_pix - rad_x, x_cent_pix + rad_x
        ymin, ymax_z = y_cent_pix - rad_y, y_cent_pix + rad_y
        # Use reduced region around the center.
        msk = (xmin < x) & (x < xmax) & (ymin < y) & (y < ymax_z)
        values = values[:, msk]

    # Final center values
    x_c, y_c = get_KDE_cent(values, gd=100)

    return x_c, y_c


def get_KDE_cent(values: np.ndarray, gd: int) -> tuple[float, float]:
    """Estimate the center coordinates using a Gaussian Kernel Density Estimation.

    :param values: Array of x and y coordinates.
    :type values: np.ndarray
    :param gd: Grid density (number of points).
    :type gd: int

    :return: Center coordinates in (x, y).
    :rtype: tuple[float, float]
    """
    xmin, ymin = values.min(1)
    xmax, ymax = values.max(1)

    # Obtain Gaussian KDE.
    kernel = stats.gaussian_kde(values)

    # Custom bandwith?
    # bdw = kernel.covariance_factor() * np.max(values.std(axis=1)) * .5
    # kernel = stats.gaussian_kde(values, bw_method=bdw / np.max(values.std(axis=1)))

    # Grid density (number of points).
    gd_c = complex(0, gd)
    # Define x,y grid.
    x_grid, y_grid = np.mgrid[xmin:xmax:gd_c, ymin:ymax:gd_c]
    positions = np.vstack([x_grid.ravel(), y_grid.ravel()])
    # Evaluate kernel in grid positions.
    k_pos = kernel(positions)
    # Coordinates of max value in x,y grid (ie: center position).
    x_c, y_c = positions.T[np.argmax(k_pos)]

    return x_c, y_c


def fdens_radius(x, y, cent):
    return np.nan


def king_radius(x, y, cent) -> tuple[float, float, float, float]:
    """ """
    x0, y0, _ = RDPCurve(x, y, cent)

    kp = KingProjectedAnalytic1D()

    def lnlike(theta):
        """Returns a value from np.inf to 0 for an exact match"""
        cd, rc, rt, fd = theta
        return y0 - (kp.evaluate(x0, cd, rc, rt) + fd)

    rc0, rt0, fd0 = 0.25 * max(x0), np.median(x0), max(0, np.median(y0[-10:]))
    cd0 = (max(y0) - fd0) / (1 - 1.0 / np.sqrt(1.0 + (rt0 / rc0) ** 2)) ** 2
    p_init = np.array([cd0, rc0, rt0, fd0])

    min_cd, max_cd = fd0, 10 * max(y0)
    min_rc, max_rc = x0[0], x0[-1]
    min_rt, max_rt = x0[0], 2 * x0[-1]
    min_fd, max_fd = fd0 * 0.1, max(y0)
    bounds = np.array(
        [(min_cd, max_cd), (min_rc, max_rc), (min_rt, max_rt), (min_fd, max_fd)]
    ).T

    # try:
    res = least_squares(lnlike, p_init, bounds=bounds)
    cd0, rc0, rt0, fd0 = res.x
    #     lkl = res.cost * 2
    # except ValueError:
    #     lkl = np.inf

    # import matplotlib.pyplot as plt

    # plt.scatter(x0, y0)
    # plt.show()
    # breakpoint()

    return cd0, rc0, rt0, fd0


def RDPCurve(x, y, cent, RDP_rings=50, rings_rm=0.1, Nmin=10, **kwargs):
    """
    Obtain the RDP using the concentric rings method.

    HARDCODED:
    RDP_rings: number of rings to (initially) try to define
    rings_rm: remove the more conflicting last X% of radii values.
    Nmin: minimum number of stars that a ring should contain. Else, expand it.
    """

    xmax = min(abs(x.max() - cent[0]), abs(x.min() - cent[0]))
    ymax = min(abs(y.max() - cent[1]), abs(y.min() - cent[1]))
    rdp_max = min(xmax, ymax) * 0.95

    # Frame limits
    xy_cent_dist = spatial.distance.cdist([cent], np.array([x, y]).T)[0]

    # Handle the case where int()==0
    # max_i = max(1, int(rings_rm * RDP_rings))
    # The +1 adds a ring accounting for the initial 0. in the array
    radii = np.linspace(0.0, rdp_max, RDP_rings + 1)  # + max_i)[:-max_i]

    # Areas and #stars for all rad values.
    rdp_radii, rdp_points, rdp_stddev = [], [], []
    l_prev, N_in_prev = np.inf, 0.0
    for lw, h in zip(*[radii[:-1], radii[1:]]):
        N_in = ((xy_cent_dist >= lw) & (xy_cent_dist < h)).sum() + N_in_prev

        # If N_in < Nmin take the next ellipse-ring (discard this lw).
        l_now = min(lw, l_prev)

        # Require that at least 'Nmin' stars are within the ellipse-ring.
        if N_in > Nmin:
            ring_area = np.pi * (h**2 - l_now**2)
            # Store RDP parameters.
            rad_med = h if l_now == 0.0 else 0.5 * (l_now + h)
            rdp_radii.append(rad_med)
            rdp_points.append(N_in / ring_area)
            rdp_stddev.append(np.sqrt(N_in) / ring_area)

            # Reset
            l_prev, N_in_prev = np.inf, 0.0

        else:
            l_prev = l_now
            N_in_prev += N_in

    rdp_points = np.array(rdp_points) / np.max(rdp_points)
    return np.array(rdp_radii), rdp_points, rdp_stddev
