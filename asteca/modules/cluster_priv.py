import warnings
import numpy as np
from scipy import spatial
from scipy import stats
import astropy.units as u
from astropy.coordinates import SkyCoord


def radec2lonlat(ra, dec):
    gc = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)
    lb = gc.transform_to("galactic")
    lon, lat = lb.l.value, lb.b.value
    return [lon, lat]


def lonlat2radec(lon, lat):
    gc = SkyCoord(l=lon * u.degree, b=lat * u.degree, frame="galactic")
    ra, dec = gc.fk5.ra.value, gc.fk5.dec.value
    return [ra, dec]


def reject_nans(data):
    """Remove nans in 'data'"""
    msk_all = []
    # Process each dimension separately
    for arr in data:
        # Identify non-nan data
        msk = ~np.isnan(arr)
        # Keep this non-nan data
        msk_all.append(msk.data)
    # Combine into a single mask
    msk_accpt = np.logical_and.reduce(msk_all)

    # Indexes that survived
    idx_clean = np.arange(data.shape[1])[msk_accpt]

    return idx_clean, data.T[msk_accpt].T


def get_Nd_dists(cents, data, dists_flag=False):
    """Obtain indexes and distances of stars to the given center"""
    # Distances to center
    dist_Nd = spatial.distance.cdist(data, cents).T[0]
    if dists_flag:
        # Return the distances
        return dist_Nd

    # Indexes that sort the distances
    d_idxs = dist_Nd.argsort()
    # Return the indexes that sort the distances
    return d_idxs


def get_5D_center(
    lon, lat, pmRA, pmDE, plx, xy_c, vpd_c, plx_c, N_clust_min, N_clust_max
):
    """
    Estimate the 5-dimensional center of a cluster:

    1. Keep only 'N_cent' stars if xy_c or plx_c are given
    2. (Re)Estimate the center in PMs (the value can be given as input)
    3. Obtain the 'N_cent' stars closest to the available center values
    4. Estimate the 5-dimensional final center using kNN

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
    vpd_c = get_pms_center(vpd_c, N_clust_min, pmRA_i, pmDE_i)

    # Get N_cent stars closest to vpd_c and given xy and/or plx centers
    lon_i, lat_i, pmRA_i, pmDE_i, plx_i = get_stars_close_center(
        lon, lat, pmRA, pmDE, plx, xy_c, vpd_c, plx_c, N_cent
    )

    # kNN center
    x_c, y_c, pmra_c, pmde_c, plx_c = get_kNN_center(
        N_clust_min, np.array([lon_i, lat_i, pmRA_i, pmDE_i, plx_i]).T
    )

    return x_c, y_c, pmra_c, pmde_c, plx_c


def filter_pms_stars(xy_c, plx_c, lon, lat, pmRA, pmDE, plx, N_cent):
    """If either xy_c or plx_c values are given, select the 'N_cent' stars
    closest to this 1D/2D/3D center, and return their proper motions.
    """
    # Create arrays with required shape
    if xy_c is None and plx_c is not None:
        cent = np.array([[plx_c]])
        data = np.array([plx]).T
    elif xy_c is not None and plx_c is None:
        cent = np.array([xy_c])
        data = np.array([lon, lat]).T
    else:
        # xy_c is not None and plx_c is not None
        cent = np.array([list(xy_c) + [plx_c]])
        data = np.array([lon, lat, plx]).T

    # Closest stars to the selected center
    idx = get_Nd_dists(cent, data)[:N_cent]
    pmRA_i, pmDE_i = pmRA[idx], pmDE[idx]

    return pmRA_i, pmDE_i


def get_pms_center(vpd_c, N_clust_min, pmRA, pmDE, N_bins=50, zoom_f=4, N_zoom=10):
    """ """
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

    if vpd_c is None and cxym is None:
        raise Exception("Could not estimate the PMs center value")

    if cxym is not None:
        cx, cy = cxym
    else:
        cx, cy = vpd_c
        warnings.warn("Could not estimate a better PMs center value")

    return [cx, cy]


def get_stars_close_center(lon, lat, pmRA, pmDE, plx, xy_c, vpd_c, plx_c, N_cent):
    """
    Distances to centers using the vpd_c and other available data
    """
    if xy_c is None and plx_c is None:
        cent = np.array([vpd_c])
        data = np.array([pmRA, pmDE]).T
    elif xy_c is None and plx_c is not None:
        cent = np.array([vpd_c + [plx_c]])
        data = np.array([pmRA, pmDE, plx]).T
    elif xy_c is not None and plx_c is None:
        cent = np.array([list(xy_c) + vpd_c])
        data = np.array([lon, lat, pmRA, pmDE]).T
    else:
        # xy_c is not None and plx_c is not None
        cent = np.array([list(xy_c) + vpd_c + [plx_c]])
        data = np.array([lon, lat, pmRA, pmDE, plx]).T

    # Closest stars to the selected center
    idx = get_Nd_dists(cent, data)[:N_cent]

    return lon[idx], lat[idx], pmRA[idx], pmDE[idx], plx[idx]


def get_kNN_center(N_clust_min, data):
    """Estimate 5D center with kNN."""

    # Better results are obtained not using the parallax data?
    data_noplx = data[:, :4]  # <-- HARDCODED

    tree = spatial.cKDTree(data_noplx)
    inx = tree.query(data_noplx, k=N_clust_min + 1)
    NN_dist = inx[0].max(1)
    # Convert to densities
    dens = 1.0 / NN_dist
    # Sort by largest density
    idxs = np.argsort(-dens)

    # # Use the single star with the largest density
    # cent = np.array([data[idxs[0]]])[0]

    # Use median of 'N_clust_min' stars with largest densities
    cent = np.median(data[idxs[:N_clust_min]], 0)

    return cent


def get_2D_center(x, y, N_max=10000):
    """Estimate the 2-dimensional center of a cluster, using only its coordinates.

    Find the KDE maximum value pointing to the center coordinates.
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
    x_cent_pix, y_cent_pix = get_XY(values, gd=50)

    # Restrict the KDE to a smaller area to improve performance
    rad_x = np.percentile(x, 60) - np.percentile(x, 40)
    rad_y = np.percentile(y, 60) - np.percentile(y, 40)
    # Generate zoom around approx center value to speed things up.
    xmin, xmax = x_cent_pix - rad_x, x_cent_pix + rad_x
    ymin, ymax_z = y_cent_pix - rad_y, y_cent_pix + rad_y
    # Use reduced region around the center.
    msk = (xmin < x) & (x < xmax) & (ymin < y) & (y < ymax_z)

    # Final center values
    x_c, y_c = get_XY(values[:, msk], gd=100)

    return x_c, y_c


def get_XY(values, gd):
    """ """
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
