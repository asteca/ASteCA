import warnings

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from scipy import spatial, stats


def radec2lonlat(ra: np.ndarray, dec: np.ndarray) -> np.ndarray:
    """Convert from RA, DEC to Galactic Longitude, Latitude.

    :param ra: Right ascension.
    :type ra: np.ndarray
    :param dec: Declination.
    :type dec: np.ndarray
    :return: Galactic longitude and latitude.
    :rtype: np.ndarray
    """


def lonlat2radec(lon: np.ndarray, lat: np.ndarray) -> list:
    """Convert from Galactic Longitude, Latitude to RA, DEC.

    :param lon: Galactic longitude.
    :type lon: np.ndarray
    :param lat: Galactic latitude.
    :type lat: np.ndarray
    :return: Right ascension and declination.
    :rtype: list
    """


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


def get_Nd_dists(
    cents: np.ndarray,  np.ndarray, dists_flag: bool = False
) -> np.ndarray:
    """Obtain indexes and distances of stars to the given center

    :param cents: Center coordinates.
    :type cents: np.ndarray
    :param  Array of data.
    :type  np.ndarray
    :param dists_flag: If True, return distances instead of indexes, defaults to False
    :type dists_flag: bool, optional
    :return: Indexes or distances of stars to the given center.
    :rtype: np.ndarray
    """


def get_5D_center(
    lon: np.ndarray,
    lat: np.ndarray,
    pmRA: np.ndarray,
    pmDE: np.ndarray,
    plx: np.ndarray,
    xy_c: list[float] | None,
    vpd_c: list[float] | None,
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
    :param pmRA: Proper motion in RA.
    :type pmRA: np.ndarray
    :param pmDE: Proper motion in DEC.
    :type pmDE: np.ndarray
    :param plx: Parallax.
    :type plx: np.ndarray
    :param xy_c: Center in (lon, lat), defaults to None
    :type xy_c: list[float] | None, optional
    :param vpd_c: Center in (pmRA, pmDE), defaults to None
    :type vpd_c: list[float] | None, optional
    :param plx_c: Center in parallax, defaults to None
    :type plx_c: float | None, optional
    :param N_clust_min: Minimum number of stars in the cluster.
    :type N_clust_min: int
    :param N_clust_max: Maximum number of stars in the cluster.
    :type N_clust_max: int
    :return: Center in (lon, lat, pmRA, pmDE, plx).
    :rtype: tuple[float, float, float, float, float]
    """


def filter_pms_stars(
    xy_c: list[float] | None,
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

    :param xy_c: Center in (lon, lat), defaults to None
    :type xy_c: list[float] | None, optional
    :param plx_c: Center in parallax, defaults to None
    :type plx_c: float | None, optional
    :param lon: Galactic Longitude.
    :type lon: np.ndarray
    :param lat: Galactic Latitude.
    :type lat: np.ndarray
    :param pmRA: Proper motion in RA.
    :type pmRA: np.ndarray
    :param pmDE: Proper motion in DEC.
    :type pmDE: np.ndarray
    :param plx: Parallax.
    :type plx: np.ndarray
    :param N_cent: Number of stars to select.
    :type N_cent: int
    :return: Proper motions of the selected stars.
    :rtype: tuple[np.ndarray, np.ndarray]
    """


def get_pms_center(
    vpd_c: list[float] | None,
    N_clust_min: int,
    pmRA: np.ndarray,
    pmDE: np.ndarray,
    N_bins: int = 50,
    zoom_f: int = 4,
    N_zoom: int = 10,
) -> list[float]:
    """Estimate the center in proper motion space.

    :param vpd_c: Center in (pmRA, pmDE), defaults to None
    :type vpd_c: list[float] | None, optional
    :param N_clust_min: Minimum number of stars in the cluster.
    :type N_clust_min: int
    :param pmRA: Proper motion in RA.
    :type pmRA: np.ndarray
    :param pmDE: Proper motion in DEC.
    :type pmDE: np.ndarray
    :param N_bins: Number of bins for the histogram, defaults to 50
    :type N_bins: int, optional
    :param zoom_f: Zoom factor, defaults to 4
    :type zoom_f: int, optional
    :param N_zoom: Number of zoom iterations, defaults to 10
    :type N_zoom: int, optional
    :return: Center in (pmRA, pmDE).
    :rtype: list[float]
    """


def get_stars_close_center(
    lon: np.ndarray,
    lat: np.ndarray,
    pmRA: np.ndarray,
    pmDE: np.ndarray,
    plx: np.ndarray,
    xy_c: list[float] | None,
    vpd_c: list[float] | None,
    plx_c: float | None,
    N_cent: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Distances to centers using the vpd_c and other available data

    :param lon: Galactic Longitude.
    :type lon: np.ndarray
    :param lat: Galactic Latitude.
    :type lat: np.ndarray
    :param pmRA: Proper motion in RA.
    :type pmRA: np.ndarray
    :param pmDE: Proper motion in DEC.
    :type pmDE: np.ndarray
    :param plx: Parallax.
    :type plx: np.ndarray
    :param xy_c: Center in (lon, lat), defaults to None
    :type xy_c: list[float] | None, optional
    :param vpd_c: Center in (pmRA, pmDE), defaults to None
    :type vpd_c: list[float] | None, optional
    :param plx_c: Center in parallax, defaults to None
    :type plx_c: float | None, optional
    :param N_cent: Number of stars to select.
    :type N_cent: int
    :return: Selected stars' (lon, lat, pmRA, pmDE, plx).
    :rtype: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]
    """


def get_kNN_center(N_clust_min: int,  np.ndarray) -> np.ndarray:
    """Estimate 5D center with kNN.

    :param N_clust_min: Minimum number of stars in the cluster.
    :type N_clust_min: int
    :param  Array of data.
    :type  np.ndarray
    :return: Center coordinates.
    :rtype: np.ndarray
    """


def get_2D_center(x: np.ndarray, y: np.ndarray, N_max: int = 10000) -> tuple[float, float]:
    """Estimate the 2-dimensional center of a cluster, using only its coordinates.

    Find the KDE maximum value pointing to the center coordinates.

    :param x: X coordinates.
    :type x: np.ndarray
    :param y: Y coordinates.
    :type y: np.ndarray
    :param N_max: Maximum number of stars to use, defaults to 10000
    :type N_max: int, optional
    :return: Center coordinates.
    :rtype: tuple[float, float]
    """


def get_XY(values: np.ndarray, gd: int) -> tuple[float, float]:
    """Estimate the center coordinates using a Gaussian KDE.

    :param values: Array of coordinates.
    :type values: np.ndarray
    :param gd: Grid density (number of points).
    :type gd: int
    :return: Center coordinates.
    :rtype: tuple[float, float]
    """
