
import numpy as np
from scipy.spatial.distance import cdist


def main(clp, x, y, mags, **kwargs):
    """
    Calculate the contamination index value. This parameter is defined as the
    ratio of field stars density over the density of stars in the cluster
    region. Uses the 'incomplete' data.

    CI ~ 0.0 -> the field contamination in the cluster region is very small.
    CI ~ 0.5 -> an equal number of field stars and cluster members are
    expected inside the cluster region.
    CI ~ 1.0 -> there are no expected cluster members inside the cluster region
    (not a good sign).
    """

    # If the cluster's area is less than 10% of the full frame's area, don't
    # estimate the CI.
    frame_area = np.ptp(x) * np.ptp(y)
    if (frame_area - clp['cl_area']) / frame_area > .1:

        # Count the total number of stars within the defined cluster region
        # (including stars with rejected photometric errors)
        dist = cdist([clp['kde_cent']], np.array([x, y]).T)[0]
        n_in_cl_reg = (dist < clp['clust_rad']).sum()

        # Estimated number of members versus maximum magnitude. For plotting.
        membvsmag = NmembVsMag(
            x, y, mags[0], clp['clust_rad'], clp['cl_area'], clp['field_dens'],
            dist)

        # Final contamination index.
        cont_index, n_memb_i, _ = CIfunc(
            n_in_cl_reg, clp['field_dens'], clp['cl_area'])
        # Used in the King profile fitting
        n_memb_i = int(round(n_memb_i))
        if n_memb_i < 2:
            print("  WARNING: The estimated number of members is < 2")
            n_memb_i = 1

        if cont_index >= 1.:
            print("  WARNING: contamination index value is very large: "
                  "{:.2f}".format(cont_index))
        else:
            print("Contamination index obtained ({:.2f})".format(cont_index))
    else:
        print("  WARNING: cluster radius is too large to obtain\n"
              "  a reliable contamination index value")
        cont_index, n_memb_i, membvsmag = np.nan, np.nan, []

    clp['cont_index'], clp['n_memb_i'], clp['membvsmag'] =\
        cont_index, n_memb_i, membvsmag
    return clp


def CIfunc(n_in_cl_reg, field_dens, area):
    """
    """
    # Estimated number of field stars in the area.
    n_fl = field_dens * area
    # Estimated number of members in the area.
    n_memb = n_in_cl_reg - n_fl

    # Star density in the cluster region.
    cl_dens = n_in_cl_reg / area
    # Contamination index.
    CI = field_dens / cl_dens

    return CI, n_memb, n_fl


def NmembVsMag(x, y, mag, clust_rad, cl_area, field_dens, cent_dists):
    """
    Number of members versus magnitude cut. Used for plotting.
    """

    area_tot = (np.nanmax(x) - np.nanmin(x)) * (np.nanmax(y) - np.nanmin(y))
    area_out = area_tot - cl_area

    membvsmag = []
    mag_ranges = np.linspace(np.nanmin(mag), np.nanmax(mag), 11)
    # handle possible nan values
    mag[np.isnan(mag)] = np.inf
    for i, mmax in enumerate(mag_ranges[1:]):
        msk = mag < mmax
        if msk.sum() > 2:
            n_in_cl_reg = (cent_dists[msk] < clust_rad).sum()
            # Use the global field density value when the maximum magnitude
            # is used.
            if i == 9:
                fdens = field_dens
            else:
                Ntot = msk.sum()
                fdens = (Ntot - n_in_cl_reg) / area_out
            n_fl = fdens * cl_area
            n_memb_i = max(0, int(round(n_in_cl_reg - n_fl)))
            membvsmag.append([mmax, n_memb_i])

    return np.array(membvsmag).T
