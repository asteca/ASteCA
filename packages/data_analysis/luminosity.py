
import numpy as np
import warnings


def mag_completeness(mmag):
    """
    Calculate the completeness level in each magnitude bin beyond the one
    with the maximum count (ie: the assumed 100% completeness limit)

    This will be used by the isochrone/synthetic cluster fitting algorithm.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mag_hist, bin_edges = np.histogram(
            mmag, 50, range=(np.nanmin(mmag), np.nanmax(mmag)))
    # Index of the bin with the maximum number of stars.
    max_indx = mag_hist.argmax(axis=0)

    # Percentage of stars in each bin beyond the maximum interval (included),
    # assuming the first one is 100%.
    comp_perc = mag_hist[max_indx:] / float(mag_hist[max_indx])

    # Store everything in a single list.
    completeness = [bin_edges, max_indx, comp_perc]

    return completeness


def main(clp, **kwargs):
    """
    Obtain the Luminosity Function for the field regions and the cluster
    region normalized to their area. Subtract the field curve from the
    cluster curve so as to clean it.

    The completeness will be used by the isochrone/synthetic cluster
    fitting algorithm.

    Uses the main magnitude **after** error rejection.
    """

    # (Main) Magnitudes of all stars AFTER error rejection.
    mmag = np.array(list(zip(*(list(zip(*clp['acpt_stars_c']))[3])))[0])

    # Get the completeness level for each magnitude bin.
    completeness = mag_completeness(mmag)

    # This is the curve for the entire observed frame, normalized to the area
    # of the cluster.
    lf_all, lf_edg_all = np.histogram(
        mmag, bins=completeness[1], range=(np.nanmin(mmag), np.nanmax(mmag)))
    x_all = np.concatenate((np.array([0.]), lf_edg_all))
    y_all = np.concatenate(
        (np.array([0.]), lf_all / clp['frame_norm'], np.array([0.])))

    # Obtain histogram for cluster region.
    mag_cl = list(zip(*list(zip(*clp['cl_region_c']))[3]))[0]
    lf_clust, lf_edg_c = np.histogram(
        mag_cl, bins=completeness[1],
        range=(np.nanmin(mag_cl), np.nanmax(mag_cl)))

    # Create arrays adding elements so plt.step will plot the first and last
    # vertical bars.
    x_cl = np.concatenate((np.array([0.]), lf_edg_c))
    y_cl = np.concatenate((np.array([0.]), lf_clust, np.array([0.])))

    # Now for field regions.
    mag_fl = []
    if clp['flag_no_fl_regs_c'] is False:

        # Extract main magnitudes for all stars in all field regions defined.
        for freg in clp['field_regions_c']:
            for star in freg:
                mag_fl.append(star[3][0])

        # Obtain histogram for field region.
        lf_field, lf_edg_f = np.histogram(
            mag_fl, bins=completeness[1],
            range=(np.nanmin(mag_fl), np.nanmax(mag_fl)))

        # Create arrays adding elements so plt.step will plot the first and
        # last vertical bars.
        x_fl = np.concatenate((np.array([0.]), lf_edg_f))
        y_fl = np.concatenate(
            (np.array([0.]), (lf_field / float(len(clp['field_regions_c']))),
             np.array([0.])))
    else:
        print("  WARNING: no field regions defined. Luminosity function\n"
              "  is not cleaned from field star contamination.")
        # Pass dummy lists.
        x_fl, y_fl = [], []

    # Pack values.
    lum_func = [x_cl, y_cl, x_fl, y_fl, x_all, y_all]

    print('LF and completeness magnitude levels obtained.')

    clp['lum_func'], clp['completeness'] = lum_func, completeness
    return clp
