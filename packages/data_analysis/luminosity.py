import numpy as np


def main(clp, **kwargs):
    """
    Obtain the Luminosity Function for the field regions and the cluster
    region normalized to their area. Subtract the field curve from the
    cluster curve so as to clean it.

    Uses the main magnitude **after** error rejection.
    """

    # (Main) Magnitudes of all stars AFTER error rejection.
    mmag = np.array(list(zip(*(list(zip(*clp['acpt_stars_c']))[3])))[0])

    # This is the curve for the entire observed frame, normalized to the area
    # of the cluster.
    lf_all, lf_edg_all = np.histogram(
        mmag, bins=50, range=(np.nanmin(mmag), np.nanmax(mmag)))
    x_all = np.concatenate((np.array([0.]), lf_edg_all))
    y_all = np.concatenate(
        (np.array([0.]), lf_all / clp['frame_norm'], np.array([0.])))

    # Obtain histogram for cluster region.
    mag_cl = list(zip(*list(zip(*clp['cl_region_c']))[3]))[0]
    lf_clust, lf_edg_c = np.histogram(
        mag_cl, bins=50, range=(np.nanmin(mag_cl), np.nanmax(mag_cl)))

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
            mag_fl, bins=50, range=(np.nanmin(mag_fl), np.nanmax(mag_fl)))

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

    print("Luminosity function obtained")

    clp['lum_func'] = lum_func
    return clp
