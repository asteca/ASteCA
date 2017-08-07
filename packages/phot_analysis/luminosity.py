
import numpy as np


def mag_completeness(mags):
    '''
    Calculate the completeness level in each magnitude bin beyond the one
    with the maximum count (ie: the assumed 100% completeness limit)
    '''
    # Number of bins.
    bins = int((max(mags) - min(mags)) / 0.1)
    # Get histogram.
    mag_hist, bin_edges = np.histogram(mags, bins)
    # Index of the bin with the maximum number of stars.
    max_indx = mag_hist.argmax(axis=0)
    # Total number of stars beyond the peak bin (included).
    total = sum(mag_hist[max_indx:])
    # Get percentages per interval beyond the maximum interval (included).
    # These values tell me the percentages of stars beyond the magnitude peak
    # that are located inside each magnitude bin. The peak magnitude bin (the
    # first one) will have the biggest percentage.
    comp_perc = [(i * 100.) / total for i in mag_hist[max_indx:]]

    # Store everything in a single list.
    completeness = [bin_edges, max_indx, comp_perc]

    return completeness


def main(clp, mags, **kwargs):
    '''
    Obtain the Luminosity Function for the field regions and the cluster
    region normalized to their area. Subtract the field curve from the
    cluster curve so as to clean it.

    The completeness will be used by the isochrone/synthetic cluster
    fitting algorithm.
    '''
    cl_region, field_regions, flag_no_fl_regs = [
        clp[_] for _ in ['cl_region', 'field_regions', 'flag_no_fl_regs']]

    # Calculate number of bins used by the histograms.
    binwidth = 0.25
    x_min, x_max = min(mags[0]) - 0.5, max(mags[0]) + 0.5
    bins_n = np.arange(int(x_min), int(x_max + binwidth), binwidth)

    # USE MAIN MAGINTUDE.
    mag_cl = zip(*zip(*cl_region)[3])[0]
    # Obtain histogram for cluster region.
    lf_clust, lf_edg_c = np.histogram(mag_cl, bins=bins_n)

    # Create arrays adding elements so plt.step will plot the first and last
    # vertical bars.
    x_cl = np.concatenate((np.array([0.]), lf_edg_c))
    y_cl = np.concatenate((np.array([0.]), lf_clust, np.array([0.])))

    # Now for field regions. USE MAIN MAGINTUDE.
    mag_fl = []
    if flag_no_fl_regs is False:

        for freg in field_regions:
            for star in freg:
                mag_fl.append(star[3][0])

        # Obtain histogram for field region.
        lf_field, lf_edg_f = np.histogram(mag_fl, bins=bins_n)

        # Create arrays adding elements so plt.step will plot the first and
        # last vertical bars.
        x_fl = np.concatenate((np.array([0.]), lf_edg_f))
        y_fl = np.concatenate((np.array([0.]),
                              (lf_field / float(len(field_regions))),
                              np.array([0.])))
    else:
        print("  WARNING: no field regions defined. Luminosity function\n"
              "  is not cleaned from field star contamination.")
        # Pass dummy lists.
        x_fl, y_fl = [], []

    # Pack values.
    lum_func = [x_cl, y_cl, x_fl, y_fl]

    # Get the completeness level for each magnitude bin. This will be used by
    # the isochrone/synthetic cluster fitting algorithm.
    mag_all = list(mag_cl) + mag_fl
    completeness = mag_completeness(mag_all)

    print('LF and completeness magnitude levels obtained.')

    clp['lum_func'], clp['completeness'] = lum_func, completeness
    return clp
