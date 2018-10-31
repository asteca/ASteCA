
import numpy as np
import warnings


def main(clp):
    """
    Calculate the completeness level in each magnitude bin beyond the one
    with the maximum count (ie: the assumed 100% completeness limit)

    This will be used by the isochrone/synthetic cluster fitting algorithm.
    """

    # TODO in place for #301
    complt_flag = False

    # (Main) Magnitudes of accepted stars after error rejection.
    mmag_acpt_c = np.array(list(zip(*(list(zip(*clp['acpt_stars_c']))[3])))[0])

    # Number of stars per bin: 10% of the total number of accepted stars.
    Nbins = min(int(mmag_acpt_c.size * .1), 50)
    if Nbins < 10:
        print("  WARNING: small number of accepted stars used in error\n"
              "  removal function.")
        Nbins = 10

    if len(clp['rjct_stars_c']) > 0:
        # (Main) Magnitudes of error rejected stars.
        mmag_rjct_c = np.array(
            list(zip(*(list(zip(*clp['rjct_stars_c']))[3])))[0])
        all_mags = np.concatenate((mmag_acpt_c, mmag_rjct_c))

        eqN_edges = histedges_equalN(all_mags, Nbins)
        # Histograms for accepted and all stars.
        h_mag_all_c, _ = mmagHist(all_mags, eqN_edges)

        perc_rmvd = 100. * (mmag_rjct_c.size / all_mags.size)

    else:
        eqN_edges = histedges_equalN(mmag_acpt_c, Nbins)
        # Histograms for accepted and all stars.
        h_mag_all_c, _ = mmagHist(mmag_acpt_c, eqN_edges)

        perc_rmvd = 0.

    # Main mag histogram.
    h_mag_acpt_c, _ = mmagHist(mmag_acpt_c, eqN_edges)

    # Estimate error removal function: percentage of stars that remain after
    # error rejection.
    err_rm_perc = h_mag_acpt_c / h_mag_all_c

    if complt_flag:

        # TODO how to combine 'err_per' with the empirical completeness
        # function when they have different bin edges?

        import pdb; pdb.set_trace()  # breakpoint 64a0b3ed //

    else:
        # Get main mag histogram.
        h_mag_acpt_c, bin_edges = mmagHist(mmag_acpt_c, 50)

        # Index of the bin with the maximum number of stars.
        max_indx = h_mag_acpt_c.argmax(axis=0)

        # Percentage of stars in each bin beyond the maximum interval
        # (included), assuming the first one is 100%.
        comp_perc = h_mag_acpt_c[max_indx:] / float(h_mag_acpt_c[max_indx])

        # Store everything in a single list.
        completeness = [bin_edges, max_indx, comp_perc]
        print('Completeness function estimated.')

    clp['completeness'], clp['err_rm_perc'] = completeness,\
        [err_rm_perc, eqN_edges, perc_rmvd]
    return clp


def mmagHist(mmag, edges):
    """
    Histogram of the main magnitude.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mag_hist, bin_edges = np.histogram(
            mmag, edges, range=(np.nanmin(mmag), np.nanmax(mmag)))

    return mag_hist, bin_edges


def histedges_equalN(x, nbin):
    """
    Edges of histogram with equal number of elements.
    Source: https://stackoverflow.com/a/39419049/1391441
    """
    npt = len(x)
    return np.interp(np.linspace(0, npt, nbin + 1),
                     np.arange(npt),
                     np.sort(x))
