
import numpy as np
from scipy.optimize import curve_fit
import warnings


def main(clp, cld_i, cld_c):
    """
    Estimate the percentage of stars lost after each process:

    * photometric analysis (completeness)
    * photometric data completeness
    * error removal

    These three processes remove stars from the observed field and need to
    be taken into account in the synthetic cluster generating algorithm.
    """

    # (Main) Magnitudes of accepted stars after error rejection.
    mmag_acpt_c = np.array(list(zip(*(list(zip(*clp['acpt_stars_c']))[3])))[0])

    phot_analy_compl = photoAnalysis(cld_i)

    phot_data_compl = photDataCompl(cld_i, cld_c)

    err_rm_data = errRemv(clp, mmag_acpt_c)

    final_compl = combineCompl(cld_i, phot_analy_compl, mmag_acpt_c)

    clp['phot_analy_compl'], clp['phot_data_compl'], clp['err_rm_data'],\
        clp['completeness'] = phot_analy_compl, phot_data_compl, err_rm_data,\
        final_compl
    return clp


def photoAnalysis(cld_i):
    """
    Stars lost throughout the photometry process.

    If no photometric analysis completeness function is given, the code will
    approximate one. It does this using the main magnitude from the incomplete
    data set (ie: all magnitudes in the input data file), fitting an
    exponential (IMF equivalent) to the observed luminosity function (ie: the
    histogram of the main magnitude). After this, the (estimated) number of
    stars lost in the photometry process is obtained for magnitude values
    beyond the bin with the maximum number of stars in the original LF. All
    bins before that one (ie: brighter stars) are considered to be 100%
    complete.

    """
    # TODO in place for #301
    complt_flag = False

    if complt_flag:
        # comp_b_edges, comp_perc
        pass
        print("Photometric analysis incompleteness function read")

    else:
        # Main magnitudes histogram (full incomplete dataset).
        h_mag_acpt_c, comp_b_edges = mmagHist(cld_i['mags'][0])

        # Index of the bin with the maximum number of stars.
        max_indx = h_mag_acpt_c.argmax(axis=0)

        # Percentage of stars in each bin beyond the maximum interval
        # (included), assuming the previous values are all 100%.
        comp_perc = np.concatenate((
            np.ones(h_mag_acpt_c[:max_indx + 1].size),
            h_mag_acpt_c[max_indx + 1:] / float(h_mag_acpt_c[max_indx])))

        try:
            # Estimate pre-completeness LF, ie: the LF with no photometric loss
            # of stars.
            def func(x, A, B, C):
                m = np.exp(-x)
                return A * m ** (-B) + C
            x, y = comp_b_edges[:(max_indx + 1)], h_mag_acpt_c[:(max_indx + 1)]
            popt, pcov = curve_fit(func, x, y)

            # import matplotlib.pyplot as plt
            # plt.plot(x, y, 'b-', label='data')
            # plt.plot(comp_b_edges, func(comp_b_edges, *popt), 'r-')
            # plt.show()

            # Estimate number of stars beyond the maximum value, using the
            # above fitted LF.
            Nstars = func(comp_b_edges[(max_indx + 2):], *popt)
            # Correct the completeness loss fraction using this "correct"
            # number of stars.
            comp_perc[(max_indx + 1):] = h_mag_acpt_c[(max_indx + 1):] / Nstars

        except RuntimeError:
            # LF could not be fit.
            print("  WARNING: complete LF could not be estimated")

        print("Photometric analysis incompleteness function estimated")

    return [comp_b_edges, comp_perc]


def photDataCompl(cld_i, cld_c):
    """
    """
    # Number of stars per bin: 10% of the total number of accepted stars,
    # between the limits [10, 50].
    Nbins = max(10, min(int(cld_i['mags'][0].size * .1), 50))
    eqN_edges = histedges_equalN(cld_i['mags'][0], Nbins)

    # Histograms for accepted and all stars.
    h_mag_i, _ = mmagHist(cld_i['mags'][0], eqN_edges)
    # Accepted main magnitudes histogram (complete dataset).
    h_mag_c, _ = mmagHist(cld_c['mags'][0], eqN_edges)

    # Percentage of stars that remain after photom completeness.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        err_rm_perc = h_mag_c / h_mag_i.astype(float)

    # TODO I think Python3 does not need this 'float'
    perc_rmvd = 100. *\
        (1. - cld_c['mags'][0].size / float(cld_i['mags'][0].size))

    print("Photometric completeness data function estimated")
    return [eqN_edges, err_rm_perc, perc_rmvd]


def errRemv(clp, mmag_acpt_c):
    """
    """
    eqN_edges, err_rm_perc, perc_rmvd =\
        [mmag_acpt_c.min(), mmag_acpt_c.max()], np.ones(1), 0.
    if len(clp['rjct_stars_c']) > 0:
        # (Main) Magnitudes of error rejected stars.
        mmag_rjct_c = np.array(
            list(zip(*(list(zip(*clp['rjct_stars_c']))[3])))[0])
        all_mags = np.concatenate((mmag_acpt_c, mmag_rjct_c))

        # Number of stars per bin: 10% of the total number of accepted stars,
        # between the limits [10, 50].
        Nbins = max(10, min(int(mmag_acpt_c.size * .1), 50))
        eqN_edges = histedges_equalN(all_mags, Nbins)

        # Histograms for accepted and all stars.
        h_mag_all_c, _ = mmagHist(all_mags, eqN_edges)

        # Accepted main magnitudes histogram (complete dataset).
        h_mag_acpt_c, _ = mmagHist(mmag_acpt_c, eqN_edges)

        # Estimate error removal function: percentage of stars that remain
        # after error rejection.
        err_rm_perc = h_mag_acpt_c / h_mag_all_c.astype(float)

        # TODO I think Python3 does not need this 'float'
        perc_rmvd = 100. * (mmag_rjct_c.size / float(all_mags.size))

    print("Error removal function estimated")
    return [eqN_edges, err_rm_perc, perc_rmvd]


def combineCompl(cld_i, phot_analy_compl, mmag_acpt_c):
    """
    """
    comp_b_edges, comp_perc = phot_analy_compl

    # Generate edges of N bins with approx equal number of elements each.
    eqN_edges = histedges_equalN(cld_i['mags'][0], 20)

    # Correct the initial LF using the completeness loss function.
    h_mag_i, _ = mmagHist(cld_i['mags'][0], comp_b_edges)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        h_mag_i_full = h_mag_i / comp_perc
    # Rebin corrected full LF according to 'eqN_edges' array.
    h_mag_i_full = rebin(comp_b_edges, h_mag_i_full, eqN_edges)

    # Accepted main magnitudes histogram (complete dataset).
    h_mag_acpt_c, _ = mmagHist(mmag_acpt_c, eqN_edges)

    # Percentage of stars that should be *removed* from the synthetic cluster.
    # This is reversed with respect to the other percentages so that
    # the completeness_rm() function is simpler.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        err_rm_perc = 1. - (h_mag_acpt_c / h_mag_i_full)

    # TODO I think Python3 does not need this 'float'
    perc_rmvd = 100. * (1. - mmag_acpt_c.size / float(cld_i['mags'][0].size))

    print("Combined completeness function estimated")
    return [eqN_edges, err_rm_perc, perc_rmvd]


def mmagHist(mmag, edges=50):
    """
    Magnitude histogram.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mag_hist, b_edges = np.histogram(
            mmag, edges, range=(np.nanmin(mmag), np.nanmax(mmag)))

    return mag_hist, b_edges


def histedges_equalN(x, nbin):
    """
    Edges of histogram with equal number of elements.
    Source: https://stackoverflow.com/a/39419049/1391441
    """
    npt = len(x)
    return np.interp(np.linspace(0, npt, nbin + 1),
                     np.arange(npt),
                     np.sort(x))


def rebin(x1, y1, x2):
    """
    Copyright (c) 2011, Joshua M. Hykes

    Source: https://github.com/jhykes/rebin

    ---------------------------------------

    Rebin histogram values y1 from old bin edges x1 to new edges x2.

    Input
    -----
     * x1 : m+1 array of old bin edges.
     * y1 : m array of old histogram values. This is the total number in
              each bin, not an average.
     * x2 : n+1 array of new bin edges.

    Returns
    -------
     * y2 : n array of rebinned histogram values.

    The rebinning algorithm assumes that the counts in each old bin are
    uniformly distributed in that bin.

    Bins in x2 that are entirely outside the range of x1 are assigned 0.
    """

    # Divide y1 by bin widths.
    #  This converts y-values from bin total to bin average over bin width.
    x1_bin_widths = np.ediff1d(x1)
    y1_ave = y1 / x1_bin_widths

    # allocating y2 vector
    n = x2.size - 1
    y2 = np.zeros(n, dtype=y1.dtype)

    i_place = np.searchsorted(x1, x2)

    # find out where x2 intersects with x1, this will determine which x2 bins
    # we need to consider
    start_pos = 0
    end_pos = n

    start_pos_test = np.where(i_place == 0)[0]
    if start_pos_test.size > 0:
        start_pos = start_pos_test[-1]

    end_pos_test = np.where(i_place == x1.size)[0]
    if end_pos_test.size > 0:
        end_pos = end_pos_test[0]

    # the first bin totally covers x1 range
    if (start_pos == end_pos - 1 and i_place[start_pos] == 0 and
            i_place[start_pos + 1] == x1.size):
        sub_edges = x1
        sub_dx = np.ediff1d(sub_edges)
        sub_y_ave = y1_ave

        y2[start_pos] = np.sum(sub_dx * sub_y_ave)

        return y2

    # the first bin overlaps lower x1 boundary
    if i_place[start_pos] == 0 and start_pos < end_pos:
        x2_lo, x2_hi = x2[start_pos], x2[start_pos + 1]
        i_lo, i_hi = i_place[start_pos], i_place[start_pos + 1]

        sub_edges = np.hstack([x1[i_lo:i_hi], x2_hi])
        sub_dx = np.ediff1d(sub_edges)
        sub_y_ave = y1_ave[i_lo: i_hi]

        y2[start_pos] = np.sum(sub_dx * sub_y_ave)

        start_pos += 1

    # the last bin overlaps upper x1 boundary
    if (i_place[end_pos] == x1.size and start_pos < end_pos):
        x2_lo, x2_hi = x2[end_pos - 1], x2[end_pos]
        i_lo, i_hi = i_place[end_pos - 1], i_place[end_pos]

        sub_edges = np.hstack([x2_lo, x1[i_lo:i_hi]])
        sub_dx = np.ediff1d(sub_edges)
        sub_y_ave = y1_ave[i_lo - 1:i_hi]

        y2[end_pos - 1] = np.sum(sub_dx * sub_y_ave)

        end_pos -= 1

    if start_pos < end_pos:
        # deal with interior bins

        # deal with whole parts of bin that are spanned
        cum_sum = np.cumsum(y1)
        running_sum = (
            cum_sum[i_place[start_pos + 1:end_pos + 1] - 2] -
            cum_sum[i_place[start_pos:end_pos] - 1])

        y2[start_pos:end_pos] += running_sum

        # deal with fractional start of bin
        p_sub_dx = x1[i_place[start_pos:end_pos]] - x2[start_pos:end_pos]
        sub_y_ave = y1_ave[i_place[start_pos:end_pos] - 1]

        y2[start_pos:end_pos] += p_sub_dx * sub_y_ave

        # deal with fractional end of bin
        p_sub_dx = x2[start_pos + 1:end_pos + 1] -\
            x1[i_place[start_pos + 1:end_pos + 1] - 1]
        sub_y_ave = y1_ave[i_place[start_pos + 1:end_pos + 1] - 1]

        y2[start_pos:end_pos] += p_sub_dx * sub_y_ave

    return y2
