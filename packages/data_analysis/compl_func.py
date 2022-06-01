
import numpy as np
from scipy.optimize import curve_fit


def main(pd, clp, cld):
    """
    Estimate the (in)completeness function.
    """
    # Magnitudes
    mags = cld['mags'][0]

    complt_flag = True
    if pd['completeness'][0] == 'n':
        complt_flag = False
        final_compl = [[], [], complt_flag]
    elif pd['completeness'][0] == 'y':
        final_compl = photoAnalysis(mags)
        final_compl += [complt_flag]
        print("(In)completeness function estimated")
    else:
        # Manual completeness function given
        final_compl = manualFunc(pd['completeness'], mags)
        final_compl += [complt_flag]
        print("(In)completeness function read")

    clp['completeness'] = final_compl
    return clp


def photoAnalysis(mags, bins=50):
    """
    If no photometric analysis (in)completeness function is given, the code
    will approximate one using the LF.
    """
    # Magnitude histogram
    mag_histo, mag_edges = np.histogram(mags, bins)

    # Index of the bin with the maximum number of stars.
    max_indx = mag_histo.argmax(axis=0)

    # Percentage of stars in each bin beyond the maximum interval
    # (included), assuming the bin with the maximum value represents 100%.
    comp_perc = np.concatenate((
        np.ones(mag_histo[:max_indx + 1].size),
        mag_histo[max_indx + 1:] / float(mag_histo[max_indx])))

    # Exponential correction
    comp_perc = expCorrection(mag_histo, mag_edges, max_indx, comp_perc)

    # Add a '0.' at the beginning to indicate that all stars with smaller
    # magnitudes than the brightest star observed should be removed by the
    # 'completeness_rm()' process.
    comp_perc = np.array([0.] + list(comp_perc))

    # Percentage of stars that should be *removed* from the synthetic cluster.
    # This is reversed so that the 'completeness_rm()' function is simpler.
    rm_perc = np.clip(1. - comp_perc, a_min=0., a_max=1.)

    return [mag_edges, rm_perc]


def expCorrection(mag_histo, mag_edges, max_indx, comp_perc):
    """
    Exponential correction to the (in)completeness function.

    Fits an exponential (IMF equivalent) to the observed luminosity function
    (ie: the histogram of the main magnitude). The (estimated) number of
    stars lost in the photometry process is obtained for magnitude values
    beyond the bin with the maximum number of stars in the original LF. All
    bins before that one (ie: brighter stars) are considered to be 100%
    complete.
    """
    try:
        # Estimate pre-completeness LF, ie: the LF with no photometric loss
        # of stars.
        def func(x, A, B, C):
            m = np.exp(-x)
            return A * m ** (-B) + C
        x, y = mag_edges[:(max_indx + 1)], mag_histo[:(max_indx + 1)]
        popt, _ = curve_fit(func, x, y)

        # Estimate number of stars beyond the maximum value, using the
        # above fitted LF.
        Nstars = func(mag_edges[(max_indx + 2):], *popt)
        # Correct the completeness loss fraction using this "correct"
        # number of stars.
        comp_perc[(max_indx + 1):] = mag_histo[(max_indx + 1):] / Nstars

    except RuntimeError:
        # LF could not be fit.
        print("  WARNING: no exponential correction applied on "
              "completeness function")

    return comp_perc


def manualFunc(completeness, mags):
    """
    The manual completeness function reports *completeness* values. We
    invert this here because the synthetic cluster function uses the
    *incompleteness* function.
    """

    edges = list(map(float, completeness[::2]))
    # Add 0% completeness towards the bright range
    rm_perc = [0.] + list(map(float, completeness[1::2]))

    # Make sure the rightmost edge contains the dimmest stars, with a 0%
    # completeness
    if mags.max() > edges[-1]:
        edges.append(mags.max() + 0.01)
        rm_perc.append(0.)

    # Invert to incompleteness values
    rm_perc = 1. - np.array(rm_perc)

    return [np.array(edges), rm_perc]
