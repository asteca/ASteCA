import numpy as np


def main(isoch_extin, completeness, rand_unif_vals):
    """
    Remove a number of stars according to the percentages of star loss found in
    the 'mag_completeness' function of the luminosity module, for the real
    observation.

    The 'bin_edges, comp_perc' arrays must be understood as follows:

    bin_edges[0] --> % of elements rmvd for mag <= comp_perc[0]
    bin_edges[1] --> % of elements rmvd for comp_perc[0] < mag <= comp_perc[1]
    etc...

    The first element of both arrays is important because it applies a
    bight end (minimum magnitude) cut.
    """
    if completeness is None:
        return np.array(isoch_extin), 0

    # 'comp_perc' means the percentage of stars that should be *REMOVED* from each
    # mag bin.
    bin_edges, comp_perc = completeness

    # If no star in the isochrone exists beyond the completeness magnitude level, skip
    # the process.
    if np.max(isoch_extin[0]) < bin_edges[0]:
        return np.array(isoch_extin), 0

    # Indexes of stars in 'isoch_extin[0]' whose main magnitude
    # value falls between the ranges given by 'bin_edges'.
    #
    # Magnitude values *below* the minimum magnitude edge will be
    # assigned the integer '0'. For each star these are the indexes that
    # point to their corresponding position inside 'bin_edges'
    c_indx = np.searchsorted(bin_edges, isoch_extin[0], side="left")

    # This mask points to the stars that should be KEPT
    msk = comp_perc[c_indx] < rand_unif_vals[: c_indx.size]
    isoch_compl = np.array(isoch_extin[:, msk])

    return isoch_compl, (~msk).sum()
