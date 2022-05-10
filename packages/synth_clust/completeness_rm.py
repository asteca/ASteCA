
import numpy as np


def main(isoch, completeness, rand_unif_vals):
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

    # 'comp_perc' means the percentage of stars that should
    # be *REMOVED* from each mag bin.
    bin_edges, comp_perc, compl_flag = completeness

    # If no completeness function was defined or no star in the isochrone
    # exists beyond the completeness magnitude level, skip the process.
    if compl_flag is False or np.max(isoch[0]) < bin_edges[0]:
        return np.array(isoch), np.array([])

    # Indexes of stars in 'isoch[0]' whose main magnitude
    # value falls between the ranges given by 'bin_edges'.
    #
    # Magnitude values *below* the minimum magnitude edge will be
    # assigned the integer '0'. For each star these are the indexes that
    # point to their corresponding position inside 'bin_edges'
    c_indx = np.searchsorted(bin_edges, isoch[0], side='left')

    # This mask points to the stars that should be KEPT
    msk = comp_perc[c_indx] < rand_unif_vals[:c_indx.size]
    isoch_compl = np.array(isoch[:, msk])

    return isoch_compl, msk

    # #
    # # DEPRECATED OLD METHOD 04/22
    # msk = None

    # # Equivalent to np.histogram(isoch_binar[0], bin_edges)[0]
    # count = np.bincount(c_indx, minlength=len(comp_perc))

    # # Round to integer and clip at '0' so there are no negative values.
    # di = np.rint(count * comp_perc).astype(int).clip(0)

    # # The stars are already shuffled in 'mass_interp', so this selection
    # # of the first 'd' elements is not removing a given type of star over
    # # any other.
    # d_i = []
    # for i, d in enumerate(di):
    #     d_i.append(np.where(c_indx == i)[0][:d])
    # d_i = np.concatenate(d_i)

    # # Remove stars pointed to by 'd_i' from *all* the sub-arrays in
    # # 'isoch_binar'.
    # isoch_compl = np.delete(isoch_binar, d_i, axis=1)

