
import numpy as np


def main(isoch_binar, completeness):
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

    if completeness[-1] is False:
        return isoch_binar

    # 'comp_perc' means the percentage of stars that should
    # be *REMOVED* from each mag bin.
    bin_edges, comp_perc = completeness[:2]

    # If stars exist in the isochrone beyond the completeness magnitude
    # level, then apply the removal of stars. Otherwise, skip it.
    if np.max(isoch_binar[0]) > bin_edges[0]:

        # Indexes of stars in 'isoch_binar[0]' whose main magnitude
        # value falls between the ranges given by 'bin_edges'.
        #
        # Magnitude values *below* the minimum magnitude edge will be
        # assigned the integer '0'.
        c_indx = np.searchsorted(bin_edges, isoch_binar[0], side='left')

        # Equivalent to np.histogram(isoch_binar[0], bin_edges)[0]
        count = np.bincount(c_indx, minlength=len(comp_perc))

        # Round to integer and clip at '0' so there are no negative values.
        di = np.rint(count * comp_perc).astype(int).clip(0)

        # The stars are already shuffled in 'mass_interp', so this selection
        # of the first 'd' elements is not removing a given type of star over
        # any other.
        d_i = []
        for i, d in enumerate(di):
            d_i.append(np.where(c_indx == i)[0][:d])
        d_i = np.concatenate(d_i)

        # Remove stars pointed to by 'd_i' from *all* the sub-arrays in
        # 'isoch_binar'.
        isoch_compl = np.delete(isoch_binar, d_i, axis=1)
        #
        # import matplotlib.pyplot as plt
        # plt.hist(isoch_binar[0], bins=bin_edges, histtype='step', label="orig")
        # plt.hist(isoch_compl[0], bins=bin_edges, histtype='step', label="new")
        # plt.hist(isoch_compl2[0], bins=bin_edges, histtype='step', ls=':',
        #          label="old")
        # plt.legend()
        # plt.show()

    else:
        isoch_compl = isoch_binar

    return isoch_compl
