
import numpy as np


def main(isoch_binar, completeness):
    """
    Remove a number of stars according to the percentages of star loss found in
    the 'mag_completeness' function of the luminosity module, for the real
    observation.
    """

    # Remember that 'comp_perc' here means the percentage of stars that should
    # be *REMOVED* from each mag range|bin.
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

        # # DEPRECATED 03/12/19 #445
        # # The minimum length is that of the 'comp_perc' list plus one,
        # # so after removing the '0' elements both lists will have the same
        # # length.
        # count = np.bincount(c_indx, minlength=len(comp_perc) + 1)[1:]
        # di = np.rint(count * comp_perc).astype(int).clip(0)
        # # Actual indexes of stars, stored in each edge range.
        # rang_indx = idxFind(len(bin_edges), c_indx)
        # # Pick a number (given by the list 'di') of random elements in
        # # each range. Those are the indexes of the elements that
        # # should be *removed* from the sub-lists.
        # d_i = indxRem(di, rang_indx, cmpl_rnd)

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


# DEPRECATED 03/12/ #445
# def indxRem(di, rang_indx, cmpl_rnd):
#     """
#     Select a fixed number (given by 'di') of random indexes in 'rang_indx'.
#     These correspond to the stars that will be removed in each magnitude
#     range.

#     Source: https://stackoverflow.com/a/46079837/1391441
#     """
#     lens = np.array([len(_) for _ in rang_indx])
#     di0 = np.minimum(lens, di)
#     invalid_mask = lens[:, None] <= np.arange(lens.max())

#     # Create a 2D random array in interval [0,1) to cover the max length of
#     # subarrays.
#     rand_nums = np.copy(cmpl_rnd[:len(lens) * lens.max()].reshape(
#         len(lens), lens.max()))

#     # For each subarray, set the invalid places to 1.0. Get argsort for each
#     # row. Those 1s corresponding to the invalid places would stay at the back
#     # because there were no 1s in the original random array. Thus, we have the
#     # indices array.
#     rand_nums[invalid_mask] = 1
#     # Slice each row of those indices array to the extent of the lengths
#     # listed in di.
#     shuffled_indx = np.argpartition(rand_nums, lens - 1, axis=1)

#     # Start a loop and slice each subarray from 'rang_indx' using those sliced
#     # indices.
#     out = []
#     for i, all_idx in enumerate(shuffled_indx):
#         if lens[i] > 0:
#             slice_idx = all_idx[:di0[i]]
#             out += rang_indx[i][slice_idx].tolist()

#     return np.asarray(out)


# def idxFind(N, c_indx):
#     """
#     Store the actual indexes of stars in the accepted edge ranges, stored in
#     each corresponding range.
#     """
#     # Reject stars in the 0th position. These are stars below the value
#     # where the completeness loss starts.
#     mask = (c_indx > 0)
#     # Keep those stars with indexes in the accepted magnitude range.
#     c_mask = c_indx[mask]
#     # Ordered indexes for the masked stars.
#     indices = np.arange(c_indx.size)[mask]
#     # Indexes that would sort 'c_mask'.
#     sorting_idx = np.argsort(c_mask, kind='mergesort')
#     # Keep only the ordered indexes that are associated with 'c_mask'.
#     ind_sorted = indices[sorting_idx]
#     # Indexes of ordered indexes (N) positioned into 'c_mask'.
#     x = np.searchsorted(
#         c_mask, range(N), side='right', sorter=sorting_idx)

#     # Store star indices into each edge range.
#     rang_indx = [ind_sorted[x[i]:x[i + 1]] for i in range(N - 1)]

#     return rang_indx
