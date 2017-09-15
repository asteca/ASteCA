
import numpy as np
import itertools


def main_orig(isoch_binar, binar_idx0, completeness):
    '''
    Remove a number of stars according to the percentages of star loss found in
    the mag_completeness function of the luminosity module, for the real
    observation.
    '''
    # If stars exist in the isochrone beyond the completeness magnitude
    # level, then apply the removal of stars. Otherwise, skip it.
    # completeness = [bin_edges, max_indx, comp_perc]
    if max(isoch_binar[0]) > completeness[0][completeness[1]]:

        # Get histogram. completeness[0] = bin_edges of the observed
        # region histogram.
        synth_mag_hist = np.histogram(isoch_binar[0], completeness[0])[0]
        pi = completeness[2]
        n1, p1 = synth_mag_hist[completeness[1]], pi[0]
        di = np.around((synth_mag_hist[completeness[1]:] -
                        (n1 / p1) * np.asarray(pi)), 0)

        # Store indexes of *all* elements in isoch_binar whose main magnitude
        # value falls between the ranges given.
        c_indx = np.searchsorted(completeness[0][completeness[1]:],
                                 isoch_binar[0], side='left')
        N = len(completeness[0][completeness[1]:])
        mask = (c_indx > 0) & (c_indx < N)
        elements = c_indx[mask]
        indices = np.arange(c_indx.size)[mask]
        sorting_idx = np.argsort(elements, kind='mergesort')
        ind_sorted = indices[sorting_idx]
        x = np.searchsorted(elements, range(N), side='right',
                            sorter=sorting_idx)
        # Indexes.
        rang_indx = [ind_sorted[x[i]:x[i + 1]] for i in range(N - 1)]

        # Pick a number (given by the list 'di') of random elements in
        # each range. Those are the indexes of the elements that
        # should be removed from the sub-lists.
        rem_indx = []
        for indx, num in enumerate(di):
            if rang_indx[indx].any() and len(rang_indx[indx]) >= num:
                rem_indx.append(np.random.choice(rang_indx[indx],
                                int(num), replace=False))
            else:
                rem_indx.append(rang_indx[indx])

        # Remove items from list.
        # itertools.chain() flattens the list of indexes and sorted()
        # with reverse=True inverts them so we don't change the
        # indexes of the elements in the lists after removing them.
        d_i = sorted(list(itertools.chain(*rem_indx)), reverse=True)
        # Remove those selected indexes from *all* sub-lists.
        isoch_compl = np.delete(np.asarray(isoch_binar), d_i, axis=1)

        # Remove stars from the binaries list that were removed by the
        # completeness process.
        # Sort list first, so smaller indexes are first.
        d_i.sort()
        binar_idx1 = np.setdiff1d(binar_idx0, d_i)
        # Correct indexes of stars after completeness removal, so they will
        # point to the actual binary systems.
        binar_idx = binar_idx1 - np.searchsorted(d_i, binar_idx1)
    else:
        isoch_compl, binar_idx = np.asarray(isoch_binar), []

    return isoch_compl, binar_idx


def idxFind(N, c_indx):
    """
    """
    # Reject stars in the 0th position. These are stars below the value
    # where the completeness loss starts.
    # TODO I believe the '& (c_indx < N)' part is not necessary.
    mask = (c_indx > 0) & (c_indx < N)
    # Keep those stars with indexes in the accepted magnitude range.
    elements = c_indx[mask]
    # Generate new ordered indexes for the masked stars.
    indices = np.arange(c_indx.size)[mask]
    # Store the indexes that would sort the 'elements' array.
    sorting_idx = np.argsort(elements, kind='mergesort')
    # ???
    ind_sorted = indices[sorting_idx]
    x = np.searchsorted(elements, range(N), side='right',
                        sorter=sorting_idx)

    # Indexes.
    rang_indx = [ind_sorted[x[i]:x[i + 1]] for i in range(N - 1)]

    return rang_indx


def remStars(isoch_binar, binar_idx0, d_i):
    """
    """
    # Remove those selected indexes from *all* sub-lists.
    isoch_compl = np.delete(isoch_binar, d_i, axis=1)

    # Remove stars from the binaries list that were removed by the
    # completeness process.
    # Sorted elements of 'binar_idx0' that are not in 'd_i'.
    # d_i.sort()  # <--- ??????
    binar_idx1 = np.setdiff1d(binar_idx0, d_i)
    # Correct indexes of stars after completeness removal, so they will
    # point to the actual binary systems.
    binar_idx = binar_idx1 - np.searchsorted(d_i, binar_idx1)

    return isoch_compl, binar_idx


def main_old2(isoch_binar, binar_idx0, completeness):
    '''
    Remove a number of stars according to the percentages of star loss found in
    the mag_completeness function of the luminosity module, for the real
    observation.
    '''

    # 'bin_edges' of the observed region histogram.
    bin_edges, max_indx, comp_perc = completeness

    # If stars exist in the isochrone beyond the completeness magnitude
    # level, then apply the removal of stars. Otherwise, skip it.
    if max(isoch_binar[0]) > bin_edges[max_indx]:

        # Synthetic cluster histogram.
        synth_mag_hist = np.histogram(isoch_binar[0], bin_edges)[0]
        # Number of stars that should be removed in each edge range beyond the
        # 'max_indx' magnitude index.
        di = np.around(
            (synth_mag_hist[max_indx:] -
             synth_mag_hist[max_indx] * comp_perc), 0)

        # Map indexes of elements in 'isoch_binar' to the edge ranges defined
        # beyond 'max_indx' (using the main magnitude)
        c_indx = np.searchsorted(
            bin_edges[max_indx:], isoch_binar[0], side='left')

        # Actual indexes of stars, stored in each edge range.
        rang_indx = idxFind(len(bin_edges[max_indx:]), c_indx)

        # Pick a number (given by the list 'di') of random elements in
        # each range. Those are the indexes of the elements that
        # should be removed from the sub-lists.
        rem_indx = []
        for indx, num in enumerate(di):
            if rang_indx[indx].any() and len(rang_indx[indx]) >= num:
                rem_indx.append(np.random.choice(rang_indx[indx],
                                int(num), replace=False))
            else:
                rem_indx.append(rang_indx[indx])
        # Remove items from list.
        # itertools.chain() flattens the list of indexes and sorted()
        # with reverse=True inverts them so we don't change the
        # indexes of the elements in the lists after removing them.
        # d_i = sorted(list(itertools.chain(*rem_indx)), reverse=True)
        d_i = list(itertools.chain(*rem_indx))

        # Remove those selected indexes from *all* sub-lists.
        isoch_compl, binar_idx = remStars(isoch_binar, binar_idx0, d_i)
    else:
        isoch_compl, binar_idx = isoch_binar, binar_idx0

    return isoch_compl, binar_idx


def indxRem(di, rang_indx):
    """
    Select a fixed number (given by 'di') of random indexes in 'rang_indx'.
    These correspond to the stars that will be removed in each magnitude
    range.

    Source: https://stackoverflow.com/a/46079837/1391441
    """
    lens = np.array([len(_) for _ in rang_indx])
    di0 = np.minimum(lens, di)
    invalid_mask = lens[:, None] <= np.arange(lens.max())
    # Create a 2D random array in interval [0,1) to cover the max. length of
    # subarrays.
    rand_nums = np.random.rand(len(lens), lens.max())
    # For each subarray, set the invalid places to 1.0. Get argsort for each
    # row. Those 1s corresponding to the invalid places would stay at the back
    # because there were no 1s in the original random array. Thus, we have the
    # indices array.
    rand_nums[invalid_mask] = 1
    # Slice each row of those indices array to the extent of the lengths
    # listed in di.
    shuffled_indx = np.argpartition(rand_nums, lens - 1, axis=1)

    # Start a loop and slice each subarray from 'rang_indx' using those sliced
    # indices.
    out = []
    for i, all_idx in enumerate(shuffled_indx):
        if lens[i] > 0:
            slice_idx = all_idx[:di0[i]]
            out += rang_indx[i][slice_idx].tolist()

    return np.asarray(out)


def main(isoch_binar, binar_idx0, completeness):
    '''
    Remove a number of stars according to the percentages of star loss found in
    the mag_completeness function of the luminosity module, for the real
    observation.
    '''
    # If stars exist in the isochrone beyond the completeness magnitude
    # level, then apply the removal of stars. Otherwise, skip it.

    # completeness = [bin_edges, max_indx, comp_perc]
    # 'bin_edges' of the observed region histogram.

    bin_edges, max_indx, comp_perc = completeness
    if np.max(isoch_binar[0]) > bin_edges[max_indx]:

        # Indexes of stars in 'isoch_binar' whose main magnitude
        # value falls between the ranges given.
        c_indx = np.searchsorted(bin_edges[max_indx:],
                                 isoch_binar[0], side='left')

        # Number of elements to keep in each range.
        # Reject elements in the '0' range, ie: below the smallest edge.
        # The minimum length is that of the 'comp_perc' list plus one,
        # so after removing the '0' element both lists will have the same
        # length.
        # Equivalent to np.histogram(isoch_binar[0], bin_edges)[0][max_indx:]
        count = np.bincount(c_indx, minlength=len(comp_perc) + 1)[1:]
        # Clip at '0' so there are no negative values.
        di = np.rint(count - count[0] * comp_perc).astype(int).clip(0)

        # Actual indexes of stars, stored in each edge range.
        rang_indx = idxFind(len(bin_edges[max_indx:]), c_indx)

        # Pick a number (given by the list 'di') of random elements in
        # each range. Those are the indexes of the elements that
        # should be removed from the sub-lists.
        d_i = indxRem(di, rang_indx)

        #
        isoch_compl, binar_idx = remStars(isoch_binar, binar_idx0, d_i)
    else:
        isoch_compl, binar_idx = np.asarray(isoch_binar), binar_idx0

    return isoch_compl, binar_idx


if __name__ == '__main__':

    import pickle
    import time

    with open('completeness.pickle', 'rb') as f:
        isoch_binar, binar_idx0, completeness = pickle.load(f)
    completeness[2] = np.asarray(completeness[2]) / completeness[2][0]

    N = 1000

    s = time.clock()
    for _ in range(N):
        main_orig(isoch_binar, binar_idx0, completeness)
    print(time.clock() - s)

    #
    s = time.clock()
    for _ in range(N):
        main(isoch_binar, binar_idx0, completeness)
    print(time.clock() - s)

    # import matplotlib.pyplot as plt
    # x_bars = np.arange(0, 16, 1)
    # plt.bar(x_bars, times_perc)
    # plt.xticks(x_bars,
    #     ('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13',
    #      '14', '15', '16'))
    # plt.show()


# def compl_func2(isoch_binar):
#     '''
#     Remove random stars beyond a given magnitude limit according to a
#     completeness decreasing function.
#     '''
#     import random as rd

#     # Magnitude value below the minumum magnitude where the completeness
#     # removal should start.
#     c_mags = 2.5

#     mag_data = isoch_binar[1]
#     max_mag = max(mag_data)
#     # Number of bins.
#     bins1 = int((max(mag_data) - min(mag_data)) / 0.2)

#     # Histogram of magnitude values.
#     mag_hist, bin_edg = np.histogram(mag_data, bins1)
#     # Index of maximum magnitude bin, c_mags mags below the max mag value.
#     max_indx = min(range(len(bin_edg)),
#         key=lambda i: abs(bin_edg[i] - (max_mag - c_mags)))
#     n1, p1 = mag_hist[max_indx], 100.
#     # Get completeness percentages.
#     a = rd.uniform(2., 4.)
#     # Get percentages of completeness for each mag bin, according to the core
#     # completeness function defined: 1 / (1 + np.exp(x - a))
#     comp_perc = [(1 / (1 + np.exp(i - a))) * 100.
#         for i in range(len(mag_hist[max_indx:]))]
#     # Number of stars that should be removed from each bin.
#     di = np.around((abs(mag_hist[max_indx:] - (n1 / p1) *
#         np.asarray(comp_perc))), 0)

#     # Store indexes of *all* elements in mag_data whose magnitude
#     # value falls between the ranges given.
#     c_indx = np.searchsorted(bin_edg[max_indx:], mag_data, side='left')
#     N = len(bin_edg[max_indx:])
#     mask = (c_indx > 0) & (c_indx < N)
#     elements = c_indx[mask]
#     indices = np.arange(c_indx.size)[mask]
#     sorting_idx = np.argsort(elements, kind='mergesort')
#     ind_sorted = indices[sorting_idx]
#     x = np.searchsorted(elements, range(N), side='right',
#         sorter=sorting_idx)
#     # Indexes.
#     rang_indx = [ind_sorted[x[i]:x[i + 1]] for i in range(N - 1)]

#     # Pick a number (given by the list 'di') of random elements in
#     # each range. Those are the indexes of the elements that
#     # should be removed from the three sub-lists.
#     rem_indx = []
#     for indx, num in enumerate(di):
#         if rang_indx[indx].any() and len(rang_indx[indx]) >= num:
#             rem_indx.append(np.random.choice(rang_indx[indx],
#                 int(num), replace=False))
#         else:
#             rem_indx.append(rang_indx[indx])

#     # Remove items from list.
#     # itertools.chain() flattens the list of indexes and sorted()
#     # with reverse=True inverts them so we don't change the
#     # indexes of the elements in the lists after removing them.
#     d_i = sorted(list(itertools.chain(*rem_indx)), reverse=True)
#     # Remove those selected indexes from all the sub-lists.
#     isoch_compl = np.delete(np.asarray(isoch_binar), d_i, axis=1)

#     return isoch_compl
