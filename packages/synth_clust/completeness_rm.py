
import numpy as np
import itertools


def main(isoch_binar, completeness):
    '''
    Remove a number of stars according to the percentages of star loss found in
    the mag_completeness function of the luminosity module, for the real
    observation.
    '''
    # If stars exist in the isochrone beyond the completeness magnitude
    # level, then apply the removal of stars. Otherwise, skip it.
    # completeness = [max_mag, bin_edges, max_indx, comp_perc]
    if max(isoch_binar[0]) > completeness[1][completeness[2]]:

        # Get histogram. completeness[1] = bin_edges of the observed
        # region histogram.
        synth_mag_hist, bin_edges = np.histogram(isoch_binar[0],
                                                 completeness[1])
        pi = completeness[3]
        n1, p1 = synth_mag_hist[completeness[2]], pi[0]
        di = np.around((synth_mag_hist[completeness[2]:] -
                        (n1 / p1) * np.asarray(pi)), 0)

        # Store indexes of *all* elements in isoch_binar whose main magnitude
        # value falls between the ranges given.
        c_indx = np.searchsorted(completeness[1][completeness[2]:],
                                 isoch_binar[0], side='left')
        N = len(completeness[1][completeness[2]:])
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
    else:
        isoch_compl = np.asarray(isoch_binar)

    return isoch_compl


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
