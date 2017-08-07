
import numpy as np
import traceback
import itertools


def rem_bad_stars(ids, x, y, mags, em, cols, ec):
    '''
    Remove stars from all lists that have too large magnitude or color
    values (or their errors) which indicates a bad photometry.
    '''
    # Set photometric range for accepted stars.
    min_lim, max_lim = -50., 50.

    # Store indexes of stars that should be removed.
    lists_arr = list(zip(*itertools.chain(mags, em, cols, ec)))
    del_indexes = [i for i, t in enumerate(lists_arr) if
                   any(e > max_lim for e in t) or any(e < min_lim for e in t)]

    # Remove stars from id list first since this are strings.
    id_clean = np.delete(np.array(ids), del_indexes)
    # Remove stars from the coordinates lists.
    x_clean, y_clean = np.delete(np.array([x, y]), del_indexes, axis=1)
    # Remove stars from the rest of the lists.
    mags_clean = np.delete(np.array(mags), del_indexes, axis=1)
    em_clean = np.delete(np.array(em), del_indexes, axis=1)
    cols_clean = np.delete(np.array(cols), del_indexes, axis=1)
    ec_clean = np.delete(np.array(ec), del_indexes, axis=1)

    return id_clean, x_clean, y_clean, mags_clean, em_clean, cols_clean,\
        ec_clean


def main(npd, id_indx, x_indx, y_indx, mag_indx, e_mag_indx, col_indx,
         e_col_indx, **kwargs):
    '''
    Get spatial and photometric data from the cluster's data file.
    '''

    data_file = npd['data_file']
    # Loads the data in 'data_file' as a list of N lists where N is the number
    # of columns. Each of the N lists contains all the data for the column.
    # If any string is found (for example 'INDEF') it is converted to 99.999.
    try:
        data = np.genfromtxt(data_file, dtype=float, filling_values=99.999,
                             unpack=True)
    except ValueError:
        print traceback.format_exc()
        raise ValueError("ERROR: could not read data input file:\n  {}\n"
                         "  Check that all rows are filled (i.e., no blank"
                         " spaces)\n  for all columns.\n".format(data_file))

    try:
        # Read coordinates data.
        x, y = data[x_indx], data[y_indx]
        # Read magnitudes.
        mags, em = [], []
        for mi, emi in zip(*[mag_indx, e_mag_indx]):
            mags.append(data[mi])
            em.append(data[emi])
        # Read colors.
        cols, ec = [], []
        for ci, eci in zip(*[col_indx, e_col_indx]):
            cols.append(data[ci])
            ec.append(data[eci])

        # Now read IDs as strings. Do this separately so numeric IDs are not
        # converted into floats by np.genfromtxt. I.e.: 190 --> 190.0
        data = np.genfromtxt(data_file, dtype=str, unpack=True)
        ids = data[id_indx]
        n_old = len(ids)
    except IndexError:
        raise IndexError("ERROR: data input file:\n  {}\n  contains "
                         "fewer columns than those given "
                         "in 'params_input.dat'.".format(data_file))

    # If any magnitude or color value (or their errors) is too large, discard
    # that star.
    ids, x, y, mags, em, cols, ec = rem_bad_stars(
        ids, x, y, mags, em, cols, ec)

    # Check if array came back empty after removal of stars with
    # bad photometry.
    if not x.size:
        raise ValueError("ERROR: no stars left after removal of those "
                         "with\n large mag/color or error values. Check "
                         "input file.")

    # Check if the range of any coordinate column is zero.
    data_names = ['x_coords', 'y_coords']
    for i, dat_lst in enumerate([x, y]):
        if min(dat_lst) == max(dat_lst):
            raise ValueError("ERROR: the range for the '{}' column\n"
                             "is zero. Check the input data format.".format(
                                 data_names[i]))
    # Check if the range of any photometric column is zero.
    data_names = ['magnitude', 'color']
    for i, dat_lst in enumerate([mags, cols]):
        for mc in dat_lst:
            if min(mc) == max(mc):
                raise ValueError(
                    "ERROR: the range for {} column {} is\nzero."
                    " Check the input data format.".format(data_names[i], i))

    print 'Data obtained from input file (N_stars: %d).' % len(ids)
    frac_reject = (float(n_old) - len(ids)) / float(n_old)
    if frac_reject > 0.05:
        print("  WARNING: {:.0f}% of stars in cluster's file were"
              " rejected.".format(100. * frac_reject))

    # Create cluster's data dictionary.
    cld = {'ids': ids, 'x': x, 'y': y, 'mags': mags, 'em': em, 'cols': cols,
           'ec': ec}
    return cld
