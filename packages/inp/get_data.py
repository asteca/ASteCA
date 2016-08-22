
import numpy as np
import traceback


def rem_bad_stars(ids, x, y, mags, em, cols, ec):
    '''
    Remove stars from all lists that have too large magnitude or color
    values (or their errors) which indicates a bad photometry.
    '''
    # Set photometric range for accepted stars.
    min_lim, max_lim = -50., 50.

    # Store indexes of stars that should be removed.
    lists_arr = zip(mags, em, cols, ec)
    del_indexes = [i for i, t in enumerate(lists_arr) if
                   any(e > max_lim for e in t) or any(e < min_lim for e in t)]

    # Remove stars from id list first since this are strings.
    id_clean = np.delete(np.array(ids), del_indexes)
    # Remove stars from the rest of the lists simultaneously.
    clean_array = np.delete(np.array([x, y, mags, em,
                            cols, ec]), del_indexes, axis=1)

    return id_clean, clean_array


def main(npd, gd_params, **kwargs):
    '''
    Get spatial and photometric data from the cluster's data file.
    '''

    data_file = npd['data_file']
    # Read indexes from input parameters.
    id_inx, x_inx, y_inx, m_inx, em_inx, c_inx, ec_inx = gd_params[:-1]

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
        # Read data columns, except IDs.
        x, y, mags, em, cols, ec = \
            data[x_inx], data[y_inx], data[m_inx], data[em_inx], data[c_inx],\
            data[ec_inx]

        # Now read IDs as strings. Do this separately so numeric IDs are not
        # converted into floats by np.genfromtxt. I.e.: 190 --> 190.0
        data = np.genfromtxt(data_file, dtype=str, unpack=True)
        ids = data[id_inx]
        n_old = len(ids)
    except IndexError:
        raise IndexError("ERROR: data input file:\n  {}\n  contains "
                         "fewer columns than those given "
                         "in 'params_input.dat'.".format(data_file))

    # If any mag or color value (or their errors) is too large, discard
    # that star.
    ids, [x, y, mags, em, cols, ec] = rem_bad_stars(
        ids, x, y, mags, em, cols, ec)

    data_names = ['x_coords', 'y_coords', 'magnitudes', 'color']
    # Check read coordinates, and photometry.
    for i, dat_lst in enumerate([x, y, mags, cols]):
        # Check if array came back empty after removal of stars with
        # bad photometry.
        if not dat_lst.size:
            raise ValueError("ERROR: no stars left after removal of those "
                             "with\n  large mag/color or error values. Check "
                             "input file.")
        # Check if the range of any photometric column, excluding errors,
        # is none.
        if min(dat_lst) == max(dat_lst):
            raise ValueError("ERROR: the range for the '{}' column\n"
                             "  is zero. Check the input data format.".format(
                                 data_names[i]))

    print 'Data obtained from input file (N_stars: %d).' % len(ids)
    frac_reject = (float(n_old) - len(ids)) / float(n_old)
    if frac_reject > 0.05:
        print("  WARNING: {:.0f}% of stars in file were"
              " rejected.".format(100. * frac_reject))

    # Create cluster's data dictionary.
    cld = {'ids': ids, 'x': x, 'y': y, 'mags': mags, 'em': em, 'cols': cols,
           'ec': ec}
    return cld
