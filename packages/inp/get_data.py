
import numpy as np
import traceback
import input_params as g


def rem_bad_stars(id_star, x_data, y_data, mag_data, e_mag, col1_data,
                  e_col1):
    '''
    Remove stars from all lists that have too large magnitude or color
    values (or their errors) which indicates a bad photometry.
    '''
    # Set photometric range for accepted stars.
    min_lim, max_lim = -50., 50.

    # Store indexes of stars that should be removed.
    lists_arr = zip(mag_data, e_mag, col1_data, e_col1)
    del_indexes = [i for i, t in enumerate(lists_arr) if
                   any(e > max_lim for e in t) or any(e < min_lim for e in t)]

    # Remove stars from id list first since this are strings.
    id_clean = np.delete(np.array(id_star), del_indexes)
    # Remove stars from the rest of the lists simultaneously.
    clean_array = np.delete(np.array([x_data, y_data, mag_data, e_mag,
                            col1_data, e_col1]), del_indexes, axis=1)

    return id_clean, clean_array


def main(data_file):
    '''
    Get spatial and photometric data from the cluster's data file.
    '''

    # Read indexes from input params.
    id_inx, x_inx, y_inx, m_inx, em_inx, c_inx, ec_inx = g.gd_params[:-1]

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
        x_data, y_data, mag_data, e_mag, col1_data, e_col1 = \
            data[x_inx], data[y_inx], data[m_inx], data[em_inx], data[c_inx],\
            data[ec_inx]

        # Now read IDs as strings. Do this separately so numeric IDs are not
        # converted into floats by np.genfromtxt. I.e.: 190 --> 190.0
        data = np.genfromtxt(data_file, dtype=str, unpack=True)
        id_star = data[id_inx]
        n_old = len(id_star)
    except IndexError:
        raise IndexError("ERROR: data input file:\n  {}\n  contains "
                         "fewer columns than those given "
                         "in 'params_input.dat'.".format(data_file))

    # If any mag or color value (or their errors) is too large, discard
    # that star.
    id_star, [x_data, y_data, mag_data, e_mag, col1_data, e_col1] = \
        rem_bad_stars(id_star, x_data, y_data, mag_data, e_mag, col1_data,
                      e_col1)

    data_names = ['x_coords', 'y_coords', 'magnitudes', 'color']
    # Check read coordinates, and photometry.
    for i, dat_lst in enumerate([x_data, y_data, mag_data, col1_data]):
        # Check if array came back empty after removal of stars with
        # bad photometry.
        if not dat_lst.size:
            raise ValueError("ERROR: no stars left after removal of those "
                             "with\n large mag/color or error values. Check "
                             "input file.")
        # Check if the range of any photometric column, excluding errors,
        # is none.
        if min(dat_lst) == max(dat_lst):
            raise ValueError("ERROR: the range for the '{}' column\n"
                             "  is zero. Check the input data format.".format(
                                 data_names[i]))

    print 'Data obtained from input file (N_stars: %d).' % len(id_star)
    frac_reject = (float(n_old) - len(id_star)) / float(n_old)
    if frac_reject > 0.05:
        print("  WARNING: {:.0f}% of stars in file were"
              " rejected.".format(100. * frac_reject))

    return id_star, x_data, y_data, mag_data, e_mag, col1_data, e_col1
