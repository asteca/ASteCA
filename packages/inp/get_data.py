
import numpy as np
from astropy.io import ascii
from collections import defaultdict, Iterable
import operator
import copy


def list_duplicates(seq):
    """
    Find and report duplicates in list.

    Source: https://stackoverflow.com/a/5419576/1391441
    """
    tally = defaultdict(list)
    for i, item in enumerate(seq):
        tally[item].append(i)
    dups = ((key, map(str, locs)) for key, locs in tally.items()
            if len(locs) > 1)
    return dups


def fill_cols(tbl, fill=np.nan, kind='f'):
    """
    In-place fill of 'tbl' columns which have dtype 'kind' with 'fill' value.

    Source: https://stackoverflow.com/a/50164954/1391441
    """
    for col in tbl.itercols():
        try:
            if col.dtype.kind == kind:
                col[...] = col.filled(fill)
        except AttributeError:
            # If there where no elements to mask, the type is 'Column' and
            # not 'MaskedColumn', so .filled() is not a valid attribute.
            pass


def dataCols(data_file, data, id_indx, x_indx, y_indx, mag_indx, e_mag_indx,
             col_indx, e_col_indx):
    """
    Separate data into appropriate columns.
    """
    try:
        # Read coordinates data.
        x, y = np.array(data.columns[x_indx]), np.array(data.columns[y_indx])
        # Read magnitudes.
        mags, em = [], []
        for mi, emi in zip(*[mag_indx, e_mag_indx]):
            mags.append(np.array(data.columns[mi]))
            em.append(np.array(data.columns[emi]))
        # Read colors.
        cols, ec = [], []
        for ci, eci in zip(*[col_indx, e_col_indx]):
            cols.append(np.array(data.columns[ci]))
            ec.append(np.array(data.columns[eci]))
        mags, cols, em, ec = np.array(mags), np.array(cols), np.array(em),\
            np.array(ec)

        ids = np.array(data.columns[id_indx])
        dups = list(list_duplicates(ids))
        if dups:
            print("ERROR: duplicated IDs found in data file:")
            for dup in dups:
                print("  ID '{}' found in lines: {}".format(dup[0],
                      ", ".join(dup[1])))
            raise ValueError("Duplicated IDs found.")

    except IndexError:
        raise IndexError("ERROR: data input file:\n  {}\n  contains "
                         "fewer columns than those given "
                         "in 'params_input.dat'.".format(data_file))

    # Check if the range of any coordinate column is zero.
    data_names = ['x_coords', 'y_coords']
    for i, dat_lst in enumerate([x, y]):
        if np.min(dat_lst) == np.max(dat_lst):
            raise ValueError("ERROR: the range for the '{}' column\n"
                             "is zero. Check the input data format.".format(
                                 data_names[i]))
    # Check if the range of any photometric column is zero.
    data_names = ['magnitude', 'color']
    for i, dat_lst in enumerate([mags, cols]):
        for mc in dat_lst:
            if np.min(mc) == np.max(mc):
                raise ValueError(
                    "ERROR: the range for {} column {} is\nzero."
                    " Check the input data format.".format(data_names[i], i))

    return ids, x, y, mags, cols, em, ec


def flatten(l):
    """
    Source: https://stackoverflow.com/a/2158532/1391441
    """
    for el in l:
        if isinstance(el, Iterable) and not isinstance(el, basestring):
            for sub in flatten(el):
                yield sub
        else:
            yield el


def main(npd, id_indx, x_indx, y_indx, mag_indx, e_mag_indx, col_indx,
         e_col_indx, **kwargs):
    '''
    Read all data from the cluster's data file.
    '''

    data_file = npd['data_file']
    try:
        # Name of IDs column (the '+ 1' is because astropy Tables' first column
        # is named 'col1', not 'col0'). Store IDs as strings.
        id_colname = 'col' + str(id_indx + 1)
        data = ascii.read(
            data_file, fill_values=[('INDEF', '0'), ('9999.99', '0')],
            converters={id_colname: [ascii.convert_numpy(np.str)]},
            format='no_header')

        # Remove not wanted columns.
        col_ids = flatten([
            id_indx, x_indx, y_indx, mag_indx, e_mag_indx, col_indx,
            e_col_indx])
        col_names_keep = ['col' + str(_ + 1) for _ in col_ids]
        for col in data.columns:
            if col not in col_names_keep:
                data.remove_column(col)

        # Remove all rows with at least one masked element.
        try:
            data_compl = data[reduce(
                operator.and_, [~data[col].mask for col in data.columns])]
        except AttributeError:
            # If there where no elements to mask, there were no bad values.
            data_compl = copy.deepcopy(data)

        # Change masked elements with 'nan' values, in place.
        fill_cols(data)

    except ascii.InconsistentTableError:
        raise ValueError("ERROR: could not read data input file:\n  {}\n"
                         "  Check that all rows are filled (i.e., no blank"
                         " spaces)\n  for all columns.\n".format(data_file))

    # Create cluster's dictionary with the *incomplete* data.
    ids, x, y, mags, cols, em, ec = dataCols(
        data_file, data, id_indx, x_indx, y_indx, mag_indx, e_mag_indx,
        col_indx, e_col_indx)
    # Check percentage of complete data.
    N_ids, m_size, c_size = ids.size, [], []
    for m in mags:
        m_size.append(np.count_nonzero(~np.isnan(m)))
    for c in cols:
        c_size.append(np.count_nonzero(~np.isnan(c)))
    N_min = min(m_size + c_size)
    print('Data lines in input file (N_stars: {}).'.format(N_ids))
    frac_reject = (N_ids - N_min) / float(N_ids)
    if frac_reject > 0.05:
        print("  WARNING: {:.0f}% of stars in input file contain\n"
              "  invalid photometric data.".format(100. * frac_reject))
    cld_i = {'ids': ids, 'x': x, 'y': y, 'mags': mags, 'em': em,
             'cols': cols, 'ec': ec}

    # Create cluster's dictionary with the *complete* data.
    ids, x, y, mags, cols, em, ec = dataCols(
        data_file, data_compl, id_indx, x_indx, y_indx, mag_indx, e_mag_indx,
        col_indx, e_col_indx)
    cld_c = {'ids': ids, 'x': x, 'y': y, 'mags': mags, 'em': em,
             'cols': cols, 'ec': ec}

    return cld_i, cld_c
