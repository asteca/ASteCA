
import numpy as np
from astropy.io import ascii
from collections import defaultdict, Iterable
import operator
import copy


def main(npd, read_mode, id_col, x_col, y_col, mag_col, e_mag_col,
         col_col, e_col_col, plx_col, e_plx_col, pmx_col, e_pmx_col,
         pmy_col, e_pmy_col, rv_col, e_rv_col, **kwargs):
    '''
    Read all data from the cluster's data file.
    '''

    data_file = npd['data_file']
    try:
        fill_msk = [
            ('', '0'), ('INDEF', '0'), ('9999.99', '0'), ('99.999', '0')]
        # Store IDs as strings.
        if read_mode == 'num':
            data = ascii.read(
                data_file, fill_values=fill_msk,
                converters={id_col: [ascii.convert_numpy(np.str)]},
                format='no_header')
        else:
            data = ascii.read(
                data_file, fill_values=fill_msk,
                converters={id_col: [ascii.convert_numpy(np.str)]})

        # Arrange column names in the proper order and shape.
        col_names = [
            id_col, x_col, y_col, mag_col, e_mag_col, col_col, e_col_col,
            [plx_col, pmx_col, pmy_col, rv_col],
            [e_plx_col, e_pmx_col, e_pmy_col, e_rv_col]]

        # Remove not wanted columns *before* removing rows with 'nan' values
        # (otherwise columns that should not be read will influence the row
        # removal).
        col_names_keep = list(flatten(col_names))
        for col in data.columns:
            if col not in col_names_keep:
                data.remove_column(col)

        # Check if there are any masked elements in the data table.
        masked_elems = 0
        for col in data.columns:
            try:
                masked_elems += data[col].mask.nonzero()[0].sum()
            except AttributeError:
                pass

        # Remove all rows with at least one masked element.
        flag_data_eq = False
        if masked_elems > 0:
            data_compl = data[reduce(
                operator.and_, [~data[col].mask for col in data.columns])]
        else:
            # If there where no elements to mask, there were no bad values.
            data_compl = copy.deepcopy(data)
            flag_data_eq = True

        # Change masked elements with 'nan' values, in place.
        fill_cols(data)

    except ascii.InconsistentTableError:
        raise ValueError("ERROR: could not read data input file:\n  {}\n"
                         "  Check that all rows are filled (i.e., no blank"
                         " spaces)\n  for all columns.\n".format(data_file))

    # Create cluster's dictionary with the *incomplete* data.
    ids, x, y, mags, cols, kine, em, ec, ek = dataCols(
        data_file, data, col_names)
    perc_compl_check(ids, mags, cols)
    cld_i = {'ids': ids, 'x': x, 'y': y, 'mags': mags, 'em': em,
             'cols': cols, 'ec': ec, 'kine': kine, 'ek': ek}

    # Create cluster's dictionary with the *complete* data.
    ids, x, y, mags, cols, kine, em, ec, ek = dataCols(
        data_file, data_compl, col_names)
    cld_c = {'ids': ids, 'x': x, 'y': y, 'mags': mags, 'em': em,
             'cols': cols, 'ec': ec, 'kine': kine, 'ek': ek}

    clp = {'flag_data_eq': flag_data_eq}

    return cld_i, cld_c, clp


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


def dataCols(data_file, data, col_names):
    """
    Separate data into appropriate columns.
    """
    try:
        ids = np.array(data[col_names[0]])
        dups = list(list_duplicates(ids))
        if dups:
            print("ERROR: duplicated IDs found in data file:")
            for dup in dups:
                print("  ID '{}' found in lines: {}".format(dup[0],
                      ", ".join(dup[1])))
            raise ValueError("Duplicated IDs found.")

        # Read coordinates data.
        x, y = np.array(data[col_names[1]]), np.array(data[col_names[2]])

        # Read magnitudes
        mags, em = [], []
        for mi, emi in zip(*[col_names[3], col_names[4]]):
            mags.append(np.array(data[mi]))
            em.append(np.array(data[emi]))
        # Read colors.
        cols, ec = [], []
        for ci, eci in zip(*[col_names[5], col_names[6]]):
            cols.append(np.array(data[ci]))
            ec.append(np.array(data[eci]))
        mags, cols, em, ec = np.array(mags), np.array(cols), np.array(em),\
            np.array(ec)

        # Read PMs, parallax, RV.
        kine, ek = [], []
        for ki, eki in zip(*[col_names[7], col_names[8]]):
            if ki is not False:
                kine.append(np.array(data[ki]))
                ek.append(np.array(data[eki]))
            else:
                kine.append(np.full(ids.size, np.nan))   # NOPE
                ek.append(np.full(ids.size, np.nan))
        kine, ek = np.array(kine), np.array(ek)

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

    return ids, x, y, mags, cols, kine, em, ec, ek


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


def perc_compl_check(ids, mags, cols):
    """
    Check percentage of complete data.
    """
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


if __name__ == '__main__':
    main()
