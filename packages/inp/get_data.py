
import numpy as np
import pandas as pd
# import operator
# from functools import reduce
from ..aux_funcs import flatten, list_duplicates


def main(
    npd, separators, sep, id_col, x_col, y_col, mag_col, e_mag_col, col_col,
    e_col_col, plx_col, e_plx_col, pmx_col, e_pmx_col, pmy_col, e_pmy_col,
        **kwargs):
    """
    Read data from the cluster's input file.
    """

    # Arrange column names in the proper order and shape.
    col_names = [
        id_col, x_col, y_col, mag_col, e_mag_col, col_col, e_col_col,
        [plx_col, pmx_col, pmy_col],
        [e_plx_col, e_pmx_col, e_pmy_col]]
    # Remove not wanted columns.
    col_names_keep = list(filter(bool, list(flatten(col_names))))

    # Read cluster's data file
    data_file = npd['data_file']
    data = readDataFile(col_names_keep, data_file, separators[sep])
    N_all = len(data)

    # Remove rows that contain missing data (excluding IDs column)
    # data = dataClean(data, col_names_keep[1:])
    data = data.dropna()

    N_final = len(data)
    print("Stars read from input file: N={}".format(N_all))
    if N_final != N_all:
        frac_reject = 1. - (float(N_final) / N_all)
        print(
            ("  WARNING: N={} ({:.1f}%) stars with no valid data were"
             + " removed").format(N_all - N_final, 100. * frac_reject))

    # Create cluster's dictionary
    ids, x, y, mags, cols, kine, em, ec, ek = dataCols(
        data_file, data, col_names)
    cld = {'ids': ids, 'x': x, 'y': y, 'mags': mags, 'em': em,
           'cols': cols, 'ec': ec, 'kine': kine, 'ek': ek}
    # Initiate empty parameters dictionary
    clp = {}

    return cld, clp


def readDataFile(
    col_names_keep, data_file, sep,
        nanvals=('INDEF', 'NAN', 'NaN', '--', 'nan')):
    """
    Read input data file.
    """
    id_col = col_names_keep[0]

    dtypes = {id_col: 'str'}
    for col in col_names_keep[1:]:
        dtypes[col] = np.float64
    data = pd.read_csv(
        data_file, usecols=col_names_keep, index_col=False,
        dtype=dtypes, na_values=nanvals, sep=sep)

    return data


# def dataClean(data, col_lst):
#     """
#     Remove rows with invalid/missing data
#     """
#     # Define data columns.
#     data_cols = list(flatten(col_lst))
#     # Check if there are any masked elements in the data.
#     masked_elems = 0
#     cmsk = []
#     for col in data_cols:
#         # Catch "AttributeError: 'Column' object has no attribute 'mask'"
#         # if column is not masked.
#         try:
#             masked_elems += data[col].mask.sum()
#             cmsk.append(~data[col].mask)
#         except AttributeError:
#             pass

#     # Remove rows with at least one masked element.
#     if masked_elems > 0:
#         data = data[reduce(operator.and_, cmsk)]

#     return data


def dataCols(data_file, data, col_names):
    """
    Separate data into appropriate columns.
    """
    try:
        ids = np.array(data[col_names[0]])
        dups = list(list_duplicates(ids))
        if dups:
            # print("ERROR: duplicated IDs found in data file:")
            # for dup in dups:
            #     print("  ID '{}' found in lines: {}".format(
            #         dup[0], ", ".join(dup[1])))
            raise ValueError("ERROR: duplicated IDs found in data file")

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
                kine.append(np.full(ids.size, np.nan))
                ek.append(np.full(ids.size, np.nan))
        kine, ek = np.array(kine), np.array(ek)

    except IndexError:
        raise IndexError("ERROR: data input file:\n  {}\n  contains "
                         "fewer columns than those given "
                         "in 'asteca.ini'".format(data_file))

    # Check if the range of any coordinate column is zero.
    data_names = ['x_coords', 'y_coords']
    for i, dat_lst in enumerate([x, y]):
        if np.nanmin(dat_lst) == np.nanmax(dat_lst):
            raise ValueError("ERROR: the range for the '{}' column\n"
                             "is zero. Check the input data format".format(
                                 data_names[i]))
    # Check if the range of any photometric column is zero.
    data_names = ['magnitude', 'color']
    for i, dat_lst in enumerate([mags, cols]):
        for mc in dat_lst:
            if np.nanmin(mc) == np.nanmax(mc):
                raise ValueError(
                    "ERROR: the range for {} column {} is\nzero."
                    " Check the input data format".format(data_names[i], i))

    return ids, x, y, mags, cols, kine, em, ec, ek


if __name__ == '__main__':
    main()
