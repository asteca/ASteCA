
import sys
import traceback
import numpy as np
import itertools
import re
from . import isochs_format
from .. import update_progress


def main(met_f_filter, age_values, evol_track, all_syst_filters):
    '''
    Stores the available isochrones of different metallicities and
    ages, according to the ranges given to these parameters.

    Returns: isoch_list, extra_pars
    where:
    isoch_list[i][j] --> i: metallicity index ; j: age index
    extra_pars[i][j] --> i: metallicity index ; j: age index

    These lists store the magnitudes and extra data for each isochrone and
    each metallicity value:

    isoch_list = [metal_1, ..., metal_M]
    metal_i = [isoch_i1, ..., isoch_iN]
    isoch_ij = [filter1, filter2, filter3, ...]

    extra_pars = [metal_1, ..., metal_M]
    metal_i = [isoch_i1, ..., isoch_iN]
    isoch_ij = [M_ini, M_act, logL/Lo, logTe, logG, mbol]

    Where 'filterX' runs through all filters defined, for all the photometric
    systems in use. The order in which they are stored follows the order
    of the 'all_syst_filters' tuple, with its first elements removed (since
    they indicate the photometric system), and flattened.

    '''

    isoch_list, extra_pars = [], []

    # Equal for all photometric systems in all sets of evolutionary tracks.
    age_format = isochs_format.cmd_age_format()
    # Depends on the evolutionary track set.
    line_start = isochs_format.cmd_line_start_format(evol_track)

    # Store here the number of age values defined in each file, for checking.
    numb_age_values = []
    # For each group of metallicity files (representing a single metallicity
    # value) in all photometric systems defined.
    met_fls_photsysts = list(zip(*met_f_filter))
    N_met_files = len(met_fls_photsysts)
    for i_met, met_fls in enumerate(met_fls_photsysts):

        # Iterate through the metallicity files stored, one per system.
        all_systs = []
        for j, met_f in enumerate(met_fls):

            metal_isoch = filters_and_extra_pars(
                evol_track, age_values, all_syst_filters, age_format,
                line_start, met_f, j)
            # Store list holding all the isochrones with the same metallicity
            # in the final isochrone list.
            all_systs.append(metal_isoch)
            numb_age_values.append([met_f, len(metal_isoch)])

        # Store data for this metallicity value. Re-arrange first.
        all_systs_r = [list(itertools.chain(*_)) for _ in zip(*all_systs)]
        isoch_list.append(all_systs_r)

        # The extra isochrone parameters are equal across photometric
        # systems, for a given metallicity and age. Thus, we read their
        # values from the *first system defined*, for this metallicity value
        # and all the ages defined, and store it in a separate list with the
        # same order as the 'isoch_list' array.
        e_pars = filters_and_extra_pars(
            evol_track, age_values, all_syst_filters, age_format, line_start,
            met_fls[0], -1)
        # Store in list.
        extra_pars.append(e_pars)

        update_progress.updt(N_met_files, i_met + 1)

    # Check that all metallicity files contain the same number of ages.
    if not checkEqual([len(_) for _ in isoch_list]):
        print("")
        for met_f, N in numb_age_values:
            print(" {} ages in: {}".format(N, met_f.split('isochrones/')[1]))
        sys.exit("ERROR: not all metallicity files contain the same number\n"
                 "of defined ages.")

    return isoch_list, extra_pars


def filters_and_extra_pars(
    evol_track, age_values, all_syst_filters, age_format, line_start,
        met_f, j):
    """
    Read the data from all the filters in all the photometric systems defined,
    and also the extra parameters for each metallicity and age value (masses,
    effective temperatures, etc.)
    """
    # Depends on the photometric system analyzed.
    l_s = isochs_format.read_line_start(met_f, line_start)
    if j >= 0:
        identif = 'filters data'
        # Column indexes for all the filters defined in this system.
        uniq_fltrs = all_syst_filters[j][1:]
        # Column numbers for the filters defined in this system.
        ids = isochs_format.girardi_filters_ids(l_s, uniq_fltrs)
    else:
        identif = 'extra parameters'
        # Column numbers for the extra parameters for this metallicity.
        ids = isochs_format.cmd_common_ids(evol_track, l_s)

    try:
        met_f_ages = readCMDFile(
            met_f, age_values, line_start, age_format, ids)
    except Exception:
        print(traceback.format_exc())
        sys.exit("Error reading {} from \n'{}'\n"
                 "metallicity file.".format(identif, met_f))

    return met_f_ages


def readCMDFile(met_f, age_values, line_start, age_format, column_ids):
    '''
    Read a given metallicity file from the CMD service, and return the
    isochrones for the ages within the age range.
    '''
    metal_isoch = []

    # Open the metallicity file.
    with open(met_f, mode="r") as f_iso:
        content = f_iso.readlines()

        # identify positions of all isochrone starting lines.
        idx = [
            i for i, line in enumerate(content) if line.startswith(line_start)]

        # Age value for all isochrones in file.
        all_ages = [
            float(re.findall(age_format, content[i - 1])[0]) for i in idx]
        # TODO hardcoded rounding. Must be modified in accordance with
        # 'met_ages_values.CMDAges()'.
        ages = np.around(np.log10(all_ages), 4)

        def appendIsoch(a, b):
            """Extract and format isochrone data."""
            block = content[a:b]
            block = np.array([map(float, _.split()) for _ in block]).T
            return block[column_ids].tolist()

        N = len(idx)
        for k, i in enumerate(idx):
            age = ages[k]
            # If age value falls inside the given range..
            if age in age_values:
                if k < N - 1:
                    metal_isoch.append(appendIsoch(i + 1, idx[k + 1] - 1))
                else:
                    # Save the last isochrone.
                    metal_isoch.append(appendIsoch(idx[-1] + 1, -1))

    return metal_isoch


def checkEqual(lst):
    return lst[1:] == lst[:-1]
