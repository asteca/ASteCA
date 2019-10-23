
import numpy as np
import itertools
import re
from . import isochs_format
from .. import update_progress


def main(met_files, ages_strs, evol_track, CMD_extra_pars, all_syst_filters):
    """
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

    """

    isoch_list, extra_pars = [], []

    # Equal for all photometric systems in all sets of evolutionary tracks.
    age_format = isochs_format.age_format(evol_track)
    # Depends on the evolutionary track set.
    line_start = isochs_format.line_start_format(evol_track)

    # For each group of metallicity files (representing a single metallicity
    # value) in all photometric systems defined.
    met_fls_photsysts = list(zip(*met_files))
    N_met_files = len(met_fls_photsysts)
    for i_met, met_fls in enumerate(met_fls_photsysts):

        # Iterate through the metallicity files stored, one per system.
        filters, mini = [], []
        for j, met_f in enumerate(met_fls):

            # Depends on the photometric system analyzed.
            l_s = isochs_format.read_line_start(met_f, evol_track, line_start)

            # Column indexes for all the filters defined in this system.
            # TODO the extra params consume a lot of memory and are not used
            # (for now?), hence the [0] to pass only the 'M_ini' index.
            uniq_fltrs = list(all_syst_filters[j][1:]) + [CMD_extra_pars[0]]

            # Column numbers for the filters defined in this system.
            ids = isochs_format.common_ids(evol_track, uniq_fltrs, l_s)

            # Store list holding the filter columns with the same metallicity.
            filters_mini = readMetFile(
                evol_track, met_f, ages_strs, line_start, age_format, ids)

            # Separate filter data from the initial mass, assuming the latter
            # is always positioned last.
            filters.append([_[:-1] for _ in filters_mini])
            # The extra '[]' is to match the shape of 'filters'.
            mini.append([[_[-1]] for _ in filters_mini])

        # Store data for this metallicity value. Re-arrange as (#ages, #cols)
        # before appending.
        isoch_list.append([
            list(itertools.chain(*_)) for _ in zip(*filters)])
        extra_pars.append([
            list(itertools.chain(*_)) for _ in zip(*mini)])

        update_progress.updt(N_met_files, i_met + 1)

    return isoch_list, extra_pars


def readMetFile(
        evol_track, met_f, ages_strs, line_start, age_format, column_ids):
    """
    Read a given metallicity file.
    """
    if evol_track[:3] == 'PAR':
        return readCMDFile(
            met_f, ages_strs, line_start, age_format, column_ids)
    else:
        # TODO in place for #275
        pass


def readCMDFile(met_f, ages_strs, line_start, age_format, column_ids):
    """
    Read a given metallicity file from the CMD service, and return the
    isochrones for the ages within the age range.
    """

    # This cant not be a numpy array since different ages have different number
    # of stars, and an array can not have an irregular shape.
    metal_isoch = []
    # Open the metallicity file.
    with open(met_f, mode="r") as f_iso:
        content = f_iso.readlines()

        # Identify the end of the last isochrone as the first line in reverse
        # order that is neither a comment nor a newline.
        for i, rl in enumerate(reversed(content)):
            if not rl.startswith('#') and rl != '\n':
                i_end = len(content) - i
                break

        # Identify positions of all isochrone starting lines.
        idx = [
            i for i, line in enumerate(content) if line.startswith(line_start)]

        # Age value for all isochrones in file.
        all_ages = [re.findall(age_format, content[i - 1])[0] for i in idx]

        def appendIsoch(a, b):
            """Extract and format isochrone data."""
            block = content[a:b]
            block = np.array([list(map(float, _.split())) for _ in block]).T
            return block[column_ids].tolist()

        for i, age in enumerate(ages_strs):
            # The isochrone stars just after this index, hence the '+1'.
            i1 = idx[all_ages.index(age)] + 1
            try:
                # Starting line of the next isochrone is the end of this one.
                i2 = idx[all_ages.index(age) + 1] - 1
            except IndexError:
                # Reached the final isochrone in the file.
                i2 = i_end

            metal_isoch.append(appendIsoch(i1, i2))

    return metal_isoch
