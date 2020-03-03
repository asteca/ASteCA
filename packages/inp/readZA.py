
import sys
import os
import numpy as np
import re
from . import isochs_format


def main(fundam_params, iso_paths, evol_track, **kwargs):
    """
    Given the paths defined in 'iso_paths' for all the photometric systems:

    1. Read paths to all the available metallicity files, and their 'z' values.
    2. Find a subset of common 'z' values across photometric systems.
    3. Check that its range contains the values in the 'params_input' file.
    4. Filter the 'z' vales to keep only those in the defined range.
    5. Read all the age values in each accepted metallicity file.
    6. Find a subset of common ages across metallicity files and phot systems.
    7. Convert to log10(age) and sort.
    7. Check that its range contains the values in the 'params_input' file.
    8. Filter the 'log(age)' vales to keep only those in the defined range.

    **IMPORTANT**
    This process assumes that the metallicity files are generated by the
    exPadova-2 script and are thus properly formatted.

    readAgevals() is the inly function that depends on the evolutionary tracks
    service (only CMD supported for now)

    Return a list of filtered metallicity files and their values, and log ages
    and their associated identifying strings.
    """
    z_rng, a_rng = fundam_params[:2]

    # Read names and values of all the metallicity files stored in the
    # isochrones paths given.
    met_vals_all, met_files = readMetvals(iso_paths)

    # Check that a non-zero subset of z values exists across photometric
    # systems.
    met_vals_all, met_files = findMetSubset(met_vals_all, met_files)

    # Check that the z range exists in the available metallicity values.
    z_rng = checkMetAge('z', z_rng, met_vals_all)

    # Filter the metallicity files by the defined z range.
    met_vals_all, met_files = filterMetvals(z_rng, met_vals_all, met_files)

    # Read all ages in all the metallicity files defined.
    age_vals_all = readAgevals(met_files, evol_track)

    # Check that a non-zero subset of log(age) values exists across
    # photometric systems.
    ages_common = findAgeSubset(age_vals_all, iso_paths)

    # To log10() and sort.
    age_vals_all, ages_strs = logAges(ages_common)

    # Check that the log(age) range exists in the available log(age)
    # values.
    a_rng = checkMetAge('log(age)', a_rng, age_vals_all)

    age_vals_all, ages_strs = filterAgevals(a_rng, age_vals_all, ages_strs)

    return met_files, met_vals_all, age_vals_all, ages_strs


def readMetvals(iso_paths):
    """
    Read the names of *all* the metallicity files stored in each isochrones'
    path, and store them along with the z values they represent.
    """

    met_vals_all, met_files = [], []
    # For each photometric system used.
    for iso_path in iso_paths:

        # List all metallicity files in folder.
        metal_files = []
        for f in os.listdir(iso_path):
            # Remove hidden files and folders possibly present in folder.
            if os.path.isfile(os.path.join(iso_path, f))\
                    and not f.startswith('.'):
                metal_files.append(f)

        # Iterate in order through all the metallicity files stored for the
        # selected set of isochrones.
        met_vals, met_fls = [], []
        for met_file in metal_files:
            # Replace underscores in file names with decimal points.
            met_vals.append(float(met_file[:-4].replace('_', '.')))
            # Store full path to file.
            met_fls.append(os.path.join(iso_path, met_file))
        # Store values and paths to metallicity files for all the photometric
        # systems defined.
        met_vals_all.append(met_vals)
        met_files.append(met_fls)

    return met_vals_all, met_files


def findMetSubset(met_vals_all, met_files):
    """
    Find a subset of common z values among the photometric systems being
    processed. If a single photometric system is used, just return all the
    z values present in its isochrones folder.
    """

    if len(met_vals_all) == 1:
        # A single photometric system is being used. Sort before returning.
        met_vals_unq, met_files_unq = zip(
            *sorted(zip(met_vals_all[0], met_files[0])))
        return list(met_vals_unq), [list(met_files_unq)]

    else:
        # If more than one photometric system is being processed.

        # Subset of equal 'z' values across photometric systems.
        common_z = sorted(list(
            set.intersection(*[set(z) for z in met_vals_all])))

        if not common_z:
            sys.exit(
                "\nERROR: no common metallicity files were found across\n" +
                "the photometric systems.")

        met_files_c = []
        for i, mets in enumerate(met_vals_all):

            # Indexes of equal 'z' values.
            idxs = [mets.index(z) for z in common_z]

            met_files_ps = []
            for c_idx in idxs:
                met_files_ps.append(met_files[i][c_idx])

            # Sort from min to max
            met_files_c.append(sorted(met_files_ps))

        return common_z, met_files_c


def checkMetAge(par, p_rng, vals_all, round_f=5):
    """
    (Semi) Agnostic (z, log(age)) parameter check.

    Given the (z, log(age)) range in the 'params_input' file, check that this
    range, or fixed value, is contained by the available values in all the
    photometric systems being used.

    """

    if len(p_rng) > 1:
        p_rng[0] = min(vals_all) if p_rng[0] == 'min' else float(p_rng[0])
        p_rng[1] = max(vals_all) if p_rng[-1] == 'max' else float(p_rng[1])
        # Check that the values are different.
        if p_rng[0] == p_rng[1]:
            print("  WARNING: min & max {} values are equal ({})".format(
                par, p_rng[0]))
            p_rng = [p_rng[0]]
        elif p_rng[0] > p_rng[1]:
            sys.exit(
                ("\nERROR: minimum {} value '{:.5f}' is larger than\n"
                 "the maximum value '{:.5f}'").format(
                    par, p_rng[0], p_rng[1]))

    elif len(p_rng) == 1:
        # A fixed value is used.
        if p_rng[0] == 'min':
            p_rng = [min(vals_all)]
        elif p_rng[0] == 'max':
            p_rng = [max(vals_all)]

    if len(p_rng) > 1:
        if p_rng[0] < min(vals_all):
            sys.exit(
                ("\nERROR: {} defined '{}' is smaller than the\n"
                 "minimum value found '{:.5f}'").format(
                     par, p_rng[0], min(vals_all)))
        if max(vals_all) < p_rng[1]:
            sys.exit(
                ("\nERROR: {} defined '{}' is larger than the\n"
                 "maximum value found '{:.5f}'").format(
                     par, p_rng[1], max(vals_all)))

    elif len(p_rng) == 1:
        # A fixed value is used.
        p0 = p_rng[0]

        if par == 'log(age)':
            # Use same decimals as in 'logAges'
            flag = np.isclose(p0, vals_all, atol=10**-round_f).any()
        else:
            flag = p0 in vals_all

        if not flag:
            idx_r = np.searchsorted(vals_all, p0)
            if idx_r == len(vals_all) or idx_r == 0:
                idx_r = idx_r - 1 if idx_r == len(vals_all) else idx_r
                sys.exit(
                    ("\nERROR: could not find matching {} for the fixed\n" +
                     "value '{}'. The closest value found is: [{}]").format(
                        par, p0, vals_all[idx_r]))
            else:
                sys.exit(
                    ("\nERROR: could not find matching {} for the fixed\n" +
                     "value '{}'. The closest values found are: " +
                     "[{}, {}]").format(
                        par, p0, vals_all[idx_r - 1], vals_all[idx_r]))

    return p_rng


def filterMetvals(z_range, met_vals_all, met_files):
    """
    Filter all the available metallicity files in all the photometric systems
    to match either the single fixed z value defined, or its range.
    """

    # A single fixed z was defined. Return its associated file.
    if len(z_range) == 1:
        idx = met_vals_all.index(z_range[0])
        met_files_unq = []
        for phot_syst in met_files:
            met_files_unq.append([phot_syst[idx]])
        return z_range, met_files_unq

    else:
        # Indexes of range limits.
        idx_l, idx_h = np.searchsorted(met_vals_all, z_range)
        idx_l = idx_l - 1 if idx_l != 0 else idx_l

        # Store only values in range, including both extremes.
        met_vals_unq, met_files_unq = [], [[] for _ in met_files]
        for i, z in enumerate(met_vals_all):
            if idx_l <= i and i <= idx_h:
                met_vals_unq.append(z)
                for j, phot_syst in enumerate(met_files):
                    met_files_unq[j].append(phot_syst[i])

        return met_vals_unq, met_files_unq


def readAgevals(met_files, evol_track):
    """
    Read all the available ages in the filtered metallicity files in all the
    photometric systems defined.
    """

    # CMD PARSEC isochrone
    if evol_track[:3] == 'PAR':
        # Format of the line that contains the age in Gyrs (generated by the
        # ezPadova-2 code)
        age_format = isochs_format.age_format(evol_track)

        # For each photometric system
        age_vals_all = []
        for i, phot_syst in enumerate(met_files):

            # For each metallicity value
            common_a = []
            for met_file in phot_syst:
                with open(met_file, mode="r") as f_iso:
                    content = f_iso.read()
                    # Find all instances of the 'age_format' regular expression
                    # in this file.
                    ages0 = re.findall(age_format, content)
                    common_a.append(ages0)
            age_vals_all.append(common_a)
    else:
        # TODO in place for #275
        pass

    return age_vals_all


def findAgeSubset(age_vals_all, iso_paths):
    """
    age_vals_all.shape = (N photom systems, N met files, N ages)
    """
    ages_common = []
    for i, ages_phot_syst in enumerate(age_vals_all):
        common_a = list(set.intersection(*[set(z) for z in ages_phot_syst]))
        if not common_a:
            sys.exit(
                "\nERROR: no common age values were found across\n" +
                "metallicity files in the '{}' photometric system.".format(
                    iso_paths[i].split('/')[-1]))
        ages_common.append(common_a)

    # Common age values across photometric systems
    ages_common = list(set.intersection(*[set(z) for z in ages_common]))
    if not ages_common:
        sys.exit(
            "\nERROR: no common age values were found across\n" +
            "the photometric systems.")

    return ages_common


def logAges(ages_common, round_f=5):
    """
    Convert ages to log10() and sort from min to max.
    """
    # Convert all values to log10 rounding to 5 decimal places.
    age_vals_all = np.round(np.log10(list(map(float, ages_common))), round_f)
    sort_i = np.argsort(age_vals_all)

    # Sort min to max
    age_vals_all = age_vals_all[sort_i]
    ages_strs = np.array(ages_common)[sort_i].astype(str)

    return age_vals_all, ages_strs


def filterAgevals(a_range, age_vals_all, ages_strs):
    """
    Filter the available age values according to the given range.
    """
    if len(a_range) == 1:
        idx = np.argmin(abs(age_vals_all - a_range[0]))
        return a_range, [ages_strs[idx]]
    else:
        # Indexes of range limits.
        idx_l, idx_h = np.searchsorted(age_vals_all, a_range)
        idx_l = idx_l - 1 if idx_l != 0 else idx_l

        # Store only values in range, including both extremes.
        age_vals_unq, ages_strs_unq = [], []
        for i, a in enumerate(age_vals_all):
            if idx_l <= i and i <= idx_h:
                age_vals_unq.append(a)
                ages_strs_unq.append(ages_strs[i])

        return age_vals_unq, ages_strs_unq