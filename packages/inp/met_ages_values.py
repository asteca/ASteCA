
import numpy as np
import sys
import os
import re
from os.path import join
from . import isochs_format


def match_ranges(met_vals_all, met_files, age_vals_all, z_range, a_range):
    '''
    Matches available metallicity and ages values with those stored in the
    ranges given to these two parameters.
    '''

    # Match metallicity values in ranges with values available.
    met_f_filter, met_values = [], []
    for i, met_vals in enumerate(met_vals_all):
        met_f_f, met_vs = [], []
        for j, met in enumerate(met_vals):
            # Store metallicity file only if it's inside the given range.
            if np.isclose(z_range, met, atol=0.00001).any():
                met_f_f.append(met_files[i][j])
                met_vs.append(met)
        met_f_filter.append(met_f_f)
        met_values.append(met_vs)

    # Match age values in ranges with values available.
    age_values = []
    for age in age_vals_all:
        # If age value falls inside the given range, store the value.
        # TODO hardcoded 0.001 tolerance, could break something sometime?
        if np.isclose(a_range, age, atol=0.001).any():
            age_values.append(age)

    if not age_values:
        txt = "No age value read:\n\n{}\n\ncould be matched to the ages" +\
            " given as input:\n\n{}"
        sys.exit(txt.format(age_vals_all, a_range))

    return met_f_filter, met_values, age_values


def get_ranges(par_ranges):
    '''
    Calculate parameter ranges to be used by the selected best fit method.
    '''
    param_ranges = []
    for param in par_ranges:
        # List of values.
        if param[0] == 'l':
            # Store ordered values.
            param[1].sort()
            param_ranges.append(np.asarray(param[1]))

        # Range of values.
        else:
            # If only one value is defined.
            if len(param[1]) == 1:
                param_ranges.append(np.asarray([param[1][0]]))
            # If min == max store single value in array.
            elif param[1][0] == param[1][1]:
                param_ranges.append(np.asarray([param[1][0]]))
            else:
                # Store range values in array.
                p_rang = np.arange(*param[1])
                # Skip if array is empty. Checker will catch this.
                if p_rang.size:
                    # Add max value if not present. Check this way to avoid
                    # floating point errors.
                    if abs(p_rang[-1] - param[1][1]) > 0.0000001:
                        p_rang = np.append(p_rang, param[1][1])
                # Store full range for this parameter.
                param_ranges.append(p_rang)

    return param_ranges


def CMDAges(met_file):
    '''
    Read all available ages in a CMD service metallicity file.
    '''

    age_format = isochs_format.cmd_age_format()

    # Open the metallicity file.
    with open(met_file, mode="r") as f_iso:
        regex = age_format  # Define regular expression.
        ages0 = re.findall(regex, f_iso.read())  # Find all instances.
        ages1 = np.asarray(list(map(float, ages0)))
        ages2 = np.log10(ages1)
        # TODO hardcoded, could have unintended consequences
        isoch_a = np.around(ages2, 3)  # Round to 3 decimals.

    return isoch_a


def get_metals(iso_paths):
    '''
    Read names of all metallicity files stored in each isochrones' path, and
    store them along with the z values they represent.
    '''

    met_vals_all, met_files = [], []
    # For each photometric system used.
    for iso_path in iso_paths:
        # List all metallicity files in folder.
        metal_files = sorted(os.listdir(iso_path))
        # Iterate in order through all the metallicity files stored for the
        # selected set of isochrones.
        met_vals, met_fls = [], []
        for met_file in metal_files:
            # *THE NAME OF THE FILE IS IMPORTANT*
            # Extract metallicity value from the name of the file.
            # Replace underscores in file names with decimal points.
            met_vals.append(float(met_file[:-4].replace('_', '.')))
            # Store full path to file.
            met_fls.append(join(iso_path, met_file))
        # Store values and paths to metallicity files for all the photometric
        # systems defined.
        met_vals_all.append(met_vals)
        met_files.append(met_fls)

    return met_vals_all, met_files


def main(iso_paths, par_ranges):
    '''
    Obtain the correct metallicities and ages used by the code.
    '''

    # Read names of all metallicity files stored in the isochrones paths given.
    # This stores all the metallicity values available.
    # Also read full paths to metallicity files.
    met_vals_all, met_files = get_metals(iso_paths)

    # Read all ages from the *first* metallicity file defined, assuming all
    # files contain the same number of ages.
    # This check (that all metallicity files contain the same number of ages)
    # happens in 'read_isochs'.
    age_vals_all = CMDAges(met_files[0][0])

    # Get parameters ranges stored in params_input.dat file.
    param_ranges = get_ranges(par_ranges)

    # Match values in metallicity and age ranges given by the user, with
    # those available in the theoretical isochrones.
    z_range, a_range = param_ranges[:2]
    met_f_filter, met_values, age_values = match_ranges(
        met_vals_all, met_files, age_vals_all, z_range, a_range)

    return param_ranges, met_f_filter, met_values, age_values, met_vals_all,\
        age_vals_all
