
import numpy as np
import os
import re
from os.path import join
from ..inp import input_params as g
import isochs_format


def match_ranges(met_vals_all, met_files, age_vals_all, z_range, a_range):
    '''
    Matches available metallicity and ages values with those stored in the
    ranges given to these two parameters.
    '''

    # Match metallicity values in ranges with values available.
    met_f_filter, met_values = [], []
    for i, met in enumerate(met_vals_all):
        # Store metallicity file only if it's inside the given range.
        if np.isclose(z_range, met, atol=0.00001).any():
            met_f_filter.append(met_files[i])
            met_values.append(met)

    # Match age values in ranges with values available.
    age_values = []
    for age in age_vals_all:
        # If age value falls inside the given range, store the value.
        if np.isclose(a_range, age, atol=0.01).any():
            age_values.append(round(age, 2))

    return met_f_filter, met_values, age_values


def get_ranges():
    '''
    Calculate parameter ranges to be used by the selected best fit method.
    '''

    par_ranges = g.ps_params[-1]

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


def get_ages(met_file):
    '''
    Read all available ages in metallicity file.
    '''

    age_format = isochs_format.girardi_age_format()

    # Open the metallicity file.
    with open(met_file, mode="r") as f_iso:
        regex = age_format  # Define regular expression.
        ages0 = re.findall(regex, f_iso.read())  # Find all instances.
        ages1 = np.asarray(map(float, ages0))  # Map to floats.
        ages2 = np.log10(ages1)  # Take log10
        isoch_a = np.around(ages2, 2)  # Round to 2 decimals.

    return isoch_a


def get_metals(iso_path):
    '''
    Read names of all metallicity files stored in isochrones path given and
    store them along with the z values they represent.
    '''

    metal_files = sorted(os.listdir(iso_path))
    # Iterate in order through all the metallicity files stored for the
    # selected set of isochrones.
    met_vals_all, met_files = [], []
    for met_file in metal_files:
        # *THE NAME OF THE FILE IS IMPORTANT*
        # Extract metallicity value from the name of the file.
        # Replace underscores in file names with decimal points.
        met_vals_all.append(float(met_file[:-4].replace('_', '.')))
        # Store full path to file.
        met_files.append(join(iso_path, met_file))

    return met_vals_all, met_files


def main(iso_path):
    '''
    Run once to obtain the correct metallicities and ages to be used
    by the code.
    '''

    # Read names of all metallicity files stored in isochrones path given.
    # I.e.: store all metallicity values available.
    # Also read full paths to metallicity files.
    met_vals_all, metal_files = get_metals(iso_path)

    # Read all ages from the *first* metallicity file defined.
    #
    # *WE ASUME ALL METALLICITY FILES HAVE THE SAME NUMBER OF AGE VALUES*
    # (that's why we use the first metallicity file stored to obtain all
    # the age values)
    #
    # I.e: all metallicity files should contain the same amount and values for
    # the ages, otherwise something *will* fail down the line.
    age_vals_all = get_ages(metal_files[0])

    # Get parameters ranges stored in params_input.dat file.
    param_ranges = get_ranges()

    # Match values in metallicity and age ranges with those available.
    z_range, a_range = param_ranges[:2]
    met_f_filter, met_values, age_values = match_ranges(
        met_vals_all, metal_files, age_vals_all, z_range, a_range)

    return param_ranges, met_f_filter, met_values, age_values
