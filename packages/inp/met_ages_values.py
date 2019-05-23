
import numpy as np
import sys
import os
import re
from os.path import join
from . import isochs_format


def main(iso_paths, best_fit_algor, par_ranges, za_steps, ga_steps, m_step):
    '''
    Obtain the correct metallicities and ages used by the code.
    '''

    # Read names of all metallicity files stored in the isochrones paths given.
    # This stores all the metallicity values available.
    # Also read full paths to metallicity files.
    met_vals_all, met_files = get_metals(iso_paths)

    # Read all ages from the *first* metallicity file defined in the *first*
    # photometric system stored, assuming all files contain the same number of
    # ages. This check (that all metallicity files contain the same number of
    # ages) happens in 'read_isochs'.
    age_vals_all = CMDAges(met_files[0][0])

    # Get parameters ranges stored in params_input.dat file.
    params_values = getParamVals(
        best_fit_algor, par_ranges, za_steps, ga_steps, m_step)

    # Match values in metallicity and age ranges given by the user, with
    # those available in the theoretical isochrones.
    z_range, a_range = params_values[:2]
    met_f_filter, met_values, age_values = match_ranges(
        met_vals_all, met_files, age_vals_all, z_range, a_range)

    return params_values, met_f_filter, met_values, age_values, met_vals_all,\
        age_vals_all


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
        # Remove hidden files possibly present in folder.
        metal_files = [_ for _ in metal_files if not _.startswith('.')]

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
        # TODO hardcoded rounding. Must be modified in accordance with
        # 'read_isochs.readCMDFile()'.
        isoch_a = np.around(np.log10(ages1), 4)

    return isoch_a


def getParamVals(best_fit_algor, par_ranges, za_steps, ga_steps, m_step):
    '''
    Obtain parameter ranges to be used by the selected best fit method.
    '''

    # Define steps for the parameter ranges according to the best fit algorithm
    # selected.
    if best_fit_algor == 'boot+GA':
        steps = za_steps + ga_steps[:2] + [m_step, ga_steps[2]]
    else:
        steps = za_steps + [None, None] + [m_step, None]

    params_values = []
    for i, param in enumerate(par_ranges):
        # If only one value is defined.
        if len(param) == 1:
            params_values.append(np.asarray([param[0]]))
        # If min == max store single value in array.
        elif param[0] == param[1]:
            params_values.append(np.asarray([param[0]]))
        else:
            # Store range values in array.
            if steps[i] is None:
                p_rang = np.asarray(param)
            else:
                p_rang = np.arange(param[0], param[1], steps[i])
            # Skip if array is empty. Checker will catch this.
            if p_rang.size:
                # Add max value if not present. Check this way to avoid
                # floating point errors.
                if abs(p_rang[-1] - param[1]) > 0.0000001:
                    p_rang = np.append(p_rang, param[1])
            # Store full range for this parameter.
            params_values.append(p_rang)

    return params_values


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
        # TODO hardcoded tolerance, could break something sometime?
        if np.isclose(a_range, age, atol=0.001).any():
            age_values.append(age)

    if not age_values:
        txt = "No age value read:\n\n{}\n\ncould be matched to the ages" +\
            " given as input:\n\n{}"
        sys.exit(txt.format(age_vals_all, a_range))

    return met_f_filter, met_values, age_values
