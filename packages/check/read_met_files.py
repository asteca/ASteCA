
import sys
import traceback
import numpy as np
from packages.inp import met_ages_values
from packages.inp import isoch_params


def check_get(pd):
    """
    Check that all metallicity files needed are in place. To save time, we
    store the data and pass it.
    """
    # Only read files if best fit method is set to run. Else pass empty list.
    pd['fundam_params'], pd['theor_tracks'], pd['plot_isoch_data'] = [], [], []
    if pd['bf_flag']:

        # Read metallicity files' names, store proper ranges for all
        # parameters, and available metallicities and ages.
        try:
            # Get parameters values defined.
            params_values, met_f_filter, met_values, age_values, met_vals_all,\
                age_vals_all = met_ages_values.main(
                    pd['iso_paths'], pd['best_fit_algor'], pd['par_ranges'],
                    pd['za_steps'], pd['N_mass'])
        except Exception:
            print(traceback.format_exc())
            sys.exit("\nERROR: error storing metallicity files.")

        # Checks.
        ranges_files_check(
            params_values, met_values, age_values, met_vals_all, age_vals_all,
            **pd)

        # Store all the accepted values for the metallicity and age, and the
        # ranges of accepted values for the rest of the fundamental parameters.
        # The 'met_values' list contains duplicated sub-lists for each
        # photometric system defined. We store only one.
        pd['fundam_params'] = [met_values[0], age_values] + params_values[2:]

        # Store all isochrones in all the metallicity files.
        pd['theor_tracks'], pd['plot_isoch_data'] = isoch_params.main(
            met_f_filter, age_values, **pd)

    return pd


def ranges_files_check(
    params_values, met_values, age_values, met_vals_all, age_vals_all,
        par_ranges, cmd_systs, all_syst_filters, **kwargs):
    """
    Various checks.
    """
    # Check that ranges are properly defined.
    p_names = [
        ['metallicity', par_ranges[0]], ['age', par_ranges[1]],
        ['extinction', par_ranges[2]], ['distance', par_ranges[3]],
        ['mass', par_ranges[4]], ['binary', par_ranges[5]]]
    for i, p in enumerate(params_values):
        if not p.size:
            sys.exit("ERROR: No values exist for '{}' range defined:\n"
                     "min={}, max={}, step={}".format(p_names[i][0],
                                                      *p_names[i][1][1]))

    # Check that metallicity and age min, max & steps values are correct.
    # Match values in metallicity and age ranges with those available.
    z_range, a_range = params_values[:2]

    err_mssg = "ERROR: one or more metallicity files in the '{}' system\n" +\
               "could not be matched to the range given.\n\n" +\
               "The defined values are:\n\n" +\
               "{}\n\nThe read values are:\n\n" +\
               "{}\n\nThe closest values found are:\n\n" +\
               "{}\n\nThe missing elements are:\n\n{}"
    # Go through the values extracted from the metallicity files present
    # in each photometric system used.
    for i, met_vs in enumerate(met_values):
        if z_range.size > len(met_vs):
            # Find missing elements.
            missing = find_missing(z_range, met_vs)
            sys.exit(err_mssg.format(
                cmd_systs[all_syst_filters[i][0]][0], z_range,
                np.array(met_vals_all[i]), np.asarray(met_vs),
                np.asarray(missing)))
    err_mssg = "ERROR: one or more isochrones could not be matched\n" +\
               "to the age range given.\n\nThe defined values are:\n\n" +\
               "{}\n\nThe read values are:\n\n" +\
               "{}\n\nThe closest values found are:\n\n" +\
               "{}\n\nThe missing elements are:\n\n{}"
    if len(a_range) > len(age_values):
        # Find missing elements.
        missing = find_missing(a_range, age_values)
        sys.exit(err_mssg.format(
            a_range, age_vals_all, np.asarray(age_values),
            np.asarray(missing)))


def find_missing(arr_large, arr_small):
    '''
    Takes two arrays of floats, compares them and returns the missing
    elements in the smaller one.
    '''
    # Convert to strings before comparing.
    s1 = [str(_) for _ in np.round(arr_small, 5)]
    s2 = [str(_) for _ in np.round(arr_large, 5)]
    # Find missing elements. Convert to float and store.
    missing = [float(_) for _ in s2 if _ not in set(s1)]

    return missing
