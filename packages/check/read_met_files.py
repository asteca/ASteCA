
import sys
import traceback
import numpy as np
from packages.inp import met_ages_values
from packages.inp import isoch_params


def find_missing(arr_large, arr_small):
    '''
    Takes two arrays of floats, compares them and returns the missing
    elements in the smaller one.
    '''
    # Convert to strings before comparing.
    s1 = [str(_) for _ in arr_small]
    s2 = [str(_) for _ in arr_large]
    # Find missing elements. Convert to float and store.
    missing = [float(_) for _ in s2 if _ not in set(s1)]

    return missing


def ranges_files_check(iso_paths, par_ranges, **kwargs):
    """
    Obtain allowed metallicities and ages. Use the first photometric
    system defined.
    *WE ASUME ALL PHOTOMETRIC SYSTEMS CONTAIN THE SAME NUMBER OF
    METALLICITY FILES*
    """
    try:
        # Get parameters values defined.
        param_ranges, met_f_filter, met_values, age_values = \
            met_ages_values.main(iso_paths, par_ranges)
    except:
        print traceback.format_exc()
        sys.exit("\nERROR: error storing metallicity files.")

    # Check that ranges are properly defined.
    m_rs, a_rs, e_rs, d_rs, mass_rs, bin_rs = par_ranges
    p_names = [['metallicity', m_rs], ['age', a_rs], ['extinction', e_rs],
               ['distance', d_rs], ['mass', mass_rs], ['binary', bin_rs]]
    for i, p in enumerate(param_ranges):
        if not p.size:
            sys.exit("ERROR: No values exist for '{}' range defined:\n"
                     "min={}, max={}, step={}".format(p_names[i][0],
                                                      *p_names[i][1][1]))

    # Check that metallicity and age min, max & steps values are correct.
    # Match values in metallicity and age ranges with those available.
    z_range, a_range = param_ranges[:2]

    err_mssg = "ERROR: one or more metallicity files\n" +\
               "could not be matched to the range given.\n\n" +\
               "The defined values are:\n\n" +\
               "{}\n\nand the values read are:\n\n" +\
               "{}\n\nThe missing elements are:\n\n{}"
    # Go through the values extracted from the metallicity files present
    # in each photometric system used.
    for met_vs in met_values:
        if len(z_range) > len(met_vs):
            # Find missing elements.
            missing = find_missing(z_range, met_vs)
            sys.exit(err_mssg.format(z_range, np.asarray(met_vs),
                     np.asarray(missing)))
    err_mssg = "ERROR: one or more isochrones could not be matched\n" +\
               "to the age range given.\n\nThe defined values are:\n\n" +\
               "{}\n\nand the values read are:\n\n" +\
               "{}\n\nThe missing elements are:\n\n{}"
    if len(a_range) > len(age_values):
        # Find missing elements.
        missing = find_missing(a_range, age_values)
        sys.exit(err_mssg.format(a_range, np.asarray(age_values),
                 np.asarray(missing)))

    return param_ranges, met_f_filter, met_values, age_values


def check_get(pd):
    """
    Check that all metallicity files needed are in place. To save time, we
    store the data and pass it.
    """
    # Only read files if best fit method is set to run. Else pass empty list.
    pd['fundam_params'], pd['theor_tracks'] = [], []
    if pd['bf_flag']:
        # Read metallicity files' names, store proper ranges for all
        # parameters, and available metallicities and ages.
        param_ranges, met_f_filter, met_values, age_values =\
            ranges_files_check(**pd)
        # Store all the accepted values for the metallicity and age, and the
        # ranges of accepted values for the rest of the fundamental parameters.
        # The 'metal_values' list contains duplicated sub-lists for each
        # photometric system defined. We store only one.
        pd['fundam_params'] = [met_values[0], age_values] + param_ranges[2:]

        # Store all isochrones in all the metallicity files.
        pd['theor_tracks'] = isoch_params.main(met_f_filter, age_values, **pd)

    return pd
