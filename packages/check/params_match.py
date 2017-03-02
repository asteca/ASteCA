
import sys
import numpy as np
from os.path import isdir


def check(IMF_name, R_V, evol_track, bin_methods, bf_flag, bf_params,
          cmd_evol_tracks, iso_paths, imf_funcs, par_ranges, bin_mr,
          ga_params, **kwargs):
    """
    Check all parameters related to the search for the best synthetic cluster
    match.
    """

    # Unpack.
    best_fit_algor, lkl_method, bin_method, N_b = bf_params
    # If best fit method is set to run.
    if bf_flag:

        # Check best fit method selected.
        if best_fit_algor not in ('brute', 'genet'):
            sys.exit("ERROR: the selected best fit method '{}' does not match"
                     " a valid input.".format(best_fit_algor))

        if best_fit_algor == 'genet':
            # Check GA input params.
            n_pop, n_gen, fdif, p_cross, cr_sel, p_mut, n_el, n_ei, n_es = \
                ga_params
            # First set of params.
            oper_dict0 = {'n_pop': n_pop, 'n_gen': n_gen, 'n_el': n_el,
                          'n_ei': n_ei, 'n_es': n_es}
            for oper in oper_dict0:
                if oper_dict0[oper] < 1:
                    sys.exit("ERROR: number must be greater than zero in\n"
                             "'{}' GA parameter; {} is set.".format(
                                 oper, oper_dict0[oper]))
            # Second set of params.
            oper_dict1 = {'fdif': fdif, 'p_cross': p_cross, 'p_mut': p_mut}
            for oper in oper_dict1:
                if oper_dict1[oper] < 0. or oper_dict1[oper] > 1.:
                    sys.exit("ERROR: GA '{}' parameter is out of the valid\n"
                             "[0., 1.] range; {} is set.".format(
                                 oper, oper_dict1[oper]))
            # Handle separately.
            if cr_sel not in ('1P', '2P'):
                sys.exit("ERROR: GA 'cr_sel' operator is not a valid choice;\n"
                         "'{}' is set.".format(cr_sel))
            # Number of solutions to pass to the nest generation (elitism)
            if n_el >= n_pop:
                sys.exit("ERROR: GA 'n_el' must be smaller than 'n_pop';\n"
                         "'{}' and '{}' are set respectively.".format(
                             n_el, n_pop))

        # Check likelihood method selected.
        if lkl_method not in ('tolstoy', 'dolphin'):
            sys.exit("ERROR: the selected likelihood method '{}' does not"
                     " match a valid input.".format(lkl_method))

        # Check binning method selected.
        if lkl_method == 'dolphin' and bin_method not in bin_methods:
            sys.exit("ERROR: the selected binning method '{}' for the 'Best"
                     "\nfit' function does not match a valid input."
                     .format(bin_method))

        # Check if /isochrones folder exists.
        for iso_path in iso_paths:
            if not isdir(iso_path):
                sys.exit(
                    "ERROR: 'Best synthetic cluster fit' function is set to"
                    " run but the folder:\n\n {}\n\ndoes not exists."
                    .format(iso_path))

        # Check selected isochrones set.
        if evol_track not in cmd_evol_tracks.keys():
            sys.exit("ERROR: the selected isochrones set ('{}') does\n"
                     "not match a valid input.".format(evol_track))

        # Check IMF defined.
        if IMF_name not in imf_funcs:
            sys.exit("ERROR: Name of IMF ({}) is incorrect.".format(IMF_name))

        # Check R_V defined.
        if R_V <= 0.:
            sys.exit("ERROR: Ratio of total to selective absorption\n"
                     "R_V ({}) must be positive defined.".format(R_V))

        # Check that no parameter range is empty.
        m_rs, a_rs, e_rs, d_rs, mass_rs, bin_rs = par_ranges
        p_names = [['metallicity', m_rs], ['age', a_rs], ['extinction', e_rs],
                   ['distance', d_rs], ['mass', mass_rs], ['binary', bin_rs]]
        if min(mass_rs[1]) == 0.:
            print("  WARNING: minimum total mass is zero in params_input"
                  " file.\n  Changed minimum mass to 10.\n")
            mass_rs[1][mass_rs[1].index(min(mass_rs[1]))] = 10.

        for i, p in enumerate(par_ranges):
            # Catch empty list.
            if not p:
                sys.exit("ERROR: Range defined for '{}' parameter is"
                         " empty.".format(p_names[i][0]))
            # Catch *almost* empty list since inp/input_params perhaps added
            # an identifier 'l' or 'r'. This prevents ranges given as
            # empty lists (ie: [] or () or {}) from passing as valid ranges.
            elif not p[1]:
                sys.exit("ERROR: Range defined for '{}' parameter is"
                         " empty.".format(p_names[i][0]))

        # Check binarity parameters.
        # See if it is a list of values or a range.
        if par_ranges[-1][0] == 'r':
            # Range: min, max, step. Store min and max in array.
            if len(par_ranges[-1][-1]) > 1:
                # More than one value
                bin_fr = np.array([par_ranges[-1][-1][0],
                                  par_ranges[-1][-1][1]])
            else:
                # Single value stored.
                bin_fr = np.array([par_ranges[-1][-1][0]])
        else:
            # List: store all values in array.
            bin_fr = np.array(par_ranges[-1][-1])
        # Check all values in array.
        for bin_fr_val in bin_fr:
            if bin_fr_val > 1.:
                sys.exit("ERROR: Binarity fraction value '{}' is out of\n"
                         "boundaries. Please select a value in the range "
                         "[0., 1.]".format(bin_fr_val))
        if not 0. <= bin_mr <= 1.:
            sys.exit("ERROR: Binary mass ratio set ('{}') is out of\n"
                     "boundaries. Please select a value in the range [0., 1.]".
                     format(bin_mr))
