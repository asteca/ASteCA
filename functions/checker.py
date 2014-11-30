# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 09:00:00 2014

@author: gabriel
"""
import pip
from os.path import join, isfile, isdir
import sys
import traceback
from subprocess import Popen, PIPE
from get_in_params import get_in_params as gip
from get_names_paths import names_paths as n_p
import get_isoch_params as isochp


def check(mypath, cl_files):
    '''
    Checks that the necessary files are in place, that the parameters stored
    in the input file are valid and that the ranges given for the cluster
    parameters are consistent with the isochrones available before moving
    on with the code.
    '''

    # Check packages installed.
    inst_packgs = pip.get_installed_distributions()
    inst_packgs_lst = ["%s" % (i.key) for i in inst_packgs]
    for pckg in ['numpy', 'matplotlib', 'scipy', 'astroml', 'scikit-learn',
        'astropy']:
        if pckg not in inst_packgs_lst:
            print "ERROR: '{}' package is not installed.".format(pckg)
            sys.exit("Install with: pip install {}".format(pckg))
    import numpy as np

    # Check if params_input.dat file exists.
    if not isfile(join(mypath, 'params_input.dat')):
        sys.exit('ERROR: params_input.dat file does not exist.')

    # Check if params_input file is properly formatted.
    try:
        # Read input parameters from params_input.dat file.
        mode, done_dir, gd_params, gh_params, cr_params, kp_flag,\
        im_flag, er_params, fr_number, pv_params, da_params, ps_params,\
        bf_params, sc_params, ga_params, rm_params, pl_params,\
        flag_move_file, axes_params = gip(mypath)
    except Exception:
        # Halt code.
        print traceback.format_exc()
        sys.exit('ERROR: params_input.dat is badly formatted.')

    # Check mode.
    if mode not in {'auto', 'semi', 'manual'}:
        sys.exit("ERROR: 'mode' value selected ({}) is not valid.".format(mode))

    if mode == 'semi':
        # Check if semi_input.dat file exists.
        semi_file = 'semi_input.dat'
        if not isfile(join(mypath, semi_file)):
            # File semi_input.dat does not exist.
            sys.exit("ERROR: 'semi' mode is set but semi_input.dat file does"
                "not exist.")

    # Check px/deg.
    if gd_params[-1] not in {'px', 'deg'}:
        sys.exit("ERROR: the coordinates given in the input file ({})"
                " are incorrect.".format(gd_params[-1]))

    # Selected CMD.
    if ps_params[1] not in {1, 2, 3, 4, 5, 6, 7}:
        sys.exit("ERROR: CMD selected ({}) is not valid.".format(ps_params[1]))

    # Output figure.
    if pl_params[0] is True:
        if pl_params[1] not in {'png', 'pdf', 'PNG', 'PDF'}:
            sys.exit("ERROR: figure output format selected ({}) is"
            "not valid.".format(pl_params[1]))

    # 2D positional histogram.
    if gh_params[0] not in {'auto', 'manual'}:
        sys.exit("ERROR: mode selected ({}) for 2D histogram"
        "is not valid.".format(gh_params[0]))

    # Radius finding function.
    if cr_params[0] not in {'auto', 'manual'}:
        sys.exit("ERROR: mode selected ({}) for radius finding"
        "function is not valid.".format(cr_params[0]))
    if cr_params[0] is 'manual' and cr_params[1] < 4:
        print ("  WARNING: number of points in manual radius\n"
        "  mode is {} < 4. A value of 4 will be used.".format(
        cr_params[1]))

    # Errors function.
    if er_params[0] not in {'emax', 'lowexp', 'eyefit'}:
        sys.exit("ERROR: mode selected ({}) for error rejecting"
        "function is not valid.".format(er_params[0]))
    if er_params[0] == 'emax' and len(er_params[1:]) < 1:
        sys.exit("ERROR: missing parameters for error rejecting function")
    if er_params[0] == 'eyefit' and len(er_params[1:]) < 3:
        sys.exit("ERROR: missing parameters for error rejecting function")
    if er_params[0] == 'lowexp' and len(er_params[1:]) < 4:
        sys.exit("ERROR: missing parameters for error rejecting function")

    # Check KDE p-value custer probability function.
    R_in_place = False
    if pv_params[0] not in {'auto', 'manual', 'skip'}:
        sys.exit("ERROR: Wrong name ({}) for KDE p-value function "
            "'mode'.".format(pv_params[0]))

    elif pv_params[0] in {'auto', 'manual'}:

        rpy2_inst, R_inst = True, True
        # Check if rpy2 package is installed.
        if 'rpy2' not in inst_packgs_lst:
            rpy2_inst = False
        # Now check for R.
        proc = Popen(["which", "R"], stdout=PIPE, stderr=PIPE)
        exit_code = proc.wait()
        if exit_code != 0:
            R_inst = False

        # If both R and rpy2 packages are installed.
        if R_inst and rpy2_inst:
            # R and rpy2 package are installed, function is good to go.
            R_in_place = True
        else:
            if R_inst and not rpy2_inst:
                R_pack = 'rpy2 is'
            if rpy2_inst and not R_inst:
                R_pack = 'R is'
            if not R_inst and not rpy2_inst:
                R_pack = 'R and rpy2 are'
            # Something is not installed and function was told to run.
            print ("  WARNING: {} not installed and the\n"
            "  'KDE p-value test' was set to run. The\n"
            "  function will be skipped.\n".format(R_pack))
    #else:
    # Something is not installed, but function was told not to run so
    # it will be skipped anyway.

    # Check decont algorithm params.
    mode_da = da_params[0]
    # Check if 'mode' was correctly set.
    if mode_da not in ['auto', 'manual', 'read', 'skip']:
        sys.exit("ERROR: Wrong name ({}) for decontamination algorithm "
            "'mode'.".format(mode_da))

    if mode_da == 'read':
        # Check if file exists.
        for cl_file in cl_files:

            # Get memb file names.
            memb_file = n_p(mypath, cl_file)[2]
            if not isfile(memb_file):
                # File does not exist.
                sys.exit("ERROR: 'read' mode was set for decontamination "
                "algorithm but the file:\n\n {}\n\ndoes not "
                "exist.".format(memb_file))

    ip_list = []
    # If best fit method is set to run.
    if bf_params[0]:

        # Unpack.
        iso_path, cmd_select, iso_select, par_ranges = ps_params
        m_rs, a_rs, e_rs, d_rs, mass_rs, bin_rs = par_ranges

        # Check that CMD is correctly set.
        if cmd_select not in {1, 2, 3, 4, 5, 6, 7}:
            sys.exit("ERROR: the stored CMD value ({}) does not match a valid"
                " selection.".format(cmd_select))

        # Check selected isochrones set.
        if iso_select not in {'MAR', 'PAR'}:
            sys.exit("ERROR: the selected isochrones set ({}) does not match a"
                " valid input.".format(iso_select))

        # Check if /isochrones folder exists.
        if not isdir(iso_path):
            sys.exit("ERROR: 'Best synthetic cluster fit' function is set to"
                " run but the folder:\n\n {}\n\ndoes not exists.".format(
                    iso_path))

        # Read names of all metallicity files stored in isochrones path given.
        # Also read full paths to metallicity files.
        met_vals_all, met_files = isochp.get_metals(iso_path)
        # Read Girardi metallicity files format.
        isoch_format = isochp.i_format(iso_select, cmd_select)
        # Read all ages in the first metallicity file: met_files[0]
        # *WE ASUME ALL METALLICITY FILES HAVE THE SAME NUMBER OF AGE VALUES*
        age_vals_all = isochp.get_ages(met_files[0], isoch_format[1])
        # Get parameters ranges stored in params_input.dat file.
        param_ranges, param_rs = isochp.get_ranges(par_ranges)
        # Check that ranges are properly defined.
        p_names = [['metallicity', m_rs], ['age', a_rs], ['extinction', e_rs],
            ['distance', d_rs], ['mass', mass_rs], ['binary', bin_rs]]
        for i, p in enumerate(param_ranges):
            if not p.size:
                sys.exit("ERROR: No values exist for {} range defined:\n\n"
                "min={}, max={}, step={}".format(p_names[i][0], *p_names[i][1]))

        # Check that metallicity and age min, max & steps values are correct.
        # Match values in metallicity and age ranges with those available.
        z_range, a_range = param_ranges[:2]
        met_f_filter, met_values, age_values = isochp.match_ranges(met_vals_all,
            met_files, age_vals_all, z_range, a_range)

        if len(z_range) > len(met_values):
            sys.exit("ERROR: one or more metallicity files could not be\n"
            "matched to the range given. The range defined was:\n\n"
            "{}\n\nand the closest available values are:\n\n"
            "{}".format(z_range, np.asarray(met_values)))
        if len(a_range) > len(age_values):
            sys.exit("ERROR: one or more isochrones could not be matched\n"
            "to the age range given. The range defined was:\n\n"
            "{}\n\nand the closest available values are:\n\n"
            "{}".format(a_range, np.asarray(age_values)))

        # Read metallicity files.
        try:
            # Store all isochrones in all the metallicity files in isoch_list.
            ip_list = isochp.ip(ps_params, bf_params[0])
        except:
            print traceback.format_exc()
            sys.exit("ERROR: unknown error reading metallicity files.")

    # Check IMF defined.
    imfs_dict = {'kroupa_1993', 'chabrier_2001', 'kroupa_2002'}
    if sc_params[0] not in imfs_dict:
        sys.exit("ERROR: Name of IMF ({}) is incorrect.".format(sc_params[0]))

    # Check that at least one photometric cluster file exists.
    if not cl_files:
        sys.exit("No photometric data files found in '/input' folder. Halting.")

    print 'Full check done.\n'
    return ip_list, R_in_place