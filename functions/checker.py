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
from get_isoch_params import ip


def check(mypath, cl_files):
    '''
    Checks that the necessary files are in place and the parameters
    given are consistent with the isochrones available before moving
    on with the code.
    '''

    # Check packages installed.
    inst_packgs = pip.get_installed_distributions()
    inst_packgs_lst = ["%s" % (i.key) for i in inst_packgs]
    for pckg in ['numpy', 'matplotlib', 'scipy']:
        if pckg not in inst_packgs_lst:
            print "FATAL: '{}' package is not installed.".format(pckg)
            sys.exit("Install with: pip install {}".format(pckg))

    # Check if params_input.dat file exists.
    if not isfile(join(mypath, 'params_input.dat')):
        sys.exit('FATAL: params_input.dat file does not exist.')

    # Check if params_input file is properly formatted.
    try:
        # Read input parameters from params_input.dat file.
        mode, done_dir, gd_params, gh_params, gc_params, cr_params, kp_flag,\
        im_flag, er_params, fr_number, pv_params, da_params, ps_params,\
        bf_params, sc_params, ga_params, rm_params, pl_params, flag_move_file,\
        axes_params = gip(mypath)
    except Exception:
        # Halt code.
        print traceback.format_exc()
        sys.exit('FATAL: params_input.dat is badly formatted.')

    # Check mode.
    if mode not in {'auto', 'semi', 'manual'}:
        sys.exit("FATAL: 'mode' value is incorrect.")

    if mode == 'semi':
        # Check if semi_input.dat file exists.
        if not isfile(join(mypath, 'semi_input.dat')):
            # File semi_input.dat does not exist.
            sys.exit("FATAL: 'semi' mode is set but semi_input.dat file does"
                " not exist.")

    # Check KDE p-value custer probability function.
    # Check if rpy2 package and R are installed, else skip get_p_value function.
    rpy2_inst, R_inst = True, True
    if 'rpy2' not in inst_packgs_lst:
        rpy2_inst = False
    # Now for R.
    proc = Popen(["which", "R"], stdout=PIPE, stderr=PIPE)
    exit_code = proc.wait()
    if exit_code != 0:
        R_inst = False

    R_in_place = False
    # If both R and rpy2 packages are installed.
    if R_inst and rpy2_inst:
        # R and rpy2 package are installed, function is good to go.
        R_in_place = True
    else:
        if pv_params[0]:
            if R_inst and not rpy2_inst:
                R_pack = 'rpy2 is'
            if rpy2_inst and not R_inst:
                R_pack = 'R is'
            if not R_inst and not rpy2_inst:
                R_pack = 'R and rpy2 are'
            # Something is not installed and function was told to run.
            print ("WARNING: {} not installed and the\n"
            "'KDE p-value test' was set to run. The function\n"
            "will be skipped.".format(R_pack))
        #else:
        # Something is not installed, but function was told not to run so
        # it will be skipped anyway.

    # Check decont algorithm params.
    mode_da = da_params[0]
    # Check if 'mode' was correctly set.
    if mode_da not in ['auto', 'manual', 'read', 'skip']:
        sys.exit("  FATAL: Wrong name ({}) for decontamination algorithm "
            "'mode'".format(mode_da))

    if mode_da == 'read':
        # Check if file exists.
        for cl_file in cl_files:

            # Get memb file names.
            memb_file = n_p(mypath, cl_file)[2]
            if not isfile(memb_file):
                # File does not exist.
                sys.exit("FATAL: 'read' mode was set for decontamination "
                "algorithm but the file:\n {}\ndoes not "
                "exist.".format(memb_file))

    # If best fit method is set to run.
    if bf_params[0]:

        # Unpack.
        iso_path, cmd_select, iso_select, m_rs, a_rs = ps_params[:-2]

        # Check that CMD is correctly set.
        if cmd_select not in {1, 2, 3, 4, 5, 6, 7}:
            sys.exit("FATAL: the stored CMD value ({}) does not match a valid"
                " selection.".format(cmd_select))

        # Check selected isochrones set.
        if iso_select not in {'MAR', 'PAR'}:
            sys.exit("FATAL: the selected isochrones set ({}) does not match a"
                " valid input.".format(iso_select))

        # Check if /isochrones folder exists.
        if not isdir(iso_path):
            sys.exit("FATAL: 'Best synthetic cluster fit' function is set to"
                " run but the folder:\n {}\ndoes not exists.".format(iso_path))

        # Check that metallicity and age min, max & steps values are correct.
        try:
            # Store all isochrones in all the metallicity files in isoch_list.
            # Store metallicity values and isochrones ages between the allowed
            # ranges in isoch_ma; extinction and distance modulus values in
            # isoch_ed.
            # isoch_list, isoch_ma, isoch_ed = ip_list
            # Only read files if best fit process is set to run.
            # bf_flag = bf_params[0]
            ip_list = ip(ps_params, bf_params[0])
        except:
            print traceback.format_exc()
            sys.exit("FATAL: error reading metallicity files.")

    # Check that at least one photometric cluster file exists.
    if not cl_files:
        sys.exit("No photometric data files found in '/input' folder. Halting.")

    print 'Check done.\n'
    return ip_list, R_in_place