# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 14:03:44 2014

@author: gabriel
"""

from os.path import join, isfile
import sys
#import re


def get_in_params(mypath):
    '''
    This function reads the input data parameters stored in the
    'ocaat_input.dat' file and returns them packaged for each function to use.
    '''

    # Store input and output dirs and path to input data file.
    input_dir = join(mypath, 'input/')
    output_dir = join(mypath, 'output/')
    data_file = join(input_dir, 'ocaat_input.dat')

    # Check if ocaat_input.dat file exists.
    if isfile(data_file):
        pass
    else:
        # Halt code.
        sys.exit('FATAL: ocaat_input.dat file does not exist. Halting code.')

    # Read data from file.
    true_lst = ('True', 'true')

    with open(data_file, "r") as f_dat:

        # Iterate through each line in the file.
        for line in f_dat:

            if not line.startswith("#") and line.strip() != '':
                reader = line.split()

                # Read folder paths where clusters are stored.
                if reader[0] == 'MO':
                    mode = str(reader[1])

                elif reader[0] == 'MF':
                    flag_move_file = True if reader[1] in true_lst else False
                    done_dir = str(reader[2])

                elif reader[0] == 'MP':
                    flag_make_plot = True if reader[1] in true_lst else False
                    plot_frmt = str(reader[2])
                    plot_dpi = int(reader[3])

                elif reader[0] == 'PD':
                    gd_params = map(int, reader[1:])

                elif reader[0] == 'CMD':
                    cmd_select = int(reader[1])

                elif reader[0] == 'PS':
                    iso_select = str(reader[1])

                elif reader[0] == 'CC':
                    gc_params0 = str(reader[1])
                    gc_params1 = float(reader[2])
                    #gc_params2 = str(reader[3])

                elif reader[0] == 'CR':
                    cr_params0 = str(reader[1])
                    cr_params1 = float(reader[2])

                elif reader[0] == 'KP':
                    kp_flag = True if reader[1] in true_lst else False

                elif reader[0] == 'ER':
                    er_params = map(float, reader[1:])
                elif reader[0] == 'GR':
                    gr_params = map(int, reader[1:])
                elif reader[0] == 'PV':
                    pv0_params = True if reader[1] in true_lst else False
                    pv1_params = str(reader[2])
                    pv2_params = int(reader[3])
                elif reader[0] == 'DA':
                    da0_params = str(reader[1])
                    da1_params = int(reader[2])

                #elif reader[0] == 'PS_p':
                    #iso_path = str(reader[1])

                elif reader[0] == 'BF':
                    bf_flag = True if reader[1] in true_lst else False
                    best_fit_algor = str(reader[2])
                    N_b = int(reader[3])
                elif reader[0] == 'RM':
                    flag_red_memb = str(reader[1])
                    min_prob = float(reader[2])
                elif reader[0] == 'PS_m':
                    m_rs = map(float, reader[1:])
                elif reader[0] == 'PS_a':
                    a_rs = map(float, reader[1:])
                elif reader[0] == 'PS_e':
                    e_rs = map(float, reader[1:])
                elif reader[0] == 'PS_d':
                    d_rs = map(float, reader[1:])
                elif reader[0] == 'SC':
                    IMF_name = str(reader[1])
                    tot_mass = float(reader[2])
                    f_bin = float(reader[3])
                    q_bin = float(reader[4])
                elif reader[0] == 'GA':
                    n_pop = int(reader[1])
                    n_gen = int(reader[2])
                    fdif = float(reader[3])
                    p_cross = float(reader[4])
                    cr_sel = str(reader[5])
                    p_mut = float(reader[6])
                    n_el = int(reader[7])
                    n_ei = int(reader[8])
                    n_es = int(reader[9])

    pl_params = [flag_make_plot, plot_frmt, plot_dpi]
    cr_params = [cr_params0, cr_params1]
    gc_params = [gc_params0, gc_params1]
    in_dirs = [input_dir, output_dir, done_dir]
    pv_params = [pv0_params, pv1_params, pv2_params]
    da_params = [da0_params, da1_params, input_dir]

    # Fix photometric system according to the CMD selected.
    #if cmd_select == 1 or cmd_select == 2:
        #sys_select = 'UBVI'
    #elif cmd_select == 3:
        #sys_select = 'WASH'

    # Fix isochrones location according to the CMD and set selected.
    if cmd_select in {1, 2}:
        text1 = 'ubvi'
    elif cmd_select in {3}:
        text1 = 'wash'
    if iso_select == 'MAR':
        text2 = 'marigo'
    elif iso_select == 'PAR':
        text2 = 'parsec'
    # Set iso_path according to the above values.
    iso_path = join(mypath + '/isochrones/' + 'iso_' + text1 + '_' + text2)

    # Fix metallicity and age step value since this is determined by the
    # stored isochrones.
    m_rs.append(0.0005)
    a_rs.append(0.05)
    ps_params = [iso_path, cmd_select, iso_select, m_rs, a_rs, e_rs, d_rs]

    bf_params = [bf_flag, best_fit_algor, N_b]
    sc_params = [IMF_name, tot_mass, f_bin, q_bin]
    ga_params = [n_pop, n_gen, fdif, p_cross, cr_sel, p_mut, n_el, n_ei, n_es]
    rm_params = [flag_red_memb, min_prob]

    # Fix magnitude and color names for the CMD axis
    if cmd_select == 1:
        x_ax, y_ax = '(B-V)', 'V'
    elif cmd_select == 2:
        x_ax, y_ax = '(V-I)', 'V'
    elif cmd_select == 3:
        x_ax, y_ax = '(C-{T_1})', '{T_1}'
    # Maximum and minimum axis values for the CMD plots.
    if cmd_select in {1, 2, 3}:
        # col_min col_max mag_min mag_max
        xy_minmax = [-1., 4., 7., 30.]
    # Store axes params.
    axes_params = [x_ax, y_ax, xy_minmax]

    return mode, in_dirs, gd_params, gc_params, cr_params, kp_flag, \
    er_params, gr_params, pv_params, da_params, ps_params, bf_params, \
    sc_params, ga_params, rm_params, pl_params, flag_move_file, axes_params