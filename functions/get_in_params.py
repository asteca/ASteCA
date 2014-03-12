# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 14:03:44 2014

@author: gabriel
"""

from os.path import join
#import re


def get_in_params(mypath):
    '''
    This function reads the input data parameters stored in the
    'ocaat_input.dat' file and returns them packaged for each function to use.
    '''

# Allows to work with columns data files.

    data_file = join(mypath, 'ocaat_input.dat')

    with open(data_file, mode="r") as f_dat:

        # Iterate through each line in the file.
        for line in f_dat:

            if not line.startswith("#") and line.strip() != '':
                reader = line.split()

                # Read folder paths where clusters are stored.
                if reader[0] == 'MO':
                    mode = str(reader[1])

                elif reader[0] == 'CP_i':
                    input_dir = str(reader[1])
                elif reader[0] == 'CP_o':
                    output_dir = str(reader[1])
                elif reader[0] == 'CP_d':
                    done_dir = str(reader[1])

                elif reader[0] == 'MP':
                    flag_make_plot = True if reader[1] == 'True' else False

                elif reader[0] == 'PD':
                    gd_params = map(int, reader[1:])

                elif reader[0] == 'CMD':
                    cmd_select = int(reader[1])

                elif reader[0] == 'PS':
                    iso_select = str(reader[1])

                elif reader[0] == 'CC':
                    gc_params = map(float, reader[1:])
                elif reader[0] == 'BR':
                    br_params = map(float, reader[1:])
                elif reader[0] == 'CR':
                    cr_params = map(float, reader[1:])
                elif reader[0] == 'ER':
                    er_params = map(float, reader[1:])
                elif reader[0] == 'GR':
                    gr_params = map(int, reader[1:])
                elif reader[0] == 'PV':
                    pv0_params = True if reader[1] == 'True' else False
                    pv1_params = str(reader[2])
                    pv2_params = int(reader[3])
                elif reader[0] == 'DA':
                    da0_params = str(reader[1])
                    da1_params = int(reader[2])
                    da2_params = int(reader[3])

                #elif reader[0] == 'PS_p':
                    #iso_path = str(reader[1])

                elif reader[0] == 'PS_m':
                    m_rs = map(float, reader[1:])
                elif reader[0] == 'PS_a':
                    a_rs = map(float, reader[1:])
                elif reader[0] == 'PS_e':
                    e_rs = map(float, reader[1:])
                elif reader[0] == 'PS_d':
                    d_rs = map(float, reader[1:])
                elif reader[0] == 'BF':
                    bf_flag = True if reader[1] == 'True' else False
                    best_fit_algor = str(reader[2])
                    N_b = int(reader[3])
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
                elif reader[0] == 'MF':
                    flag_move_file = True if reader[1] == 'True' else False

                #elif reader[0] == 'XA':
                    #x_ax = re.search(r'"(.*)"', line).groups()[0]
                #elif reader[0] == 'YA':
                    #y_ax = re.search(r'"(.*)"', line).groups()[0]

                elif reader[0] == 'MM':
                    xy_minmax = map(float, reader[1:])

    in_dirs = [input_dir, output_dir, done_dir]
    pv_params = [pv0_params, pv1_params, pv2_params]
    da_params = [da0_params, da1_params, da2_params, input_dir]

    # Fix photometric system according to the CMD selected.
    if cmd_select == 1 or cmd_select == 2:
        sys_select = 'UBVI'
    elif cmd_select == 3:
        sys_select = 'WASH'

    # Fix isochrones location according to the CMD and set selected.
    if sys_select == 'UBVI':
        text1 = 'ubvi'
    elif sys_select == 'WASH':
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
    ps_params = [iso_path, sys_select, iso_select, m_rs, a_rs, e_rs, d_rs]

    bf_params = [bf_flag, best_fit_algor, N_b]
    sc_params = [IMF_name, tot_mass, f_bin, q_bin]
    ga_params = [n_pop, n_gen, fdif, p_cross, cr_sel, p_mut, n_el, n_ei, n_es]

    # Fix magnitude and color names for the CMD axis
    if cmd_select == 1:
        x_ax, y_ax = '(B-V)', 'V'
    elif cmd_select == 2:
        x_ax, y_ax = '(V-I)', 'V'
    elif cmd_select == 3:
        x_ax, y_ax = '(C-{T_1})', '{T_1}'
    axes_params = [x_ax, y_ax, xy_minmax]

    return mode, in_dirs, gd_params, gc_params, br_params, cr_params, \
    er_params, gr_params, pv_params, da_params, ps_params, bf_params, \
    sc_params, ga_params, flag_make_plot, flag_move_file, axes_params