# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 14:03:44 2014

@author: gabriel
"""

from os.path import join
import re


def char_remove(in_lst):
    '''
    Correctly convert input data for parameters ranges to lists.
    '''
    lst = []
    # If input list is empty, return empty list.
    if in_lst[1:]:
        l0 = []
        if in_lst[1][0] in {'[', '(', '{'}:
            # Remove non-numeric characters and append numbers as floats.
            l0.append([float(i) for i in re.findall('[0-9.]+',
                      str(in_lst[1:]))])
            # Store indicating that this is a list of values.
            lst = ['l', map(float, l0[0])]
        else:
            # Store indicating that this is a range of values.
            lst = ['r', map(float, in_lst[1:4])]

    return lst


def init(mypath, pars_f_path):
    '''
    This function reads the input data parameters stored in the
    'params_input_XX.dat' file and returns them packaged for each function to
    use.
    '''

    # Declare globals.
    global up_flag, mode, done_dir, gd_params, gh_params, cr_params, kp_flag,\
        im_flag, er_params, fr_number, pv_params, da_params, ps_params, \
        bf_params, sc_params, ga_params, rm_params, pl_params, \
        flag_move_file, axes_params, flag_back_force

    # Accept these variations of the 'true' flag.
    true_lst = ('True', 'true', 'TRUE')

    # Read data from file.
    with open(pars_f_path, "r") as f_dat:

        # Iterate through each line in the file.
        for l, line in enumerate(f_dat):

            if not line.startswith("#") and line.strip() != '':
                reader = line.split()

                # Updater.
                if reader[0] == 'UP':
                    up_flag = True if reader[1] in true_lst else False

                # Mode.
                elif reader[0] == 'MO':
                    mode = str(reader[1])

                # Input data parameters.
                elif reader[0] == 'PD':
                    gd_params = map(int, reader[1:])
                elif reader[0] == 'PX':
                    gd_params.append(str(reader[1]))
                elif reader[0] == 'CMD':
                    cmd_select = int(reader[1])

                # Output parameters.
                elif reader[0] == 'FB':
                    flag_back_force = True if reader[1] in true_lst else False
                elif reader[0] == 'MP':
                    flag_make_plot = True if reader[1] in true_lst else False
                    plot_frmt = str(reader[2])
                    plot_dpi = int(reader[3])
                elif reader[0] == 'MF':
                    flag_move_file = True if reader[1] in true_lst else False
                    done_dir = str(reader[2])

                # Structure functions parameters.
                elif reader[0] == 'CH':
                    gh_params0 = str(reader[1])
                    gh_params1 = float(reader[2])
                elif reader[0] == 'CR':
                    cr_params = str(reader[1])
                elif reader[0] == 'KP':
                    kp_flag = True if reader[1] in true_lst else False
                elif reader[0] == 'GR':
                    try:
                        fr_number = int(reader[1])
                    except:
                        fr_number = str(reader[1])

                # Photometric functions parameters.
                elif reader[0] == 'ER':
                    er_params = [str(reader[1])] + map(float, reader[2:])
                elif reader[0] == 'IM':
                    im_flag = True if reader[1] in true_lst else False
                elif reader[0] == 'PV':
                    pv0_params = str(reader[1])
                    pv1_params = int(reader[2])
                elif reader[0] == 'DA':
                    da0_params = str(reader[1])
                    da1_params = int(reader[2])

                # Reduce membership parameters.
                elif reader[0] == 'RM':
                    mode_red_memb = str(reader[1])
                    local_bin = str(reader[2])
                    min_prob = float(reader[3])

                # Cluster parameters assignation.
                elif reader[0] == 'BF':
                    bf_flag = True if reader[1] in true_lst else False
                    best_fit_algor = str(reader[2])
                    lkl_method = str(reader[3])
                    bin_method = str(reader[4])
                    N_b = int(reader[5])
                elif reader[0] == 'PS':
                    iso_select = str(reader[1])

                # Synthetic cluster parameters.
                elif reader[0] == 'SC':
                    IMF_name = str(reader[1])
                    m_high = float(reader[2])
                elif reader[0] == 'PS_m':
                    m_rs = char_remove(reader)
                elif reader[0] == 'PS_a':
                    a_rs = char_remove(reader)
                elif reader[0] == 'PS_e':
                    e_rs = char_remove(reader)
                elif reader[0] == 'PS_d':
                    d_rs = char_remove(reader)
                elif reader[0] == 'TM':
                    mass_rs = char_remove(reader)
                elif reader[0] == 'BI':
                    bin_mr = float(reader[1])
                    # Use [1:] since the mass ratio is located before the
                    # range.
                    bin_rs = char_remove(reader[1:])

                # Genetic algorithm parameters.
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

                else:
                    # Get parameters file name from path.
                    pars_f_name = pars_f_path.split('/')[-1]
                    print ("  WARNING: Unknown '{}' ID found in line {}\n"
                           "  of '{}' file.\n").format(
                           reader[0], l + 1, pars_f_name)

    # Pack params in lists.
    pl_params = [flag_make_plot, plot_frmt, plot_dpi]
    gh_params = [gh_params0, gh_params1]
    cr_params = [cr_params]
    pv_params = [pv0_params, pv1_params]
    da_params = [da0_params, da1_params]

    # Fix isochrones location according to the CMD and set selected.
    text1, text2 = 'none', 'none'
    text1 = 'parsec' + iso_select[-2:]
    if cmd_select in {1, 2, 3}:
        text2 = 'ubvrijhk'
    elif cmd_select in {4}:
        text2 = 'washington'
    elif cmd_select in {5, 6, 7}:
        text2 = '2mass'
    elif cmd_select in {8, 9}:
        text2 = 'sloan'
    elif cmd_select in {10, 11, 12}:
        text2 = 'stroemgren'
    # Set iso_path according to the above values.
    iso_path = join(mypath + '/isochrones/' + text1 + '_' + text2)

    # Fix magnitude and color names for the CMD axis.
    # m_1 is the y axis magnitude, m_2 is the magnitude used to obtain the
    # color index and the third value in each key indicates how the color
    # is to be formed, e.g: '12' means (m_1 - m_2)
    cmds_dic = {1: ('V', 'B', 21), 2: ('V', 'I', 12), 3: ('V', 'U', 21),
                4: ('{T_1}', 'C', 21), 5: ('J', 'H', 12), 6: ('H', 'J', 21),
                7: ('K', 'H', 21), 8: ('g', 'u', 21), 9: ('g', 'r', 12),
                10: ('y', 'b', 21), 11: ('y', 'v', 21), 12: ('y', 'u', 21)}
    m_1, m_2, m_ord = cmds_dic[cmd_select]
    # Store axes params.
    axes_params = [m_1, m_2, m_ord]

    # Store photometric system params in lists.
    par_ranges = [m_rs, a_rs, e_rs, d_rs, mass_rs, bin_rs]
    ps_params = [iso_path, cmd_select, iso_select, par_ranges]

    # Store GA params in lists.
    bf_params = [bf_flag, best_fit_algor, lkl_method, bin_method, N_b]
    sc_params = [IMF_name, m_high, bin_mr]
    ga_params = [n_pop, n_gen, fdif, p_cross, cr_sel, p_mut, n_el, n_ei, n_es]
    rm_params = [mode_red_memb, local_bin, min_prob]
