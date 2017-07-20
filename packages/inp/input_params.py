
import re
import CMD_phot_systs


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
            if len(in_lst[1:4]) < 3:
                # Not enough values to define a range. Store as list of values.
                lst = ['l', map(float, in_lst[1:4])]
            else:
                # Store indicating that this is a range of values.
                lst = ['r', map(float, in_lst[1:4])]

    return lst


def main(mypath, pars_f_path):
    '''
    This function reads the input data parameters stored in the
    'params_input_XX.dat' file and returns a dictionary.
    '''

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

                # Set global mode (i.e, for all clusters processed).
                elif reader[0] == 'MO':
                    mode = str(reader[1])

                # Input data parameters.
                elif reader[0] == 'PI':
                    id_coords = reader[1:]
                elif reader[0] == 'PM':
                    id_mags = reader[1:]
                elif reader[0] == 'PC':
                    id_cols = reader[1:]

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
                    radius_method = str(reader[1])
                elif reader[0] == 'KP':
                    kp_flag = True if reader[1] in true_lst else False
                elif reader[0] == 'GR':
                    try:
                        fr_number = int(reader[1])
                    except ValueError:
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

                # Membership based removal parameters.
                elif reader[0] == 'RM':
                    mode_fld_clean = str(reader[1])
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
                    evol_track = str(reader[1])

                # Synthetic cluster parameters.
                elif reader[0] == 'MM':
                    try:
                        max_mag = float(reader[1])
                    except ValueError:
                        max_mag = str(reader[1])
                elif reader[0] == 'SC':
                    IMF_name = str(reader[1])
                    m_high = float(reader[2])
                elif reader[0] == 'RV':
                    R_V = float(reader[1])
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
                    print("  WARNING: Unknown '{}' ID found in line {}\n"
                          "  of '{}' file.\n").format(
                        reader[0], l + 1, pars_f_name)

    # Pack params in lists.
    pl_params = [flag_make_plot, plot_frmt, plot_dpi]
    gh_params = [gh_params0, gh_params1]
    pv_params = [pv0_params, pv1_params]
    da_params = [da0_params, da1_params]

    # Dictionary with data on the CMD service photometric systems.
    cmd_systs = CMD_phot_systs.main()

    # Map evolutionary tracks selection to proper names, and name of the folder
    # where they should be stored.
    cmd_evol_tracks = {
        # 'PAR12C': ('parsec12C', 'COLIBRI PR16'),
        'PAR12': ('parsec12', 'PARSEC v1.2S'),
        'PAR10': ('parsec10', 'PARSEC v1.0'),
        'PAR11': ('parsec11', 'PARSEC v1.1'),
        'MAR08A': ('marigo08A', 'Marigo (2008, Case A)'),
        'MAR08B': ('marigo08B', 'Marigo (2008, Case B)'),
        'MAR08': ('marigo08', 'Marigo (2008)'),
        'GIR02': ('girardi02', 'Girardi (2002)')}

    # Photometric system parameters.
    par_ranges = [m_rs, a_rs, e_rs, d_rs, mass_rs, bin_rs]

    # Store GA params in lists.
    bf_params = [best_fit_algor, lkl_method, bin_method, N_b]
    ga_params = [n_pop, n_gen, fdif, p_cross, cr_sel, p_mut, n_el, n_ei, n_es]
    rm_params = [mode_fld_clean, local_bin, min_prob]

    # Accepted field stars removal methods.
    fld_rem_methods = ('local', 'n_memb', 'mp_05', 'top_h', 'man', 'skip')
    # Accepted binning methods.
    bin_methods = ('fixed', 'auto', 'fd', 'doane', 'scott', 'rice', 'sqrt',
                   'sturges', 'knuth', 'blocks')
    # Accepted IMF functions.
    imf_funcs = ('chabrier_2001_exp', 'chabrier_2001_log', 'kroupa_1993',
                 'kroupa_2002')

    pd = {
        'up_flag': up_flag, 'mode': mode, 'done_dir': done_dir,
        'id_coords': id_coords, 'id_mags': id_mags, 'id_cols': id_cols,
        'radius_method': radius_method,
        'evol_track': evol_track, 'par_ranges': par_ranges,
        'cmd_evol_tracks': cmd_evol_tracks, 'bf_flag': bf_flag,
        'cmd_systs': cmd_systs, 'fld_rem_methods': fld_rem_methods,
        'bin_methods': bin_methods, 'imf_funcs': imf_funcs,
        'max_mag': max_mag, 'IMF_name': IMF_name, 'm_high': m_high,
        'R_V': R_V, 'bin_mr': bin_mr,
        'gh_params': gh_params,
        'kp_flag': kp_flag, 'im_flag': im_flag,
        'er_params': er_params, 'fr_number': fr_number,
        'pv_params': pv_params, 'da_params': da_params,
        'bf_params': bf_params, 'ga_params': ga_params,
        'rm_params': rm_params, 'pl_params': pl_params,
        'flag_move_file': flag_move_file, 'flag_back_force': flag_back_force}

    # Return parameters dictionary.
    return pd
