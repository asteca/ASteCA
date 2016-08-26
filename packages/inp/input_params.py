
import re


def CMD_phot_systs_filts():
    '''
    Dictionary that stores the names and column names for each
    filter defined in each photometric system, as presented in the CMD
    Girardi et al. service: http://stev.oapd.inaf.it/cgi-bin/cmd

    Capitalization of the filter names matters!
    '''
    cmd_systs = {
        '0': ('2mass_spitzer_wise', (
              'J', 'H', 'Ks', 'IRAC_3.6', 'IRAC_4.5', 'IRAC_5.8', 'IRAC_8.0',
              'MIPS_24', 'MIPS_70', 'MIPS_160', 'W1', 'W2', 'W3', 'W4')),
        '1': ('2mass', ('J', 'H', 'Ks')),
        '2': ('ogle_2mass_spitzer', (
              'U', 'B', 'V', 'I', 'J', 'H', 'Ks', 'IRAC_3.6', 'IRAC_4.5',
              'IRAC_5.8', 'IRAC_8.0', 'MIPS_24', 'MIPS_70', 'MIPS_160')),
        '3': ('2mass_spitzer_wise_washington_ddo51', (
              'J', 'H', 'Ks', '[3.6]', '[4.5]', '[5.8]', '[8.0]', 'W1', 'W2',
              'W3', 'W4', 'C', 'M', 'T1', 'T2', 'DDO51_finf')),
        '4': ('ubvrijhk', ('U', 'B', 'V', 'R', 'I', 'J', 'H', 'K')),
        '5': ('bessell', ('UX', 'BX', 'B', 'V', 'R', 'I', 'J', 'H', 'K', 'L',
                          "L'", 'M')),
        '6': ('akari', (
              'IRC_N2', 'IRC_N3', 'IRC_N4', 'MIRS_S7', 'MIRS_S9w', 'MIRS_S11',
              'MIRL_L15', 'MIRL_L18W', 'MIRL_L24', 'FIS_N60', 'FIS_WIDE-L',
              'FIS_WIDE-S', 'FIS_N160')),
        '7': ('batc', (
              'U', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'm',
              'n', 'o', 'p', 't')),
        '8': ('megacam_wircam', ('u*', "g'", "r'", "i'", "z'", 'Y', 'J', 'H',
                                 'Ks')),
        '9': ('wircam', ('Y', 'J', 'H', 'Ks')),
        '10': ('megacam', ('u*', "g'", "r'", "i'", "z'")),
        '11': ('ciber', (
               'l06', 'l07', 'l08', 'l09', 'l10', 'l11', 'l12', 'l13', 'l14',
               'l15', 'l16', 'l17', 'l18', 'l19', 'l20', 'I', 'H', 'NBS',
               'DSSI2', 'J', 'H', 'Ks')),
        '12': ('dcmc', ('I', 'J', 'H', 'Ks')),
        '13': ('decam', ('DECam-u', 'DES-g', 'DES-r', 'DES-i', 'DES-z',
                         'DES-Y')),
        '14': ('decam_vista', ('DECam-u', 'DES-g', 'DES-r', 'DES-i', 'DES-z',
                               'DES-Y', 'Z', 'Y', 'J', 'H', 'Ks')),
        '15': ('denis', ('I', 'J', 'Ks')),
        '16': ('dmc14', (
               'band1', 'band2', 'band3', 'band4', 'band5', 'band6', 'band7',
               'band8', 'band9', 'band10', 'band11', 'band12', 'band13',
               'band14')),
        '17': ('dmc15', (
               'band1', 'band2', 'band3', 'band4', 'band5', 'band6', 'band7',
               'band8', 'band9', 'band10', 'band11', 'band12', 'band13',
               'band14', 'band15')),
        '18': ('eis', (
               'u_wfi841', 'b_wfi842', 'v_wfi843', 'r_wfi844', 'i_wfi845',
               'z_wfi846', 'j_sofi', 'h_sofi', 'k_sofi')),
        '19': ('wfi', ('U841', 'U877', 'B842', 'B878', 'V843', 'R844', 'I845',
                       'I879')),
        '20': ('wfi_sofi', ('U841', 'U877', 'B842', 'B878', 'V843', 'R844',
                            'I845', 'I879', 'j_sofi', 'h_sofi', 'k_sofi')),
        '21': ('wfi2', ('U_35060', 'U_38', 'B_123', 'B', 'V', 'R', 'I_EIS',
                        'I')),
        '22': ('galex', ('FUV', 'NUV', 'U', 'B', 'V')),
        '23': ('galex_sloan', ('FUV', 'NUV', 'u', 'g', 'r', 'i', 'z')),
        '24': ('gaia', ('G', 'G_BP', 'G_RP')),
        '25': ('UVbright', (
               'F140LP', 'F175Wfoc', 'F275Wfoc', 'F200LP', 'F218W', 'F225W',
               'F275W', 'F300X', 'F336W', 'F25QTZ', 'FUV', 'NUV', 'UVW2',
               'UVM2', 'UVW1', 'u', 'U', 'B', 'V')),
        '26': ('acs_hrc', (
               'F220W', 'F250W', 'F330W', 'F344N', 'F435W', 'F475W', 'F550M',
               'F555W', 'F606W', 'F625W', 'F658N', 'F660N', 'F775W', 'F814W',
               'F850LP', 'F892N')),
        '27': ('acs_wfc', (
               'F435W', 'F475W', 'F502N', 'F550M', 'F555W', 'F606W', 'F625W',
               'F658N', 'F660N', 'F775W', 'F814W', 'F850LP')),
        '28': ('nicmosab', ('F110W', 'F160W', 'F205W')),
        '29': ('nicmosst', ('F110W', 'F160W', 'F205W')),
        '30': ('nicmosvega', ('F110W', 'F160W', 'F205W', 'F187N', 'F190N')),
        '31': ('stis', ('F25QTZ',)),
        '32': ('wfc3ir', ('F105W', 'F110W', 'F125W', 'F140W', 'F160W')),
        '33': ('wfc3uvis1', (
               'F200LP1', 'F218W1', 'F225W1', 'F275W1', 'F336W1', 'F350LP1',
               'F390W1', 'F438W1', 'F475W1', 'F555W1', 'F600LP1', 'F606W1',
               'F625W1', 'F775W1', 'F814W1', 'F850LP1')),
        '34': ('wfc3uvis2', (
               'F200LP2', 'F218W2', 'F225W2', 'F275W2', 'F336W2', 'F350LP2',
               'F390W2', 'F438W2', 'F475W2', 'F555W2', 'F600LP2', 'F606W2',
               'F625W2', 'F775W2', 'F814W2', 'F850LP2')),
        '35': ('wfc3_wideverywide', (
               'F200LP1', 'F218W1', 'F225W1', 'F275W1', 'F300X', 'F336W',
               'F350LP', 'F390W', 'F438W', 'F475W', 'F475X', 'F555W', 'F600LP',
               'F606W', 'F625W', 'F775W', 'F814W', 'F850LP', 'F105W', 'F110W',
               'F125W', 'F140W', 'F160W')),
        '36': ('wfc3_verywide', ('F200LP1', 'F300X', 'F350LP', 'F475X',
                                 'F600LP', 'F850LP')),
        '37': ('wfc3_medium', (
               'F390M', 'F410M', 'FQ422M', 'F467M', 'F547M', 'F621M', 'F689M',
               'F763M', 'F845M', 'F098M', 'F127M', 'F139M', 'F153M')),
        '38': ('wfc3_wide', (
               'F218W1', 'F225W1', 'F275W1', 'F336W', 'F390W', 'F438W',
               'F475W', 'F555W', 'F606W', 'F625W', 'F775W', 'F814W', 'F105W',
               'F110W', 'F125W', 'F140W', 'F160W')),
        '39': ('wfpc2', (
               'F170W', 'F218W', 'F255W', 'F300W', 'F336W', 'F380W', 'F439W',
               'F450W', 'F555W', 'F569W', 'F606W', 'F622W', 'F675W', 'F702W',
               'F791W', 'F814W', 'F850LP')),
        '40': ('int_wfc', ('u', 'g', 'r', 'i', 'z')),
        '41': ('iphas', ('gR', 'gI', 'Ha')),
        '42': ('kepler', ('Kepler', 'g', 'r', 'i', 'z', 'DDO51_finf')),
        '43': ('kepler_2mass', ('Kepler', 'g', 'r', 'i', 'z', 'DDO51_finf',
                                'J', 'H', 'Ks')),
        '44': ('lbt_lbc', (
               'blue_Uspec', 'blue_U', 'blue_B', 'blue_g', 'blue_V', 'blue_r',
               'red_V', 'red_i', 'red_I', 'red_r', 'red_R', 'red_z', 'red_Y')),
        '45': ('lsst', ('u', 'g', 'r', 'i', 'z', 'Y')),
        '46': ('noao_ctio_mosaic2', ('u', 'U', 'B', 'g', 'V', 'i', 'I', 'r',
                                     'R', 'z')),
        '47': ('ogle', ('U', 'B', 'V', 'I')),
        '48': ('panstarrs1', ('gP1', 'rP1', 'iP1', 'zP1', 'yP1', 'wP1')),
        '49': ('splus', ('uJAVA', 'F378', 'F395', 'F410', 'F430', 'gSDSS',
                         'F515', 'rSDSS', 'F660', 'iSDSS', 'F861', 'zSDSS')),
        '50': ('sloan', ('u', 'g', 'r', 'i', 'z')),
        '51': ('sloan_2mass', ('u', 'g', 'r', 'i', 'z', 'J', 'H', 'Ks')),
        '52': ('sloan_ukidss', ('u', 'g', 'r', 'i', 'z', 'Z', 'Y', 'J', 'H',
                                'K')),
        '53': ('swift_uvot', ('UVW2', 'UVM2', 'UVW1', 'u')),
        '54': ('spitzer', ('IRAC_3.6', 'IRAC_4.5', 'IRAC_5.8', 'IRAC_8.0',
                           'MIPS_24', 'MIPS_70', 'MIPS_160')),
        '55': ('stroemgren', ('V', 'u', 'v', 'b', 'y', 'Hb_w', 'Hb_n')),
        '56': ('suprimecam', ('B', 'V', 'R', 'I', 'g', 'r', 'i', 'z')),
        '57': ('TESS_2mass_kepler', ('TESS', 'J', 'H', 'Ks', 'Kepler', 'g',
                                     'r', 'i', 'z', 'DDO51_finf')),
        '58': ('tycho2', ('B_T', 'V_T')),
        '59': ('ukidss', ('Z', 'Y', 'J', 'H', 'K')),
        '60': ('visir', ('F606W', 'F814W', 'Ks', 'PAH1')),
        '61': ('vista', ('Z', 'Y', 'J', 'H', 'Ks')),
        '62': ('vphas', ('u', 'g', 'r', 'i', 'z', 'Ha')),
        '63': ('vst_omegacam', ('u', 'g', 'r', 'i', 'z')),
        '64': ('vilnius', ('U', 'P', 'X', 'Y', 'Z', 'V', 'S')),
        '65': ('washington', ('C', 'M', 'T1', 'T2', 'B', 'V', 'R', 'I')),
        '66': ('washington_ddo51', ('C', 'M', 'T1', 'T2', 'DDO51_finf',
                                    'DDO51_f7.9', 'DDO51_f2.9'))
    }

    return cmd_systs


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
                    id_coords = map(int, reader[1:-1]) + [reader[-1]]
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

                # Membership based removal parameters.
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
                    evol_track = str(reader[1])

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
                    print("  WARNING: Unknown '{}' ID found in line {}\n"
                          "  of '{}' file.\n").format(
                        reader[0], l + 1, pars_f_name)

    # Pack params in lists.
    pl_params = [flag_make_plot, plot_frmt, plot_dpi]
    gh_params = [gh_params0, gh_params1]
    cr_params = [cr_params]
    pv_params = [pv0_params, pv1_params]
    da_params = [da0_params, da1_params]

    # Dictionary with data on the CMD service photometric systems.
    cmd_systs = CMD_phot_systs_filts()

    # Map evolutionary tracks selection to proper names, and name of the folder
    # where they should be stored.
    cmd_evol_tracks = {
        'PAR12C': ('parsec12C', 'COLIBRI PR16'),
        'PAR12': ('parsec12', 'PARSEC v1.2S'),
        'PAR10': ('parsec10', 'PARSEC v1.0'),
        'PAR11': ('parsec11', 'PARSEC v1.1'),
        'MAR08A': ('marigo08A', 'Marigo (2008, Case A)'),
        'MAR08B': ('marigo08B', 'Marigo (2008, Case B)'),
        'MAR08': ('marigo08', 'Marigo (2008)'),
        'GIR02': ('girardi02', 'Girardi (2002)')}

    # Store photometric system params in lists.
    par_ranges = [m_rs, a_rs, e_rs, d_rs, mass_rs, bin_rs]

    # Store GA params in lists.
    bf_params = [best_fit_algor, lkl_method, bin_method, N_b]
    sc_params = [m_high, bin_mr]
    ga_params = [n_pop, n_gen, fdif, p_cross, cr_sel, p_mut, n_el, n_ei, n_es]
    rm_params = [mode_red_memb, local_bin, min_prob]

    # Define tuple of accepted binning methods.
    bin_methods = ('blocks', 'knuth', 'scott', 'freedman', 'sturges',
                   'sqrt', 'bb')
    # Accepted IMF functions.
    imf_funcs = ('chabrier_2001_exp', 'chabrier_2001_log', 'kroupa_1993',
                 'kroupa_2002')

    # Store all read parameters in dictionary.
    pd = {
        'up_flag': up_flag, 'mode': mode, 'done_dir': done_dir,
        'id_coords': id_coords, 'id_mags': id_mags, 'id_cols': id_cols,
        'evol_track': evol_track, 'par_ranges': par_ranges,
        'cmd_evol_tracks': cmd_evol_tracks, 'bf_flag': bf_flag,
        'cmd_systs': cmd_systs, 'bin_methods': bin_methods,
        'imf_funcs': imf_funcs, 'IMF_name': IMF_name,
        'gh_params': gh_params,
        'cr_params': cr_params, 'kp_flag': kp_flag, 'im_flag': im_flag,
        'er_params': er_params, 'fr_number': fr_number,
        'pv_params': pv_params, 'da_params': da_params,
        'bf_params': bf_params, 'sc_params': sc_params,
        'ga_params': ga_params, 'rm_params': rm_params,
        'pl_params': pl_params, 'flag_move_file': flag_move_file,
        'flag_back_force': flag_back_force}

    # Return parameters dictionary 'pd'
    return pd
