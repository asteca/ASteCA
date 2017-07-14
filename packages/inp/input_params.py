
import re


def CMD_phot_systs_filts():
    '''
    Dictionary that stores the names and column names for each
    filter defined in each photometric system, as presented in the CMD
    Girardi et al. service: http://stev.oapd.inaf.it/cgi-bin/cmd

    Capitalization of the filter names matters!

    Last update: v3.0 (13th July 2017)
    '''
    cmd_systs = {
        # 2MASS + Spitzer (IRAC+MIPS)
        '0': ('2mass_spitzer', ('J', 'H', 'Ks', '[3.6]', '[4.5]', '[5.8]',
              '[8.0]', '[24]', '[70]', '[160]'), ()),
        # 2MASS + Spitzer (IRAC+MIPS) + WISE
        '1': ('2mass_spitzer_wise', ('J', 'H', 'Ks', '[3.6]', '[4.5]', '[5.8]',
              '[8.0]', '[24]', '[70]', '[160]', 'W1', 'W2', 'W3', 'W4'),
              (12329.79, 16395.59, 21522.05, 35242.46, 44540.45, 56095.61,
               77024.24, 233597.72, 693695.94, 1537007.38, 33159.26,
               45611.97, 107878.20, 219085.73)),
        # 2MASS JHKs
        '2': ('2mass', ('J', 'H', 'Ks'), (12329.79, 16395.59, 21522.05)),
        # OGLE + 2MASS + Spitzer (IRAC+MIPS)
        '3': ('ogle_2mass_spitzer', (
              'U', 'B', 'V', 'I', 'J', 'H', 'Ks', 'IRAC_3.6', 'IRAC_4.5',
              'IRAC_5.8', 'IRAC_8.0', 'MIPS_24', 'MIPS_70', 'MIPS_160'),
              (3694.87, 4424.79, 5443.37, 8320.74, 12329.79, 16395.59,
               21522.05, 35242.46, 44540.45, 56095.61, 77024.24, 233597.73,
               693695.94, 1537007.38)),
        # 2MASS+Spitzer+WISE+Washington+DDO51
        '4': ('2mass_spitzer_wise_washington_ddo51', (
              'J', 'H', 'Ks', '[3.6]', '[4.5]', '[5.8]', '[8.0]', 'W1', 'W2',
              'W3', 'W4', 'C', 'M', 'T1', 'T2', 'DDO51_finf'), (12329.79,
              16395.59, 21522.05, 35242.46, 44540.45, 56095.61, 77024.24,
              33159.26, 45611.97, 107878.20, 219085.73, 3982.34, 5120.46,
              6420.73, 8077.89, 5143.24)),
        # UBVRIJHK (cf. Maiz-Apellaniz 2006 + Bessell 1990)
        '5': ('ubvrijhk', ('U', 'B', 'V', 'R', 'I', 'J', 'H', 'K'),
              (3641.89, 4460.62, 5501.70, 6557.09, 8036.57, 12314.46, 16369.53,
               21937.19)),
        # UBVRIJHKLMN (cf. Bessell 1990 + Bessell &amp; Brett 1988)
        '6': ('bessell', ('UX', 'BX', 'B', 'V', 'R', 'I', 'J', 'H', 'K', 'L',
              "L'", 'M'), (3627.36, 4478.30, 4459.58, 5500.02, 6498.73,
              8010.27, 12269.18, 16324.85, 21874.56, 34518.32, 38017.59,
              47174.98)),
        # AKARI
        '7': ('akari', (
              'IRC_N2', 'IRC_N3', 'IRC_N4', 'MIRS_S7', 'MIRS_S9w', 'MIRS_S11',
              'MIRL_L15', 'MIRL_L18W', 'MIRL_L24', 'FIS_N60', 'FIS_WIDE-L',
              'FIS_WIDE-S', 'FIS_N160'), (23187.43, 31640.68, 42882.73,
              70587.31, 84874.40, 103435.73, 154595.31, 181052.36, 228013.25,
              639610.06, 1430084.25, 803536.31, 1603268.38)),
        # BATC
        '8': ('batc', ('U', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',
              'k', 'm', 'n', 'o', 'p', 't'), (3680.68, 3382.66, 3908.97,
              4198.95, 4554.59, 4872.27, 5248.67, 5783.22, 6072.04, 6709.81,
              7009.71, 7526.40, 8023.69, 8516.15, 9171.69, 9722.47, 6596.29)),
        # CFHT Megacam + Wircam (all ABmags)
        '9': ('megacam_wircam', ('u*', "g'", "r'", "i'", "z'", 'Y', 'J', 'H',
              'Ks'), (3880.09, 4906.14, 6251.88, 7671.34, 8850.62, 10244.67,
              12523.13, 16222.89, 21386.52)),
        # CFHT Wircam
        '10': ('wircam', ('Y', 'J', 'H', 'Ks'), (10244.67, 12523.13, 16222.89,
               21386.52)),
        # CFHT/Megacam u*g'r'i'z'
        '11': ('megacam', ('u*', "g'", "r'", "i'", "z'"), (3880.09, 4906.14,
               6251.88, 7671.34, 8850.62)),
        # CIBER
        '12': ('ciber', (
               'l06', 'l07', 'l08', 'l09', 'l10', 'l11', 'l12', 'l13', 'l14',
               'l15', 'l16', 'l17', 'l18', 'l19', 'l20', 'I', 'H', 'NBS',
               'DSSI2', 'J', 'H', 'Ks'), (5998.55, 6993.81, 7990.56, 8992.58,
               9991.60, 10990.98, 11992.04, 12992.07, 13991.63, 14992.48,
               15989.14, 16986.70, 17987.32, 18988.57, 19990.05, 10993.81,
               15414.76, 8530.00, 8337.74, 12329.79, 16395.59, 21522.05)),
        # DECAM (ABmags)
        '13': ('decam', ('DECam-u', 'DES-g', 'DES-r', 'DES-i', 'DES-z',
               'DES-Y'), (3610.03, 4791.97, 6405.90, 7814.44, 9241.46,
               10079.67)),
        # DECAM ugrizY (ABmags) + VISTA ZYJHKs (Vegamags)
        '14': ('decam_vista', ('DECam-u', 'DES-g', 'DES-r', 'DES-i', 'DES-z',
               'DES-Y', 'Z', 'Y', 'J', 'H', 'Ks'), (3610.03, 4791.97, 6405.90,
               7814.44, 9241.46, 10079.67, 8762.58, 10239.64, 12543.63,
               16372.50, 21342.85)),
        # DENIS
        '15': ('denis', ('I', 'J', 'Ks'), (7927.50, 12307.14, 21488.90)),
        # DMC 14 filters
        '16': ('dmc14', (
               'band1', 'band2', 'band3', 'band4', 'band5', 'band6', 'band7',
               'band8', 'band9', 'band10', 'band11', 'band12', 'band13',
               'band14'), (4223.63, 4609.60, 5022.91, 5413.68, 5770.05,
               6252.21, 6682.12, 7063.48, 7328.10, 7712.16, 7975.48, 8294.39,
               8603.58, 9088.43)),
        # DMC 15 filters
        '17': ('dmc15', (
               'band1', 'band2', 'band3', 'band4', 'band5', 'band6', 'band7',
               'band8', 'band9', 'band10', 'band11', 'band12', 'band13',
               'band14', 'band15'), (4223.65, 4609.59, 5022.89, 5413.67,
               5770.04, 6089.84, 6410.44, 6682.12, 7063.47, 7328.11, 7712.16,
               7975.49, 8294.39, 8603.59, 9088.43)),
        # ESO/EIS (WFI UBVRIZ + SOFI JHK)
        '18': ('eis', (
               'u_wfi841', 'b_wfi842', 'v_wfi843', 'r_wfi844', 'i_wfi845',
               'z_wfi846', 'j_sofi', 'h_sofi', 'k_sofi'), (3684.67, 4636.89,
               5387.70, 6492.74, 8599.65, 9597.25, 12302.85, 16426.90,
               21558.65)),
        # ESO/WFI
        '19': ('wfi', ('U841', 'U877', 'B842', 'B878', 'V843', 'R844', 'I845',
               'I879'), (3688.80, 3622.17, 4646.41, 4652.75, 5390.38,
               6491.33, 8600.62, 8036.86)),
        # ESO/WFI2
        '20': ('wfi2', ('U_35060', 'U_38', 'B_123', 'B', 'V', 'R', 'I_EIS',
               'I'), (3469.36, 3656.71, 4617.26, 4620.10, 5385.72, 6489.33,
               8024.69, 8589.78)),
        # GALEX FUV+NUV (Vegamag) + SDSS ugriz (ABmags)
        '21': ('galex_sloan', ('FUV', 'NUV', 'u', 'g', 'r', 'i', 'z'),
               (1722.15, 2556.69, 3587.16, 4769.82, 6179.99, 7485.86,
                8933.86)),
        # GALEX FUV+NUV + Johnson's UBV (Maiz-Apellaniz version), all Vegamags
        '22': ('galex', ('FUV', 'NUV', 'U', 'B', 'V'), (1722.15, 2556.69,
               3641.89, 4460.62, 5501.70)),
        # Gaia's G, G_BP and G_RP (Vegamags)
        '23': ('gaia', ('G', 'G_BP', 'G_RP'), (6418.68, 5328.58, 7799.00)),
        # HST+GALEX+Swift/UVOT UV filters
        '24': ('UVbright', (
               'F140LP', 'F175Wfoc', 'F275Wfoc', 'F200LP', 'F218W', 'F225W',
               'F275W', 'F300X', 'F336W', 'F25QTZ', 'FUV', 'NUV', 'UVW2',
               'UVM2', 'UVW1', 'u', 'U', 'B', 'V'), (3156.54, 3556.64, 3082.65,
               6078.51, 2498.78, 2572.12, 2832.00, 3157.55, 3382.75, 2862.35,
               1722.15, 2556.69, 3291.97, 2454.10, 3265.49, 3526.14, 3641.89,
               4460.62, 5501.70)),
        # HST/ACS HRC
        '25': ('acs_hrc', (
               'F220W', 'F250W', 'F330W', 'F344N', 'F435W', 'F475W', 'F550M',
               'F555W', 'F606W', 'F625W', 'F658N', 'F660N', 'F775W', 'F814W',
               'F850LP', 'F892N'), (2644.66, 2954.20, 3395.80, 3433.21,
               4399.30, 4844.32, 5585.36, 5381.44, 5928.98, 6301.95, 6584.81,
               6599.48, 7658.64, 8092.18, 9130.28, 8915.88)),
        # HST/ACS WFC
        '26': ('acs_wfc', (
               'F435W', 'F475W', 'F502N', 'F550M', 'F555W', 'F606W', 'F625W',
               'F658N', 'F660N', 'F775W', 'F814W', 'F850LP')),
        '27': ('wfi_sofi', ('U841', 'U877', 'B842', 'B878', 'V843', 'R844',
                            'I845', 'I879', 'j_sofi', 'h_sofi', 'k_sofi')),
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
        '50': ('sloan', ('u', 'g', 'r', 'i', 'z'),
               (3587.16, 4769.82, 6179.99, 7485.86, 8933.86)),
        '51': ('sloan_2mass', ('u', 'g', 'r', 'i', 'z', 'J', 'H', 'Ks')),
        '52': ('sloan_ukidss', ('u', 'g', 'r', 'i', 'z', 'Z', 'Y', 'J', 'H',
                                'K')),
        '53': ('swift_uvot', ('UVW2', 'UVM2', 'UVW1', 'u')),
        '54': ('spitzer', ('IRAC_3.6', 'IRAC_4.5', 'IRAC_5.8', 'IRAC_8.0',
                           'MIPS_24', 'MIPS_70', 'MIPS_160')),
        '55': ('stroemgren', ('V', 'u', 'v', 'b', 'y', 'Hb_w', 'Hb_n'),
               (5523.22, 3473.64, 4115.47, 4673.30, 5478.71, 4880.05,
                4855.42)),
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
    cr_params = [cr_params]
    pv_params = [pv0_params, pv1_params]
    da_params = [da0_params, da1_params]

    # Dictionary with data on the CMD service photometric systems.
    cmd_systs = CMD_phot_systs_filts()

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
    rm_params = [mode_red_memb, local_bin, min_prob]

    # Define tuple of accepted binning methods.
    bin_methods = ('blocks', 'knuth', 'bb', 'auto', 'fd', 'doane', 'scott',
                   'rice', 'sturges', 'sqrt', 'man')
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
        'imf_funcs': imf_funcs, 'IMF_name': IMF_name, 'm_high': m_high,
        'R_V': R_V, 'bin_mr': bin_mr,
        'gh_params': gh_params,
        'cr_params': cr_params, 'kp_flag': kp_flag, 'im_flag': im_flag,
        'er_params': er_params, 'fr_number': fr_number,
        'pv_params': pv_params, 'da_params': da_params,
        'bf_params': bf_params,
        'ga_params': ga_params, 'rm_params': rm_params,
        'pl_params': pl_params, 'flag_move_file': flag_move_file,
        'flag_back_force': flag_back_force}

    # Return parameters dictionary 'pd'
    return pd
