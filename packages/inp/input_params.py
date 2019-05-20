
# import re
from . import CMD_phot_systs


# DEPRECATED 19/04/17
# def char_remove(in_lst):
#     '''
#     Correctly convert input data for parameters ranges to lists.
#     '''
#     lst = []
#     # If input list is empty, return empty list.
#     if in_lst[1:]:
#         l0 = []
#         if in_lst[1][0] in {'[', '(', '{'}:
#             # Remove non-numeric characters and append numbers as floats.
#             l0.append([
#                 float(i) for i in re.findall('[0-9.]+', str(in_lst[1:]))])
#             # Store indicating that this is a list of values.
#             lst = ['l', list(map(float, l0[0]))]
#         else:
#             if len(in_lst[1:4]) < 3:
#                 # Not enough values to define a range. Store as list of values.
#                 lst = ['l', list(map(float, in_lst[1:4]))]
#             else:
#                 # Store indicating that this is a range of values.
#                 lst = ['r', list(map(float, in_lst[1:4]))]

#     return lst


def main(mypath, pars_f_path):
    '''
    This function reads the input data parameters stored in the
    'params_input_XX.dat' file and returns a dictionary.
    '''

    # Accept these variations of the 'true' flag.
    true_lst = ('y', 'Y', 'yes', 'Yes', 'YES')

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
                    run_mode = str(reader[1])
                elif reader[0] == 'NV':
                    nanvals = [_.replace(',', '') for _ in reader[1:]]

                # Input data parameters.
                elif reader[0] == 'MR':
                    read_mode = str(reader[1])
                elif reader[0] == 'PI':
                    id_coords = reader[1:]
                elif reader[0] == 'PM':
                    id_mags = reader[1:]
                elif reader[0] == 'PC':
                    id_cols = reader[1:]
                elif reader[0] == 'PK':
                    id_kinem = reader[1:]

                # Output parameters.
                elif reader[0] == 'MP':
                    flag_make_plot = reader[1:]
                elif reader[0] == 'PF':
                    plot_frmt = str(reader[1])
                    plot_dpi = int(reader[2])
                elif reader[0] == 'MD':
                    flag_move_file = True if reader[1] in true_lst else False

                # Structure functions parameters.
                elif reader[0] == 'CH':
                    center_bw = float(reader[1])
                elif reader[0] == 'CR':
                    radius_method = str(reader[1])
                    fdens_method = str(reader[2])
                elif reader[0] == 'KP':
                    kp_flag = True if reader[1] in true_lst else False
                elif reader[0] == 'GR':
                    try:
                        fr_number = int(reader[1])
                    except ValueError:
                        fr_number = str(reader[1])

                # Data analysis functions parameters.
                elif reader[0] == 'ER':
                    err_max = reader[1:]
                elif reader[0] == 'AD':
                    ad_runs = int(reader[1])
                    ad_k_comb = True if reader[2] in true_lst else False
                elif reader[0] == 'PP':
                    plx_flag = True if reader[1] in true_lst else False
                    plx_chains = int(reader[2])
                    plx_runs = int(reader[3])
                    pms_flag = True if reader[4] in true_lst else False
                    pms_chains = int(reader[5])
                    pms_runs = int(reader[6])

                # Decontamination algorithm parameters
                elif reader[0] == 'DA':
                    da_algor = reader[1]
                    bayesda_runs = int(reader[2])
                    fixedda_port = float(reader[3])
                    readda_idcol = int(reader[4])
                    readda_mpcol = int(reader[5])

                elif reader[0] == 'DW':
                    bayesda_weights = list(map(float, reader[1:]))

                # Cluster region field stars removal.
                elif reader[0] == 'FR':
                    fld_clean_mode = str(reader[1])
                    fld_clean_bin = str(reader[2])
                    fld_clean_prob = float(reader[3])

                # Cluster parameters assignation.
                elif reader[0] == 'CF':
                    best_fit_algor = str(reader[1])
                    # TODO extend this param to 'brute force'
                    hmax = float(reader[2])

                # Ranges for the fundamental parameters
                elif reader[0] == 'RZ':
                    z_range = list(map(float, reader[1:]))
                elif reader[0] == 'RA':
                    a_range = list(map(float, reader[1:]))
                elif reader[0] == 'RE':
                    e_range = list(map(float, reader[1:]))
                elif reader[0] == 'RD':
                    d_range = list(map(float, reader[1:]))
                elif reader[0] == 'RM':
                    m_range = list(map(float, reader[1:]))
                elif reader[0] == 'RB':
                    b_range = list(map(float, reader[1:]))

                # ptemcee algorithm parameters.
                elif reader[0] == 'PT0':
                    init_mode_ptm = reader[1]
                    popsize_ptm = int(float(reader[2]))
                    maxiter_ptm = int(float(reader[3]))
                elif reader[0] == 'PT1':
                    ntemps = reader[1]
                    nwalkers_ptm = int(float(reader[2]))
                    nburn_ptm = int(float(reader[3]))
                    nsteps_ptm = int(float(reader[4]))
                    tmax_ptm = reader[5]
                    pt_adapt = True if reader[6] in true_lst else False
                elif reader[0] == 'PT2':
                    N_conv = float(reader[1])
                    tol_conv = float(reader[2])
                elif reader[0] == 'PTZ':
                    pt_z_prior = [reader[1]] + list(map(float, reader[2:]))
                elif reader[0] == 'PTA':
                    pt_a_prior = [reader[1]] + list(map(float, reader[2:]))
                elif reader[0] == 'PTE':
                    pt_e_prior = [reader[1]] + list(map(float, reader[2:]))
                elif reader[0] == 'PTD':
                    pt_d_prior = [reader[1]] + list(map(float, reader[2:]))
                elif reader[0] == 'PTM':
                    pt_m_prior = [reader[1]] + list(map(float, reader[2:]))
                elif reader[0] == 'PTB':
                    pt_b_prior = [reader[1]] + list(map(float, reader[2:]))

                # TODO not finished yet
                # # ABC algorithm parameters.
                # elif reader[0] == 'AB':
                #     nwalkers_abc = int(float(reader[1]))
                #     nburn_abc = float(reader[2])
                #     nsteps_abc = int(float(reader[3]))
                #     priors_abc = reader[4]

                # TODO not finished yet
                # # emcee algorithm parameters.
                # elif reader[0] == 'EM':
                #     nwalkers_emc = int(float(reader[1]))
                #     nburn_emc = int(float(reader[2]))
                #     N_burn_emc = int(float(reader[3]))
                #     nsteps_emc = int(float(reader[4]))
                #     emcee_a = float(reader[5])
                #     priors_emc = reader[6]

                # Bootstrap parameters.
                elif reader[0] == 'BT':
                    hperc_btstrp = float(reader[1])
                    N_pop_btstrp = int(reader[2])
                    N_step_btstrp = int(reader[3])

                # Genetic algorithm parameters.
                elif reader[0] == 'GA':
                    N_pop = int(reader[1])
                    N_gen = int(reader[2])
                    fit_diff = float(reader[3])
                    cross_prob = float(reader[4])
                    cross_sel = str(reader[5])
                    mut_prob = float(reader[6])
                    N_el = int(reader[7])
                    N_ei = int(reader[8])
                    N_es = int(reader[9])
                elif reader[0] == 'GS':
                    ga_steps = list(map(float, reader[1:]))

                # Likelihood function
                elif reader[0] == 'LK':
                    lkl_method = str(reader[1])
                    lkl_binning = str(reader[2])
                    lkl_weight = str(reader[3])

                # Synthetic clusters parameters
                elif reader[0] == 'ET':
                    evol_track = str(reader[1])
                elif reader[0] == 'PS':
                    za_steps = list(map(float, reader[1:]))
                elif reader[0] == 'MF':
                    IMF_name = str(reader[1])
                    m_high = float(reader[2])
                    m_sample_flag = True if reader[3] in true_lst else False
                    N_IMF_interp = int(reader[4])
                    m_step = float(reader[5])
                elif reader[0] == 'BI_m':
                    bin_mr = float(reader[1])
                elif reader[0] == 'RV':
                    R_V = float(reader[1])
                elif reader[0] == 'MM':
                    try:
                        max_mag = float(reader[1])
                    except ValueError:
                        max_mag = str(reader[1])

                else:
                    # Get parameters file name from path.
                    pars_f_name = pars_f_path.split('/')[-1]
                    print("  WARNING: Unknown '{}' ID found in line {}\n"
                          "  of '{}' file.\n".format(
                              reader[0], l + 1, pars_f_name))

    # Accepted coordinate units
    coord_accpt = ('px', 'deg')
    # Accepted read modes
    read_mode_accpt = ('nam', 'num')
    # Accepted decontamination algorithms.
    da_algors_accpt = ('bayes', 'fixed', 'read', 'skip')
    # Accepted field stars removal methods.
    fld_rem_methods = ('local', 'n_memb', 'mp_05', 'top_h', 'man', 'all')
    # Accepted binning methods.
    bin_methods = ('fixed', 'auto', 'fd', 'doane', 'scott', 'rice', 'sqrt',
                   'sturges', 'knuth', 'blocks', 'blocks_max')
    bin_weights = ('mean', 'median', 'max')
    # Likelihood methods.
    lkl_methods = (
        'tolstoy', 'duong', 'dolphin', 'mighell', 'kdeKL', 'dolphin_kde')
    # Accepted IMF functions.
    imf_funcs = ('chabrier_2001_exp', 'chabrier_2001_log', 'kroupa_1993',
                 'kroupa_2002')
    # Optimizing algorithm
    # TODO 'brute', 'emcee', 'abc'
    optimz_algors = ('ptemcee', 'boot+GA', 'n')
    # Accepted forms of priors.
    bayes_priors = ('u', 'g')

    priors_ptm = [
        pt_z_prior, pt_a_prior, pt_e_prior, pt_d_prior, pt_m_prior, pt_b_prior]

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

    # Dictionary with data on the CMD service photometric systems.
    cmd_systs = CMD_phot_systs.main()

    par_ranges = [z_range, a_range, e_range, d_range, m_range, b_range]

    pd = {
        'up_flag': up_flag, 'run_mode': run_mode, 'nanvals': nanvals,
        'read_mode': read_mode,
        'id_coords': id_coords, 'id_mags': id_mags, 'id_cols': id_cols,
        'id_kinem': id_kinem,
        'flag_make_plot': flag_make_plot, 'plot_frmt': plot_frmt,
        'plot_dpi': plot_dpi,
        'flag_move_file': flag_move_file,
        'center_bw': center_bw, 'radius_method': radius_method,
        'fdens_method': fdens_method,
        'kp_flag': kp_flag, 'fr_number': fr_number, 'err_max': err_max,
        'ad_runs': ad_runs, 'ad_k_comb': ad_k_comb, 'plx_flag': plx_flag,
        'plx_chains': plx_chains, 'plx_runs': plx_runs, 'pms_flag': pms_flag,
        'pms_chains': pms_chains, 'pms_runs': pms_runs,
        # Decontamination algorithm parameters.
        'da_algor': da_algor, 'bayesda_runs': bayesda_runs,
        'fixedda_port': fixedda_port, 'readda_idcol': readda_idcol,
        'readda_mpcol': readda_mpcol, 'bayesda_weights': bayesda_weights,
        # Cluster region field stars removal parameters.
        'fld_clean_mode': fld_clean_mode, 'fld_clean_bin': fld_clean_bin,
        'fld_clean_prob': fld_clean_prob,
        # Best fit parameters.
        'best_fit_algor': best_fit_algor, 'hmax': hmax, 'N_conv': N_conv,
        'tol_conv': tol_conv, 'lkl_method': lkl_method,
        'lkl_binning': lkl_binning, 'lkl_weight': lkl_weight,
        # Synthetic cluster parameters
        'evol_track': evol_track, 'za_steps': za_steps,
        'max_mag': max_mag, 'IMF_name': IMF_name, 'm_high': m_high,
        'm_sample_flag': m_sample_flag, 'N_IMF_interp': N_IMF_interp,
        'm_step': m_step, 'R_V': R_V, 'bin_mr': bin_mr,
        # parameters ranges
        'par_ranges': par_ranges,
        # ptemcee algorithm parameters.
        'init_mode_ptm': init_mode_ptm, 'popsize_ptm': popsize_ptm,
        'maxiter_ptm': maxiter_ptm,
        'ntemps': ntemps, 'nwalkers_ptm': nwalkers_ptm, 'nburn_ptm': nburn_ptm,
        'nsteps_ptm': nsteps_ptm, "pt_adapt": pt_adapt, 'tmax_ptm': tmax_ptm,
        'priors_ptm': priors_ptm,
        # # ABC algorithm parameters.
        # 'nwalkers_abc': nwalkers_abc, 'nburn_abc': nburn_abc,
        # 'nsteps_abc': nsteps_abc, 'priors_abc': priors_abc,
        # # emcee algorithm parameters.
        # 'nwalkers_emc': nwalkers_emc, 'nburn_emc': nburn_emc,
        # "N_burn_emc": N_burn_emc, 'nsteps_emc': nsteps_emc,
        # "emcee_a": emcee_a, 'priors_emc': priors_emc,
        # Bootstrap parameters
        'hperc_btstrp': hperc_btstrp, 'N_pop_btstrp': N_pop_btstrp,
        'N_step_btstrp': N_step_btstrp,
        # Genetic algorithm parameters.
        'N_pop': N_pop, 'N_gen': N_gen, 'fit_diff': fit_diff,
        'cross_prob': cross_prob, 'cross_sel': cross_sel, 'mut_prob': mut_prob,
        'N_el': N_el, 'N_ei': N_ei, 'N_es': N_es, 'ga_steps': ga_steps,
        # Fixed accepted parameter values and photometric systems.
        'read_mode_accpt': read_mode_accpt, 'coord_accpt': coord_accpt,
        'da_algors_accpt': da_algors_accpt, 'fld_rem_methods': fld_rem_methods,
        'bin_methods': bin_methods, 'bin_weights': bin_weights,
        'imf_funcs': imf_funcs, 'lkl_methods': lkl_methods,
        'optimz_algors': optimz_algors, 'bayes_priors': bayes_priors,
        'cmd_evol_tracks': cmd_evol_tracks, 'cmd_systs': cmd_systs}

    # Return parameters dictionary.
    return pd
