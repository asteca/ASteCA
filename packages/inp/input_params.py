
from . import CMD_phot_systs


def main(pars_f_path):
    """
    This function reads the input data parameters stored in the
    'params_input_XX.dat' file and returns a dictionary.
    """

    # Accept these variations of the 'true' flag.
    true_lst = ('y', 'Y', 'yes', 'Yes', 'YES', 'true', 'True', 'TRUE')

    # Read data from file.
    with open(pars_f_path, "r") as f_dat:

        manual_struct, trim_frame_range = [], []
        # Iterate through each line in the file.
        for ln, line in enumerate(f_dat):

            if not line.startswith("#") and line.strip() != '':
                reader = line.split()

                # Input data parameters.
                if reader[0] == 'I0':
                    read_mode = str(reader[1])
                elif reader[0] == 'I1':
                    id_ids = reader[1]
                    id_xdata = reader[2]
                    id_ydata = reader[3]
                    coords = reader[4]
                    project = True if reader[5] in true_lst else False
                elif reader[0] == 'I2':
                    id_mags = reader[1:]
                elif reader[0] == 'I3':
                    id_cols = reader[1:]
                elif reader[0] == 'I4':
                    id_kinem = reader[1:]

                # Input data processing
                elif reader[0] == 'I5':
                    nanvals = [_.replace(',', '') for _ in reader[1:]]
                elif reader[0] == 'I6':
                    trim_frame_range.append([
                        reader[1], list(map(float, reader[2:]))])

                # Structure functions parameters.
                elif reader[0] == 'S0':
                    manual_struct.append(reader[1:])
                elif reader[0] == 'S1':
                    center_bw = float(reader[1])
                    mirror_flag = True if reader[2] in true_lst else False
                    NN_dd = int(reader[3])
                    fdens_method = str(reader[4])
                    nsteps_rad = int(reader[5])
                elif reader[0] == 'S2':
                    kp_ndim = int(reader[1])
                    kp_nchains = int(reader[2])
                    kp_nruns = int(reader[3])
                    kp_nburn = float(reader[4])
                    rt_max_f = float(reader[5])
                elif reader[0] == 'S3':
                    kp_emcee_moves = [_.strip() for _ in line[3:].split(';')]
                elif reader[0] == 'S4':
                    try:
                        fr_number = int(reader[1])
                    except ValueError:
                        fr_number = str(reader[1])

                # Data analysis functions parameters.
                elif reader[0] == 'E0':
                    err_max = reader[1:]
                elif reader[0] == 'A0':
                    ad_runs = int(reader[1])

                # Decontamination algorithm parameters
                elif reader[0] == 'D0':
                    da_algor = reader[1]
                    bayesda_runs = int(reader[2])
                    bayesda_dflag = reader[3:]

                # Cluster region field stars removal.
                elif reader[0] == 'F0':
                    fld_clean_mode = str(reader[1])
                    fld_clean_bin = str(reader[2])
                    fld_clean_prob = float(reader[3])

                # Parallax & PMS analysis
                elif reader[0] == 'P0':
                    plx_bayes_flag = True if reader[1] in true_lst else False
                    plx_offset = float(reader[2])
                    plx_chains = int(reader[3])
                    plx_runs = int(reader[4])
                    plx_burn = float(reader[5])
                    flag_plx_mp = True if reader[6] in true_lst else False
                elif reader[0] == 'P1':
                    PM_KDE_std = float(reader[1])
                    cosDE_flag = True if reader[2] in true_lst else False

                # Synthetic clusters parameters
                elif reader[0] == 'R0':
                    synth_rand_seed = str(reader[1])
                    evol_track = str(reader[2])
                    IMF_name = str(reader[3])
                    bin_mr = float(reader[4])
                    try:
                        max_mag = float(reader[5])
                    except ValueError:
                        max_mag = str(reader[5])

                # Ranges for the fundamental parameters
                elif reader[0] == 'RZ':
                    z_range = reader[1:]
                elif reader[0] == 'RA':
                    a_range = reader[1:]
                elif reader[0] == 'RE':
                    e_range = list(map(float, reader[1:]))
                elif reader[0] == 'RR':
                    dr_range = reader[1:]
                elif reader[0] == 'RV':
                    R_V = float(reader[1])
                elif reader[0] == 'RD':
                    d_range = list(map(float, reader[1:]))
                elif reader[0] == 'RM':
                    m_range = list(map(float, reader[1:]))
                elif reader[0] == 'RB':
                    b_range = list(map(float, reader[1:]))
                elif reader[0] == 'RS':
                    bs_range = reader[1:]

                # Cluster parameters assignation.
                elif reader[0] == 'B0':
                    best_fit_algor = str(reader[1])
                    mins_max = float(reader[2])
                    save_trace_flag = True if reader[3] in true_lst else False

                # Shared ptemcee parameters.
                elif reader[0] == 'B1':
                    nsteps_mcee = int(float(reader[1]))
                    nwalkers_mcee = int(float(reader[2]))
                    nburn_mcee = float(reader[3])
                    pt_ntemps = reader[4]
                    pt_tmax = reader[5]
                    pt_adapt = True if reader[6] in true_lst else False

                # Priors
                elif reader[0] == 'BZ':
                    z_prior = [reader[1]] + list(map(float, reader[2:]))
                elif reader[0] == 'BA':
                    a_prior = [reader[1]] + list(map(float, reader[2:]))
                elif reader[0] == 'BE':
                    e_prior = [reader[1]] + list(map(float, reader[2:]))
                elif reader[0] == 'BR':
                    dr_prior = [reader[1]] + list(map(float, reader[2:]))
                elif reader[0] == 'BV':
                    rv_prior = [reader[1]] + list(map(float, reader[2:]))
                elif reader[0] == 'BD':
                    d_prior = [reader[1]] + list(map(float, reader[2:]))
                elif reader[0] == 'BM':
                    m_prior = [reader[1]] + list(map(float, reader[2:]))
                elif reader[0] == 'BB':
                    b_prior = [reader[1]] + list(map(float, reader[2:]))
                elif reader[0] == 'BS':
                    bs_prior = [reader[1]] + list(map(float, reader[2:]))

                # Likelihood function
                elif reader[0] == 'B2':
                    lkl_method = str(reader[1])
                    lkl_binning = str(reader[2])
                    lkl_manual_bins = reader[3:]

                # Output parameters.
                elif reader[0] == 'O0':
                    flag_make_plot = reader[1:]
                elif reader[0] == 'O1':
                    # plot_frmt = str(reader[1])
                    # plot_dpi = int(reader[2])
                    plot_style = str(reader[1])
                    # TODO this flag is hidden for now
                    mirror_flag = True

                else:
                    # Get parameters file name from path.
                    pars_f_name = pars_f_path.split('/')[-1]
                    print("  WARNING: Unknown '{}' ID found in line {}\n"
                          "  of '{}' file.\n".format(
                              reader[0], ln + 1, pars_f_name))

    # Accepted coordinate units
    coord_accpt = ('px', 'deg')
    # Accepted read modes
    read_mode_accpt = ('nam', 'num')
    # Accepted decontamination algorithms.
    da_algors_accpt = ('bayes', 'read', 'skip')
    # Accepted field stars removal methods.
    fld_rem_methods = ('local', 'n_memb', 'mp_05', 'top_h', 'man', 'all')
    # Accepted binning methods.
    bin_methods = (
        'optm', 'fixed', 'auto', 'fd', 'doane', 'scott', 'rice', 'sqrt',
        'sturges', 'knuth', 'blocks', 'blocks-max', 'manual')
    # Likelihood methods.
    # 'duong', 'kdeKL', 'dolphin_kde'
    lkl_methods = (
        'tremmel', 'dolphin', 'mighell', 'tolstoy', 'isochfit')
    # Accepted IMF functions.
    imf_funcs = ('chabrier_2001_exp', 'chabrier_2001_log', 'kroupa_1993',
                 'kroupa_2002', 'salpeter_1955')
    # Optimizing algorithm
    optimz_algors = ('ptemcee', 'read', 'n', 'synth_gen')
    # Accepted forms of priors.
    bayes_priors = ('u', 'g')
    plots_names = ('A1', 'A2', 'A3', 'B1', 'B2', 'B3', 'C1', 'C2', 'C3', 'D0',
                   'D1', 'D2', 'D3', 's')

    priors_mcee = [z_prior, a_prior, e_prior, d_prior, m_prior, b_prior]

    # Map the selected evolutionary tracks to their proper folder names,
    # as given by the ezPADOVA-2 package.
    all_evol_tracks = {
        'PAR12+CS_37': 'parsec12_37', 'PAR12+CS_35': 'parsec12_35',
        'PAR12+CS_07': 'parsec12_07', 'PAR12+CPR16': 'parsec12_16',
        'PAR12+No': 'parsec12_No',
    }

    # HARDCODED AND IMPORTANT
    # If the CMD isochrones change, this needs to change too.
    # Names of the "extra" columns in the CMD service isochrones.
    CMD_extra_pars = (
        'Mini', 'int_IMF', 'Mass', 'logL', 'logTe', 'logg', 'label',
        'mbolmag')

    # Dictionary with data on the CMD service photometric systems.
    cmd_systs = CMD_phot_systs.main()

    par_ranges = [z_range, a_range, e_range, d_range, m_range, b_range]

    pd = {
        # Input data parameters
        'read_mode': read_mode, 'id_ids': id_ids, 'id_xdata': id_xdata,
        'id_ydata': id_ydata, 'coords': coords, 'project': project,
        'id_mags': id_mags, 'id_cols': id_cols, 'id_kinem': id_kinem,

        # Input data processing
        'nanvals': nanvals, 'trim_frame_range': trim_frame_range,

        # Structure functions parameters
        'manual_struct': manual_struct, 'center_bw': center_bw,
        'mirror_flag': mirror_flag, 'NN_dd': NN_dd,
        'fdens_method': fdens_method, 'nsteps_rad': nsteps_rad,
        'kp_ndim': kp_ndim, 'kp_nchains': kp_nchains, 'kp_nruns': kp_nruns,
        'kp_nburn': kp_nburn, 'rt_max_f': rt_max_f, 'fr_number': fr_number,
        'kp_emcee_moves': kp_emcee_moves,

        #
        'err_max': err_max, 'ad_runs': ad_runs,

        # Decontamination algorithm parameters.
        'da_algor': da_algor, 'bayesda_runs': bayesda_runs,
        'bayesda_dflag': bayesda_dflag,

        # Plx & PMs parameters.
        'plx_bayes_flag': plx_bayes_flag, 'plx_offset': plx_offset,
        'plx_chains': plx_chains, 'plx_runs': plx_runs, 'plx_burn': plx_burn,
        'flag_plx_mp': flag_plx_mp,
        'PM_KDE_std': PM_KDE_std, 'cosDE_flag': cosDE_flag,

        # Cluster region field stars removal parameters.
        'fld_clean_mode': fld_clean_mode, 'fld_clean_bin': fld_clean_bin,
        'fld_clean_prob': fld_clean_prob,

        # Synthetic cluster parameters
        'synth_rand_seed': synth_rand_seed, 'par_ranges': par_ranges,
        'evol_track': evol_track, 'IMF_name': IMF_name, 'bin_mr': bin_mr,
        'R_V': R_V, 'max_mag': max_mag,

        # Best fit parameters.
        'best_fit_algor': best_fit_algor, 'mins_max': mins_max,
        'save_trace_flag': save_trace_flag,
        # ptemcee algorithm parameters.
        'nsteps_mcee': nsteps_mcee, 'nwalkers_mcee': nwalkers_mcee,
        'nburn_mcee': nburn_mcee, 'priors_mcee': priors_mcee,
        'pt_ntemps': pt_ntemps, "pt_adapt": pt_adapt, 'pt_tmax': pt_tmax,
        'lkl_method': lkl_method, 'lkl_binning': lkl_binning,
        'lkl_manual_bins': lkl_manual_bins,

        # Fixed accepted parameter values and photometric systems.
        'read_mode_accpt': read_mode_accpt, 'coord_accpt': coord_accpt,
        'da_algors_accpt': da_algors_accpt, 'fld_rem_methods': fld_rem_methods,
        'bin_methods': bin_methods,
        'imf_funcs': imf_funcs, 'lkl_methods': lkl_methods,
        'optimz_algors': optimz_algors, 'bayes_priors': bayes_priors,
        'all_evol_tracks': all_evol_tracks, 'CMD_extra_pars': CMD_extra_pars,
        'cmd_systs': cmd_systs,

        # Output
        'flag_make_plot': flag_make_plot, 'plot_style': plot_style,
        'plots_names': plots_names
    }

    return pd
