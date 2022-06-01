
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
from . import add_version_plot
from . import mp_mcmc_cnvrg
from . import tracePlot
from . import prep_plots
from . prep_plots import figsize_x, figsize_y, grid_x, grid_y


def main(npd, pd, clp, td):
    """
    Make D1 block plots.
    """
    fig = plt.figure(figsize=(figsize_x, figsize_y))
    gs = gridspec.GridSpec(grid_y, grid_x)

    fit_pars = clp['isoch_fit_params']

    # Title and version for the plot.
    m_bf, s_bf = divmod(fit_pars['bf_elapsed'], 60)
    h_bf, m_bf = divmod(m_bf, 60)

    xf, yf = .02, .999
    add_version_plot.main(x_fix=xf, y_fix=yf)

    p_str = (
        "chains={:.0f}, burn={:.2f}, steps={:.0f},"
        " adapt={}").format(
            pd['nwalkers_mcee'], pd['nburn_mcee'], fit_pars['N_steps'][-1],
            pd['pt_adapt'])

    xt, yt = .5, 1.005
    plt.suptitle(
        ("{} | {:.0f}h{:.0f}m").format(p_str, h_bf, m_bf), x=xt, y=yt,
        fontsize=11)

    # Trace plots
    min_max_p = prep_plots.param_ranges(
        td['fundam_params'], clp['varIdxs'], fit_pars['pars_chains'])
    trace = fit_pars['mcmc_trace']
    best_sol = fit_pars['mean_sol']
    traceplot_args = (
        fit_pars['acorr_t'], fit_pars['med_at_c'], fit_pars['mcmc_ess'])
    post_trace, pre_trace = fit_pars['pars_chains'], fit_pars['pars_chains_bi']

    # pl_param_chain: Parameters sampler chains.
    par_list = ['metal', 'age', 'beta', 'ext', 'dr', 'rv', 'dist']
    for p in par_list:
        args = [
            p, gs, best_sol, min_max_p, traceplot_args,
            trace, clp['varIdxs'], post_trace, pre_trace]
        tracePlot.plot(0, *args)

    # Parallel Coordinates plot
    # mp_mcmc_cnvrg.plot(2, *args)

    # pl_MAP_lkl: Parameters half of pdfs.
    args = [
        gs, fit_pars['N_steps'], fit_pars['lkl_mean_steps'],
        fit_pars['lkl_steps']]
    mp_mcmc_cnvrg.plot(0, *args)

    # if pd['best_fit_algor'] == 'ptemcee':
    # pl_betas: Betas vs steps.
    args = [gs, fit_pars['Tmax'], fit_pars['N_steps'], fit_pars['betas_pt']]
    mp_mcmc_cnvrg.plot(2, *args)

    # pl_Tswaps: Tswaps AFs vs steps.
    args = [gs, fit_pars['N_steps'], fit_pars['tswaps_afs']]
    mp_mcmc_cnvrg.plot(3, *args)

    # pl_MAF: Parameters evolution of MAF.
    maf_steps = fit_pars['maf_allT']
    # if pd['best_fit_algor'] == 'ptemcee':
    # elif pd['best_fit_algor'] == 'emcee':
    #     maf_steps = fit_pars['maf_steps']
    args = [gs, fit_pars['N_steps'], maf_steps]
    mp_mcmc_cnvrg.plot(1, *args)

    # pl_tau
    args = [gs, fit_pars['N_steps'], fit_pars['tau_autocorr']]
    mp_mcmc_cnvrg.plot(4, *args)

    # TODO re-implement when/if code is fixed
    # # pl_mESS
    # args = [
    #     'mESS', gs, fit_pars['mESS'],
    #     fit_pars['minESS'],
    #     fit_pars['mESS_epsilon']]
    # mp_mcmc_cnvrg.plot(7, *args)

    # pl_lags
    args = [gs, clp['varIdxs'], fit_pars['acorr_function'], par_list]
    mp_mcmc_cnvrg.plot(5, *args)

    # pl_GW
    args = [gs, clp['varIdxs'], fit_pars['geweke_z'], par_list]
    mp_mcmc_cnvrg.plot(6, *args)

    # pl_tau_histo
    args = [gs, fit_pars['all_taus']]
    mp_mcmc_cnvrg.plot(7, *args)

    # Generate output file.
    fig.tight_layout()
    plt.savefig(join(
        npd['output_subdir'], str(npd['clust_name']) + '_D1_'
        + pd['best_fit_algor'] + npd['ext']))
    # Close to release memory.
    plt.clf()
    plt.close("all")
