
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
from . import add_version_plot
from . import mp_best_fit1_GA
from . import mp_best_fit1_mcmc
from . import prep_plots


def main(npd, pd, isoch_fit_params, fit_params_r, fit_errors_r, **kwargs):
    '''
    Make D1 block plots.
    '''
    if 'D1' in pd['flag_make_plot'] and pd['bf_flag']:
        fig = plt.figure(figsize=(30, 30))
        gs = gridspec.GridSpec(12, 12)
        add_version_plot.main(y_fix=.999)

        # TODO DEPRECATED (DELETE)
        # # Get special axis ticks for metallicity.
        # xp_min, xp_max = min_max_p[0]
        # # The max number of characters in the axis '30', is HARD-CODED.
        # # Add values to the end of this list.
        # min_max_p.append(prep_plots.BestTick(xp_min, xp_max, 30))

        # Best fitting process plots.
        if pd['best_fit_algor'] == 'genet':
            min_max_p = prep_plots.param_ranges(
                pd['best_fit_algor'], pd['fundam_params'])
            l_min_max = prep_plots.likl_y_range(
                pd['best_fit_algor'], isoch_fit_params['lkl_old'])

            args = [
                # pl_lkl: Likelihood evolution for the GA.
                gs, l_min_max, isoch_fit_params['lkl_old'],
                isoch_fit_params['model_done'],
                isoch_fit_params['new_bs_indx'], pd['N_pop'], pd['N_gen'],
                pd['fit_diff'], pd['cross_prob'], pd['cross_sel'],
                pd['mut_prob'], pd['N_el'], pd['N_ei'], pd['N_es'],
                pd['N_bootstrap']
            ]
            mp_best_fit1_GA.plot(0, *args)

        if pd['best_fit_algor'] in ('brute', 'genet'):
            arglist = [
                # pl_2_param_dens: Param vs param solutions scatter map.
                [gs, 'age-metal', min_max_p, fit_params_r, fit_errors_r,
                 isoch_fit_params['model_done']],
                [gs, 'dist-ext', min_max_p, fit_params_r, fit_errors_r,
                 isoch_fit_params['model_done']],
                [gs, 'ext-age', min_max_p, fit_params_r, fit_errors_r,
                 isoch_fit_params['model_done']],
                [gs, 'mass-binar', min_max_p, fit_params_r, fit_errors_r,
                 isoch_fit_params['model_done']],
                # pl_lkl_scatt: Parameter likelihood density plot.
                [gs, '$z$', min_max_p, fit_params_r, fit_errors_r,
                 isoch_fit_params['model_done']],
                [gs, '$log(age)$', min_max_p, fit_params_r, fit_errors_r,
                 isoch_fit_params['model_done']],
                [gs, '$E_{{(B-V)}}$', min_max_p, fit_params_r, fit_errors_r,
                 isoch_fit_params['model_done']],
                [gs, '$(m-M)_o$', min_max_p, fit_params_r, fit_errors_r,
                 isoch_fit_params['model_done']],
                [gs, '$M\,(M_{{\odot}})$', min_max_p, fit_params_r,
                 fit_errors_r, isoch_fit_params['model_done']],
                [gs, '$b_{{frac}}$', min_max_p, fit_params_r, fit_errors_r,
                 isoch_fit_params['model_done']]
            ]
            for n, args in enumerate(arglist, 1):
                mp_best_fit1_GA.plot(n, *args)

        if pd['best_fit_algor'] == 'emcee':
            nwalkers, nburn, nsteps, mcmc_param = pd['nwalkers_emc'],\
                pd['nburn_emc'], isoch_fit_params['nsteps_emc'],\
                pd['emcee_a']
        elif pd['best_fit_algor'] == 'ptemcee':
            nwalkers, nburn, nsteps, mcmc_param = pd['nwalkers_ptm'],\
                pd['nburn_ptm'], isoch_fit_params['nsteps_ptm'],\
                pd['pt_adapt']
        elif pd['best_fit_algor'] == 'abc':
            nwalkers, nburn, nsteps, mcmc_param = pd['nwalkers_abc'],\
                isoch_fit_params['nburn_abc'],\
                isoch_fit_params['nsteps_abc'], None

        if pd['best_fit_algor'] in ('abc', 'ptemcee', 'emcee'):
            min_max_p = prep_plots.param_ranges(
                pd['best_fit_algor'], pd['fundam_params'],
                isoch_fit_params['varIdxs'], isoch_fit_params['pars_chains'])

            # pl_2_param_dens: Param vs param density map.
            for p2 in [
                    'metal-age', 'metal-ext', 'metal-dist', 'metal-mass',
                    'metal-binar', 'age-ext', 'age-dist', 'age-mass',
                    'age-binar', 'ext-dist', 'ext-mass', 'ext-binar',
                    'dist-mass', 'dist-binar', 'mass-binar']:

                # Limits for the 2-dens plots.
                min_max_p2 = prep_plots.p2_ranges(p2, min_max_p)

                args = [p2, gs, min_max_p2, isoch_fit_params['varIdxs'],
                        isoch_fit_params['mcmc_trace']]
                mp_best_fit1_mcmc.plot(0, *args)

            # pl_param_pf: Parameters probability functions.
            for p in ['metal', 'age', 'ext', 'dist', 'mass', 'binar']:
                args = [p, gs, min_max_p, fit_params_r, fit_errors_r,
                        isoch_fit_params['varIdxs'],
                        isoch_fit_params['map_sol'],
                        isoch_fit_params['mcmc_trace']]
                mp_best_fit1_mcmc.plot(1, *args)

            # Parallel Coordinates plot
            # http://benalexkeen.com/parallel-coordinates-in-matplotlib/
            # https://stackoverflow.com/questions/8230638/parallel-coordinates-plot-in-matplotlib
            # mp_best_fit1_mcmc.plot(2, *args)

            # pl_MAP_lkl: Parameters half of pdfs.
            args = [
                'MAP lkl', gs, isoch_fit_params['prob_mean'],
                isoch_fit_params['map_lkl'], isoch_fit_params['map_lkl_final']]
            mp_best_fit1_mcmc.plot(3, *args)

            # pl_MAF: Parameters evolution of MAF.
            args = [
                'MAF', gs, pd['best_fit_algor'], isoch_fit_params['maf_steps']]
            mp_best_fit1_mcmc.plot(4, *args)

            # pl_param_chain: Parameters sampler chains.
            for p in ['metal', 'age', 'ext', 'dist', 'mass', 'binar']:
                args = [
                    p, gs, pd['best_fit_algor'], fit_params_r, min_max_p,
                    nwalkers, nburn, nsteps, mcmc_param,
                    isoch_fit_params['mcmc_trace'],
                    isoch_fit_params['mcmc_elapsed'],
                    isoch_fit_params['varIdxs'],
                    isoch_fit_params['pars_chains_bi'],
                    isoch_fit_params['pars_chains'],
                    isoch_fit_params['autocorr_time'],
                    isoch_fit_params['max_at_5c'],
                    isoch_fit_params['min_at_5c'],
                    isoch_fit_params['pymc3_ess']]
                mp_best_fit1_mcmc.plot(5, *args)

            # pl_tau
            args = [
                'Tau', gs, isoch_fit_params['N_steps_conv'],
                isoch_fit_params['N_conv'], isoch_fit_params['tol_conv'],
                isoch_fit_params['tau_index'],
                isoch_fit_params['tau_autocorr']]
            mp_best_fit1_mcmc.plot(6, *args)
            # pl_mESS
            args = [
                'mESS', gs, isoch_fit_params['mESS'],
                isoch_fit_params['minESS'],
                isoch_fit_params['mESS_epsilon']]
            mp_best_fit1_mcmc.plot(7, *args)
            # pl_lags
            args = [
                'lags', gs, isoch_fit_params['varIdxs'],
                isoch_fit_params['emcee_acorf']]
            mp_best_fit1_mcmc.plot(8, *args)
            # pl_GW
            args = [
                'Geweke', gs, isoch_fit_params['varIdxs'],
                isoch_fit_params['geweke_z']]
            mp_best_fit1_mcmc.plot(9, *args)

        # if pd['best_fit_algor'] == 'abc':
        #     min_max_p = prep_plots.param_ranges(
        #         pd['best_fit_algor'], pd['fundam_params'],
        #         isoch_fit_params['varIdxs'], isoch_fit_params['pars_chains'])

        #     # pl_2_param_dens: Param vs param density map.
        #     for p2 in [
        #             'metal-age', 'metal-ext', 'metal-dist', 'metal-mass',
        #             'metal-binar', 'age-ext', 'age-dist', 'age-mass',
        #             'age-binar', 'ext-dist', 'ext-mass', 'ext-binar',
        #             'dist-mass', 'dist-binar', 'mass-binar']:

        #         # Limits for the 2-dens plots.
        #         min_max_p2 = prep_plots.p2_ranges(p2, min_max_p)

        #         args = [p2, gs, min_max_p2, isoch_fit_params['varIdxs'],
        #                 isoch_fit_params['mcmc_trace']]
        #         mp_best_fit1_mcmc.plot(0, *args)

        #     # pl_param_pf: Parameters probability functions.
        #     for p in ['metal', 'age', 'ext', 'dist', 'mass', 'binar']:
        #         args = [p, gs, min_max_p, fit_params_r, fit_errors_r,
        #                 isoch_fit_params['varIdxs'],
        #                 isoch_fit_params['map_sol'],
        #                 isoch_fit_params['mcmc_trace']]
        #         mp_best_fit1_mcmc.plot(1, *args)

        #     # # pl_pdf_half: Parameters half of pdfs.
        #     # args = ['Halves', gs, isoch_fit_params['mcmc_halves']]
        #     # mp_best_fit1_mcmc.plot(2, *args)

        #     # pl_MAP_lkl: Parameters half of pdfs.
        #     args = ['MAP lkl', gs, isoch_fit_params['map_lkl'],
        #             isoch_fit_params['map_lkl_final']]
        #     mp_best_fit1_mcmc.plot(3, *args)

        #     # pl_MAF: Parameters evolution of MAF.
        #     args = [
        #         'MAF', gs, pd['best_fit_algor'],
        #         isoch_fit_params['maf_steps']]
        #     mp_best_fit1_mcmc.plot(4, *args)

        #     # pl_param_chain: Parameters sampler chains.
        #     for p in ['metal', 'age', 'ext', 'dist', 'mass', 'binar']:
        #         args = [
        #             p, gs, pd['best_fit_algor'], fit_params_r, min_max_p,
        #             nwalkers, nburn, nsteps, aparam,
        #             isoch_fit_params['mcmc_trace'],
        #             isoch_fit_params['mcmc_elapsed'],
        #             isoch_fit_params['varIdxs'],
        #             isoch_fit_params['pars_chains_bi'],
        #             isoch_fit_params['pars_chains'],
        #             isoch_fit_params['autocorr_time'],
        #             isoch_fit_params['max_at_5c'],
        #             isoch_fit_params['min_at_5c'],
        #             isoch_fit_params['pymc3_ess']]
        #         mp_best_fit1_mcmc.plot(5, *args)

        #     # pl_tau
        #     args = [
        #         'Tau', gs, isoch_fit_params['N_steps_conv'],
        #         isoch_fit_params['N_conv'], isoch_fit_params['tol_conv'],
        #         isoch_fit_params['tau_index'],
        #         isoch_fit_params['tau_autocorr']]
        #     mp_best_fit1_mcmc.plot(6, *args)
        #     # pl_mESS
        #     args = [
        #         'mESS', gs, isoch_fit_params['mESS'],
        #         isoch_fit_params['minESS'],
        #         isoch_fit_params['mESS_epsilon']]
        #     mp_best_fit1_mcmc.plot(7, *args)
        #     # pl_lags
        #     args = [
        #         'lags', gs, isoch_fit_params['varIdxs'],
        #         isoch_fit_params['emcee_acorf']]
        #     mp_best_fit1_mcmc.plot(8, *args)
        #     # pl_GW
        #     args = [
        #         'Geweke', gs, isoch_fit_params['varIdxs'],
        #         isoch_fit_params['geweke_z']]
        #     mp_best_fit1_mcmc.plot(9, *args)

        # Generate output file.
        fig.tight_layout()
        plt.savefig(join(
            npd['output_subdir'], str(npd['clust_name']) + '_D1_' +
            pd['best_fit_algor'] + '.' + pd['plot_frmt']),
            dpi=pd['plot_dpi'], bbox_inches='tight')
        # Close to release memory.
        plt.clf()
        plt.close("all")

        print("<<Plots for D1 block created>>")
    else:
        print("<<Skip D1 block plot>>")
