
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
import add_version_plot
import mp_best_fit1
import prep_plots


def main(
        npd, pd, isoch_fit_params, fit_params_r, fit_errors_r, **kwargs):
    '''
    Make D1 block plots.
    '''
    # flag_make_plot = pd['pl_params'][0]
    if pd['pl_params'][0] and pd['bf_flag']:
        best_fit_algor, lkl_method, bin_method, N_b = pd['bf_params']

        fig = plt.figure(figsize=(30, 25))
        gs = gridspec.GridSpec(10, 12)
        if best_fit_algor == 'genet':
            add_version_plot.main()
        else:
            add_version_plot.main(y_fix=.8)

        min_max_p = prep_plots.param_ranges(pd['fundam_params'])
        # Get special axis ticks for metallicity.
        xp_min, xp_max = min_max_p[0]
        # The max number of characters in the axis '30', is HARD-CODED.
        # Add values to the end of this list.
        min_max_p.append(prep_plots.BestTick(xp_min, xp_max, 30))

        # Best fitting process plots for GA.
        if best_fit_algor == 'genet':

            lkl_old, new_bs_indx, model_done = isoch_fit_params[1:]
            l_min_max = prep_plots.likl_y_range(lkl_old)

            arglist = [
                # pl_ga_lkl: Likelihood evolution for the GA.
                [gs, l_min_max, lkl_old, model_done, new_bs_indx,
                 pd['ga_params'], N_b],
                # pl_2_param_dens: Param vs param solutions scatter map.
                [gs, 'age-metal', min_max_p, fit_params_r, fit_errors_r,
                 model_done],
                [gs, 'dist-ext', min_max_p, fit_params_r, fit_errors_r,
                 model_done],
                [gs, 'metal-dist', min_max_p, fit_params_r, fit_errors_r,
                 model_done],
                [gs, 'mass-binar', min_max_p, fit_params_r, fit_errors_r,
                 model_done]
            ]
            for n, args in enumerate(arglist):
                mp_best_fit1.plot(n, *args)

        arglist = [
            # pl_lkl_scatt: Parameter likelihood density plot.
            [gs, '$z$', min_max_p, fit_params_r, fit_errors_r, model_done],
            [gs, '$log(age)$', min_max_p, fit_params_r, fit_errors_r,
             model_done],
            [gs, '$E_{{(B-V)}}$', min_max_p, fit_params_r, fit_errors_r,
             model_done],
            [gs, '$(m-M)_o$', min_max_p, fit_params_r, fit_errors_r,
             model_done],
            [gs, '$M\,(M_{{\odot}})$', min_max_p, fit_params_r, fit_errors_r,
             model_done],
            [gs, '$b_{{frac}}$', min_max_p, fit_params_r, fit_errors_r,
             model_done]
        ]
        for n, args in enumerate(arglist):
            mp_best_fit1.plot(n + 5, *args)

        # Generate output file.
        fig.tight_layout()
        pl_fmt, pl_dpi = pd['pl_params'][1:3]
        plt.savefig(
            join(npd['output_subdir'], str(npd['clust_name']) +
                 '_D1.' + pl_fmt), dpi=pl_dpi, bbox_inches='tight')
        # Close to release memory.
        plt.clf()
        plt.close()

        print("<<Plots for D1 block created>>")
