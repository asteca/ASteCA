
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
from . import add_version_plot
from . import cornerPlot
from . import prep_plots
from . prep_plots import figsize_x, figsize_y, grid_x, grid_y


def main(npd, pd, isoch_fit_params, isoch_fit_errors, **kwargs):
    """
    Make D2 block plots (corner plot).
    """
    if 'D2' in pd['flag_make_plot'] and pd['best_fit_algor'] != 'n':
        fig = plt.figure(figsize=(figsize_x, figsize_y))
        gs = gridspec.GridSpec(grid_y, grid_x)

        # Title and version for the plot.
        m_bf, s_bf = divmod(isoch_fit_params['bf_elapsed'], 60)
        h_bf, m_bf = divmod(m_bf, 60)

        xf, yf = .02, .999
        add_version_plot.main(x_fix=xf, y_fix=yf)

        msol = 'MAP'
        # pl_2_param_dens: Param vs param density map.
        for p2 in [
                'metal-age', 'metal-ext', 'metal-dist', 'metal-mass',
                'metal-binar', 'age-ext', 'age-dist', 'age-mass',
                'age-binar', 'ext-dist', 'ext-mass', 'ext-binar',
                'dist-mass', 'dist-binar', 'mass-binar']:

            # Limits for the 2-dens plots.
            min_max_p = prep_plots.param_ranges(
                pd['fundam_params'], isoch_fit_params['varIdxs'],
                isoch_fit_params['pars_chains'])
            min_max_p2 = prep_plots.p2_ranges(p2, min_max_p)

            trace = isoch_fit_params['mcmc_trace']
            args = [p2, gs, min_max_p2, isoch_fit_params['varIdxs'], trace]
            cornerPlot.plot(0, *args)

        par_list = ['metal', 'age', 'ext', 'dist', 'mass', 'binar']

        # pl_param_pf: Parameters probability functions.
        for p in par_list:
            args = [
                p, gs, min_max_p, isoch_fit_params['varIdxs'],
                isoch_fit_params['mean_sol'], isoch_fit_params['map_sol'],
                isoch_fit_params['median_sol'], isoch_fit_params['mode_sol'],
                isoch_fit_params['param_r2'], isoch_fit_params['pardist_kde'],
                trace, msol]
            cornerPlot.plot(1, *args)

        # Generate output file.
        fig.tight_layout()
        plt.savefig(join(
            npd['output_subdir'], str(npd['clust_name']) + '_D2_' +
            pd['best_fit_algor']), bbox_inches='tight')
        print("<<Plots for D2 block created>>")

        plt.clf()
        plt.close("all")
    else:
        print("<<Skip D2 block plot>>")
