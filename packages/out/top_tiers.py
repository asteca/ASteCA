
from os.path import join
from time import strftime
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from .._version import __version__
from ..best_fit.best_fit_synth_cl import synth_cl_plot
from ..best_fit import imf
import top_tiers_plot
import prep_plots


def top_tiers_file(output_subdir, clust_name, best_model, top_tiers):
    '''
    Create output file containing top tiers.
    '''

    # Output file name.
    out_file_name = join(output_subdir, clust_name + '_top_tiers.dat')
    now_time = strftime("%Y-%m-%d %H:%M:%S")

    out_data_file = open(out_file_name, 'w')
    out_data_file.write("#\n\
# [ASteCA {}]\n\
#\n\
# Created: {}\n\
#\n\
# The 'Top tiers' models present a difference of <5% with the best likelihood\n\
# value found for the 'Best model' and one or more parameters differing\n\
# >10% with the best value assigned for said parameter(s).\n\
#\n\n\
# Best model\n\
#\n\
# Lkl               z   log(age)    d (kpc)     E(B-V)     M (Mo)       b_fr\n\
{:<10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}\n\n\
# Top tiers\n\
#\n\
# Lkl               z   log(age)    d (kpc)     E(B-V)     M (Mo)       b_fr\n"
    .format(__version__, now_time, *best_model))
    out_data_file.close()

    # Write values to file.
    with open(out_file_name, "a") as f_out:
        for mod in top_tiers:
            f_out.write('''{:<10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}'''
                        .format(*mod))
            f_out.write('\n')


def plot_top_tiers(top_tiers_flo, output_subdir, clust_name, mags,
                   cols, ip_list, err_lst, completeness, pd):
    '''
    Plot all top tiers.
    '''
    e_max, bin_mr, cmd_sel, pl_params = pd['er_params'][1],\
        pd['sc_params'][2], pd['ps_params'][1], pd['pl_params']
    # Obtain mass distribution using the selected IMF.
    st_dist_mass = imf.main(pd['sc_params'])

    # figsize(x1, y1), GridSpec(y2, x2) --> To have square plots: x1/x2 =
    # y1/y2 = 2.5
    fig = plt.figure(figsize=(25, 10))  # create the top-level container
    gs = gridspec.GridSpec(4, 10)      # create a GridSpec object
    # Add version number to top left.
    ver = '[ASteCA ' + __version__ + ']'
    x_coord = 0.95 - (len(__version__) - 6) * 0.001
    plt.figtext(x_coord, .986, ver, fontsize=9, color='#585858')

    phot_x, phot_y = prep_plots.ax_data(mags, cols)
    x_ax, y_ax, x_ax0, y_axis = prep_plots.ax_names(pd['axes_params'])
    x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd = prep_plots.diag_limits(
        y_axis, phot_x, phot_y)

    synth_cls = []
    for mod in top_tiers_flo:
        # Call function to generate synthetic cluster.
        shift_isoch, synth_clst = synth_cl_plot(
            ip_list, [mod], err_lst, completeness, st_dist_mass,
            e_max, bin_mr, cmd_sel)
        # Store plotting parameters.
        synth_cls.append([gs, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd,
                         x_ax, y_ax, synth_clst, mod, shift_isoch])

    # Generate all plots.
    for n, args in enumerate(synth_cls, 1):
        top_tiers_plot.plot(n, *args)

    fig.tight_layout()
    # Generate output file for each data file.
    pl_fmt, pl_dpi = pl_params[1:3]
    plt.savefig(join(output_subdir, clust_name + '_top_tiers.' + pl_fmt),
                dpi=pl_dpi)
    # Close to release memory.
    plt.clf()
    plt.close()


def main(npd, cld, pd, err_lst, completeness, isoch_fit_params, **kwargs):
    '''
    Obtain top tier models, produce data file and output image.
    '''

    clust_name, output_subdir = npd['clust_name'], npd['output_subdir']
    ip_list, bf_flag = pd['ip_list'], pd['bf_params'][0]
    mags, cols = cld['mags'], cld['cols']

    if bf_flag:
        all_models = isoch_fit_params[-1]

        # Sort all models/solutions by their (minimum) likelihood values.
        all_m_sort = sorted(zip(*all_models), key=lambda x: x[1])
        best_mod, best_lik = all_m_sort[0][0], all_m_sort[0][1]
        best_lik_5 = best_lik * 0.05

        top_tiers_str, top_tiers_flo = [], []
        # Compare every model/solution with the best one.
        for mod in all_m_sort[1:]:
            # If the likelihood is closer than 5% to the min likelihood.
            if (mod[1] - best_lik) < best_lik_5:
                # If any parameter has a difference larger than 10%.
                mod_flag, p_mark = False, []
                for i, p in enumerate(mod[0]):
                    if abs(p - best_mod[i]) > best_mod[i] * 0.1:
                        mod_flag = True
                        p_mark.append(str(p) + '*')
                    else:
                        p_mark.append(str(p))
                if mod_flag:
                    top_tiers_str.append(['{:.3f}'.format(mod[1])] + p_mark)
                    top_tiers_flo.append(mod[0])
            # Exit when 10 models are stored.
            if len(top_tiers_str) == 10:
                break

        # Check if at least one model was stored.
        if top_tiers_str:

            # Store best model values as strings.
            best_model = ['{:.3f}'.format(best_lik)] + map(str, best_mod)

            # Create output .dat file.
            top_tiers_file(output_subdir, clust_name, best_model,
                           top_tiers_str)

            # Create output image if flag is 'true'.
            flag_make_plot = pd['pl_params'][0]
            if flag_make_plot:
                try:
                    plot_top_tiers(
                        top_tiers_flo, output_subdir, clust_name, mags,
                        cols, ip_list, err_lst, completeness, pd)
                    print('Top tier models saved to file and plotted.')
                except:
                    import traceback
                    print traceback.format_exc()
                    print("  ERROR: top tiers plot could not be generated.")
            else:
                # If top tiers were found but no plot is produced.
                print("Top tier models saved to file.")
