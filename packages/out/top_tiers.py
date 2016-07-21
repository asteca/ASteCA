
from os.path import join
from time import strftime
from .._version import __version__
from ..inp import input_params as g
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


def plot_top_tiers(top_tiers_flo, output_subdir, clust_name, mag_data,
                   col_data, ip_list, err_lst, completeness):
    '''
    Plot all top tiers.
    '''
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    # Obtain mass distribution using the selected IMF.
    st_dist_mass = imf.main()

    # figsize(x1, y1), GridSpec(y2, x2) --> To have square plots: x1/x2 =
    # y1/y2 = 2.5
    fig = plt.figure(figsize=(25, 10))  # create the top-level container
    gs = gridspec.GridSpec(4, 10)      # create a GridSpec object
    # Add version number to top left.
    ver = '[ASteCA ' + __version__ + ']'
    x_coord = 0.95 - (len(__version__) - 6) * 0.001
    plt.figtext(x_coord, .986, ver, fontsize=9, color='#585858')

    phot_x, phot_y = prep_plots.ax_data(mag_data, col_data)
    x_ax, y_ax, x_ax0, y_axis = prep_plots.ax_names()
    x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd = prep_plots.diag_limits(
        y_axis, phot_x, phot_y)

    synth_cls = []
    for mod in top_tiers_flo:
        # Call function to generate synthetic cluster.
        shift_isoch, synth_clst = synth_cl_plot(
            ip_list, [mod], err_lst, completeness, st_dist_mass)
        # Store plotting parameters.
        synth_cls.append([gs, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd,
                         x_ax, y_ax, synth_clst, mod, shift_isoch])

    # Generate all plots.
    for n, args in enumerate(synth_cls, 1):
        top_tiers_plot.plot(n, *args)

    fig.tight_layout()
    # Generate output file for each data file.
    pl_fmt, pl_dpi = g.pl_params[1:3]
    plt.savefig(join(output_subdir, clust_name + '_top_tiers.' + pl_fmt),
                dpi=pl_dpi)
    # Close to release memory.
    plt.clf()
    plt.close()


def main(clust_name, output_subdir, mag_data, col_data, ip_list,
         err_lst, completeness, bf_return):
    '''
    Obtain top tier models, produce data file and output image.
    '''

    bf_flag = g.bf_params[0]
    if bf_flag:

        isoch_fit_params = bf_return[0]
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
            if g.pl_params[0]:
                try:
                    plot_top_tiers(
                        top_tiers_flo, output_subdir, clust_name, mag_data,
                        col_data, ip_list, err_lst, completeness)
                    print 'Top tier models saved to file and plotted.'
                except:
                    print("  ERROR: top tiers plot could not be generated.")
            else:
                # If top tiers were found but no plot is produced.
                print "Top tier models saved to file."
