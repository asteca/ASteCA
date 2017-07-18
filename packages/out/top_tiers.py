
from os.path import join
from time import strftime
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from .._version import __version__
from ..synth_clust import synth_cl_plot
from ..synth_clust import imf
import top_tiers_plot
import prep_plots


def top_tiers_file(output_subdir, clust_name, best_model, top_tiers, X_lkl,
                   X_par):
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
# The 'Top tiers' models present a difference of <{}% with the best \
likelihood\n# value found for the 'Best model' and one or more parameters \
differing\n# >{}% with the best value assigned for said parameter(s).\n\
#\n\n\
# Best model\n\
#\n\
# Lkl               z   log(age)    E(B-V)      dist_m     M (Mo)       b_fr\n\
{:<10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}\n\n\
# Top tiers\n\
#\n\
# Lkl               z   log(age)    E(B-V)      dist_m     M (Mo)       b_fr\n"
                        .format(__version__, now_time, int(X_lkl * 100),
                                int(X_par * 100), *best_model))
    out_data_file.close()

    # Write values to file.
    with open(out_file_name, "a") as f_out:
        for mod in top_tiers:
            f_out.write('''{:<10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}'''
                        .format(*mod))
            f_out.write('\n')


def plot_top_tiers(pd, top_tiers_flo, output_subdir, clust_name, mags,
                   cols, isoch_fit_params, ext_coefs, N_fc, err_lst,
                   completeness, max_mag_syn):
    '''
    Plot all top tiers.
    '''
    e_max = pd['er_params'][1]
    # Obtain mass distribution using the selected IMF.
    st_dist_mass = imf.main(pd['IMF_name'], pd['m_high'])

    # figsize(x1, y1), GridSpec(y2, x2) --> To have square plots: x1/x2 =
    # y1/y2 = 2.5
    fig = plt.figure(figsize=(25, 10))  # create the top-level container
    gs = gridspec.GridSpec(4, 10)      # create a GridSpec object
    # Add version number to top left.
    ver = '[ASteCA ' + __version__ + ']'
    x_coord = 0.95 - (len(__version__) - 6) * 0.001
    plt.figtext(x_coord, .986, ver, fontsize=9, color='#585858')

    x_ax, y_ax, y_axis = prep_plots.ax_names(pd['filters'], pd['colors'])

    synth_cls = []
    for mod in top_tiers_flo:
        # Call function to generate synthetic cluster.
        shift_isoch, synth_clst = synth_cl_plot.main(
            e_max, pd['bin_mr'], pd['fundam_params'], pd['theor_tracks'],
            isoch_fit_params, err_lst, completeness, max_mag_syn,
            st_dist_mass, pd['R_V'], ext_coefs, N_fc)

        # TODO using first magnitude and color defined
        first_col, first_mag = synth_clst[0][0][1], synth_clst[0][0][0]
        x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd = prep_plots.diag_limits(
            y_axis, first_col, first_mag)

        # Store plotting parameters.
        synth_cls.append([gs, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd,
                         x_ax, y_ax, synth_clst, mod, shift_isoch])

    # Generate all plots.
    for n, args in enumerate(synth_cls, 1):
        top_tiers_plot.plot(n, *args)

    fig.tight_layout()
    # Generate output file for each data file.
    pl_fmt, pl_dpi = pd['pl_params'][1:3]
    plt.savefig(join(output_subdir, clust_name + '_D3.' + pl_fmt),
                dpi=pl_dpi)
    # Close to release memory.
    plt.clf()
    plt.close()


def main(npd, cld, pd, err_lst, ext_coefs, N_fc, completeness, max_mag_syn,
         isoch_fit_params, **kwargs):
    '''
    Obtain top tier models, produce data file and output image.
    '''

    if pd['bf_flag']:
        clust_name, output_subdir = npd['clust_name'], npd['output_subdir']
        all_models = isoch_fit_params[-1]

        # Sort all models/solutions by their (minimum) likelihood values.
        all_m_sort = sorted(zip(*all_models), key=lambda x: x[1])
        best_mod, best_lik = all_m_sort[0][0], all_m_sort[0][1]
        # 'X_lkl' is the percentage that determines how close to the best
        # likelihood a 'top tier' must be.
        X_lkl = 0.1
        best_lik_X = best_lik * X_lkl

        # 'X_par' is the percentage that determines how far away from the best
        # match parameter value a 'top tier's parameter value must be.
        X_par = 0.1
        top_tiers_str, top_tiers_flo = [], []
        # Compare every model/solution with the best one.
        for mod in all_m_sort[1:]:
            # If the likelihood is closer than X% to the min likelihood.
            if (mod[1] - best_lik) < best_lik_X:
                # If any parameter has a difference larger than X_par%.
                mod_flag, p_mark = False, []
                for i, p in enumerate(mod[0]):
                    if abs(p - best_mod[i]) > best_mod[i] * X_par:
                        mod_flag = True
                        p_mark.append(str(p) + '*')
                    else:
                        p_mark.append(str(p))
                if mod_flag:
                    if mod[0] not in top_tiers_flo:
                        top_tiers_str.append(
                            ['{:.3f}'.format(mod[1])] + p_mark)
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
                           top_tiers_str, X_lkl, X_par)
            print("Top tier models saved to file.")

            # Create output image if flag is 'true'.
            flag_make_plot = pd['pl_params'][0]
            if flag_make_plot:
                try:
                    plot_top_tiers(
                        pd, top_tiers_flo, output_subdir, clust_name,
                        cld['mags'], cld['cols'], isoch_fit_params, ext_coefs,
                        N_fc, err_lst, completeness, max_mag_syn)
                    print("<<Plots from 'D3' block created>>")
                except:
                    import traceback
                    print traceback.format_exc()
                    print("  ERROR: top tiers plot could not be generated.")
