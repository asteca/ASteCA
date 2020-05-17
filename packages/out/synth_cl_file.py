
import numpy as np
from ..synth_clust import synth_cluster


def main(clp, npd, best_fit_algor, fundam_params, filters, colors,
         theor_tracks, R_V, m_ini_idx, binar_flag, **kwargs):
    """
    Create output data file with stars in the best fit synthetic cluster found
    by the 'Best Fit' function.
    """

    clp['synth_clst_plot'], clp['binar_idx_plot'], clp['shift_isoch'],\
        clp['synthcl_Nsigma'] = [np.array([]) for _ in range(4)]

    if best_fit_algor != 'n':
        isoch_moved, synth_clst, sigma, extra_pars, synthcl_Nsigma =\
            synth_cl_plot(
                fundam_params, clp['isoch_fit_params'],
                clp['isoch_fit_errors'], theor_tracks, clp['completeness'],
                clp['max_mag_syn'], clp['st_dist_mass'], R_V, clp['ext_coefs'],
                clp['N_fc'], clp['err_pars'], m_ini_idx, binar_flag)

        # If cluster is not empty.
        if synth_clst.any():
            # Prepare data.
            e_mags_cols = sigma.T
            binar_idx, ini_mass = extra_pars[0], extra_pars[2]
            # Prepare header.
            hdr = ['ID  '] + [f[1] + '   ' for f in filters]
            hdr += ['(' + c[1].replace(',', '-') + ')   ' for c in colors]
            hdr += ['e_' + f[1] + '   ' for f in filters]
            hdr += ['e_(' + c[1].replace(',', '-') + ')   ' for c in colors]
            # TODO storing only 'M_ini' for now
            hdr = ''.join(hdr) + 'Mini\n'

            # Save best fit synthetic cluster found to file.
            with open(npd['synth_file_out'], "w") as f_out:
                f_out.write(hdr)
                for i, st in enumerate(synth_clst):
                    if binar_idx[i] > 1.:
                        ID = '2' + str(i)
                    else:
                        ID = '1' + str(i)
                    f_out.write("{:<8}".format(ID))
                    for _ in st:
                        f_out.write("{:<8.4f}".format(_))
                    for _ in e_mags_cols[i]:
                        f_out.write("{:<8.4f}".format(_))
                    # TODO storing only 'M_ini' for now
                    # f_out.write("{:<8.4f}\n".format(*extra_pars[i]))
                    f_out.write("{:<8.4f}\n".format(ini_mass[i]))

            # Save for plotting.
            clp['synth_clst_plot'], clp['binar_idx_plot'],\
                clp['shift_isoch'], clp['synthcl_Nsigma'] = synth_clst,\
                binar_idx, isoch_moved, synthcl_Nsigma

        else:
            print("  ERROR: empty synthetic cluster could not be saved\n"
                  "  to file")

    return clp


def synth_cl_plot(
    fundam_params, isoch_fit_params, isoch_fit_errors, theor_tracks,
    completeness, max_mag_syn, st_dist_mass, R_V, ext_coefs, N_fc, err_pars,
        m_ini_idx, binar_flag, N_models=1000):
    """
    Generate a synthetic cluster and a N(mu, sigma) region for plotting, from
    'N_models' models.
    """

    # DEPRECATED May 2020
    # if best_fit_algor == 'boot+GA':
    #     # Use ML fit values for all parameters.
    #     model = isoch_fit_params['map_sol']
    # elif best_fit_algor in ('ptemcee', 'emcee'):
    # Use mean fit values for all parameters.
    model = isoch_fit_params['mean_sol']

    # Generate a model with the "best" fitted parameters.
    model_var = np.array(model)[isoch_fit_params['varIdxs']]

    # Pack common args.
    syntClustArgs = (
        fundam_params, isoch_fit_params['varIdxs'], theor_tracks,
        completeness, max_mag_syn, st_dist_mass, R_V, ext_coefs,
        N_fc, err_pars, m_ini_idx, binar_flag)

    isoch_moved, synth_clust, sigma, extra_pars = setSynthClust(
        model_var, True, *syntClustArgs)

    # Generate the synthetic clusters required to plot the 1-sigma zone.
    p_vals, nancount = [], 0
    for i, p in enumerate(model):
        std = isoch_fit_errors[i][-1]
        if not np.isnan(std):
            vals = [p, std]
            p_vals.append(vals)
        else:
            nancount += 1

    # Check if any parameter has an uncertainty attached.
    if nancount == 6:
        synthcl_Nsigma = np.array([])
        return isoch_moved, synth_clust, sigma, extra_pars, synthcl_Nsigma

    # Generate 'N_models' random models.
    models = []
    for par in p_vals:
        models.append(np.random.normal(par[0], par[1], N_models))

    # Generate the synthetic clusters from the sampled parameters.
    synthcl_Nsigma = [[] for _ in range(sum(N_fc))]
    for model in np.array(models).T:
        # Synthetic cluster
        sc = setSynthClust(model, True, *syntClustArgs)[1].T
        if sc.any():
            for i, photom_dim in enumerate(sc):
                synthcl_Nsigma[i] += list(photom_dim)
    synthcl_Nsigma = np.array(synthcl_Nsigma)

    return isoch_moved, synth_clust, sigma, extra_pars, synthcl_Nsigma


def setSynthClust(
    model, extra_pars_flag, fundam_params, varIdxs, theor_tracks,
    completeness, max_mag_syn, st_dist_mass, R_V, ext_coefs, N_fc, err_pars,
        m_ini_idx, binar_flag):
    """
    Generate synthetic cluster given by 'model'.
    """
    synth_clust, sigma, extra_pars, isoch_moved, mass_dist, isoch_binar,\
        isoch_compl = synth_cluster.main(
            fundam_params, varIdxs, model, theor_tracks, completeness,
            max_mag_syn, st_dist_mass, R_V, ext_coefs, N_fc, err_pars,
            m_ini_idx, binar_flag, extra_pars_flag)

    return isoch_moved, synth_clust, sigma, extra_pars
