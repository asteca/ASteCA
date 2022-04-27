
import numpy as np
from ..best_fit.bf_common import getSynthClust
from . add_errors import getSigmas


def main(clp, pd, td):
    """
    Generate the arrays 'synth_cl_phot', 'synth_cl_sigma', 'binar_idx'.

    synth_cl_phot  : a sample of the "best" synthetic cluster fit. Used for
    plotting and storing to file.
    synth_cl_sigma : uncertainties attached to the above sample. Used for
    the output file
    binar_idx      : indexes of binary systems. Used for plotting.
    """

    if pd['best_fit_algor'] == 'n':
        return clp

    # Use the selected solution values for all the parameters.
    model = clp['isoch_fit_params'][pd['D3_sol'] + '_sol']

    # Generate isochrone, synthetic cluster (with uncertainties), and sigma
    # values for the "best" fitted parameters.
    # Model with the "best" fitted parameters.
    model_var = np.array(model)[clp['varIdxs']]
    # shape: (N_dim, N_stars)
    synth_clust = getSynthClust(model_var, *clp['syntClustArgs'], False)[0]
    # Get uncertainties
    sigma = []
    for i, popt_mc in enumerate(clp['err_lst']):
        sigma.append(getSigmas(synth_clust[0], popt_mc))

    # Save for plotting and storing.
    clp['synth_cl_phot'], clp['synth_cl_sigma'] = synth_clust, sigma

    # Mask that points to *binary* systems. 'synth_clust[-1]' is the
    # mass of binary systems, where a value of 0 indicates a single
    # system.
    clp['binar_idx'] = synth_clust[-1] > 0.

    return clp
