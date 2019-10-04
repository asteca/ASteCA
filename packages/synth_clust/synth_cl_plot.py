
import numpy as np
from . import move_isochrone
from . import synth_cluster
from . import zaWAverage


def main(best_fit_algor, fundam_params, isoch_fit_params, synthcl_args):
    """
    Generate shifted isochrone and synthetic cluster for plotting.
    """
    if best_fit_algor == 'boot+GA':
        # Use ML fit values for all parameters.
        model = isoch_fit_params['map_sol']
    elif best_fit_algor == 'ptemcee':
        # Use mean fit values for all parameters.
        model = isoch_fit_params['mean_sol']

    theor_tracks = synthcl_args[0]
    R_V, ext_coefs, N_fc = synthcl_args[6:9]
    E_BV, dm = model[2], model[3]

    # Generate a model with the fitted parameters only.
    model_var = np.array(model)[isoch_fit_params['varIdxs']]

    # Generate shifted best fit isochrone.
    isochrone, _ = zaWAverage.main(
        theor_tracks, fundam_params, isoch_fit_params['varIdxs'], model_var)
    shift_isoch = move_isochrone.main(
        isochrone, E_BV, dm, R_V, ext_coefs, N_fc)

    # Generate best fit synthetic cluster.
    synth_clst = synth_cluster.main(isochrone, model[2:], *synthcl_args[1:])

    return shift_isoch, synth_clst
