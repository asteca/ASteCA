
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
    elif best_fit_algor in ('ptemcee', 'emcee'):
        # Use mean fit values for all parameters.
        model = isoch_fit_params['mean_sol']

    # Generate a model with the "best" fitted parameters.
    model_var = np.array(model)[isoch_fit_params['varIdxs']]
    # Generate best fit synthetic cluster.
    synth_clst = synth_cluster.main(
        fundam_params, isoch_fit_params['varIdxs'], model_var, *synthcl_args)

    # Generate shifted best fit isochrone.
    theor_tracks = synthcl_args[0]
    R_V, ext_coefs, N_fc = synthcl_args[6:9]
    E_BV, dm = model[2], model[3]
    _, z_model, a_model, ml, mh, al, ah = synth_cluster.properModel(
        fundam_params, model_var, isoch_fit_params['varIdxs'])
    isochrone = zaWAverage.main(
        theor_tracks, fundam_params, z_model, a_model, ml, mh, al, ah)
    shift_isoch = move_isochrone.main(
        isochrone, E_BV, dm, R_V, ext_coefs, N_fc)

    return shift_isoch, synth_clst
