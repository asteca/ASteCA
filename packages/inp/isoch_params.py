
import numpy as np
import read_isochs


def interp_isoch(isochrone, N=1500):
    '''
    Interpolate extra values for all the parameters in the theoretic
    isochrones.
    '''
    t, xp = np.linspace(0, 1, N), np.linspace(0, 1, len(isochrone[0]))
    # Store isochrone's interpolated values.
    isoch_inter = np.asarray([np.interp(t, xp, _) for _ in isochrone])

    return isoch_inter


def main(pd):
    '''
    Read isochrones and parameters if best fit function is set to run.
    '''
    ip_list = []
    # Only read files of best fit method is set to run.
    if pd['bf_flag']:
        # Print info about tracks.
        print("Process {} theoretical isochrones".format(
            pd['cmd_evol_tracks'][pd['evol_track']][1]))

        for syst in pd['all_syst_filters']:
            print("in the '{}' photometric system.\n".format(
                pd['cmd_systs'][syst[0]][0]))

        # Get isochrones and their parameter values.
        isoch_list = read_isochs.main(pd['met_f_filter'], pd['age_values'])

        # Interpolate extra points into all the isochrones.
        isochs_interp = [[] for _ in isoch_list]
        for i, _ in enumerate(isoch_list):
            for isoch in _:
                isochs_interp[i].append(interp_isoch(isoch))

        # Pack parameterss.
        param_values = [pd['met_values'], pd['age_values']] +\
            pd['param_ranges'][2:]
        ip_list = [isochs_interp, param_values]

        # Obtain number of models in the solutions space.
        lens = [len(_) for _ in param_values]
        total = reduce(lambda x, y: x * y, lens, 1)
        print(
            "Number of values per parameter:\n"
            "  {} metallicity values (z),\n"
            "  {} age values (per z),\n"
            "  {} reddening values,\n"
            "  {} distance values,\n"
            "  {} mass values,\n"
            "  {} binary fraction values.".format(*lens))
        print("  = {:.1e} approx total models.\n".format(total))

    pd['ip_list'] = ip_list
    return pd
