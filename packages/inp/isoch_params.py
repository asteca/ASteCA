
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


def main(pd, param_ranges, met_f_filter, met_values, age_values):
    '''
    Read isochrones and parameters if best fit function is set to run.
    '''

    pd['isochs_theor'], pd['param_values'] = [], []
    # Only read files of best fit method is set to run.
    if pd['bf_flag']:
        # Print info about tracks.
        print("Process {} theoretical isochrones".format(
            pd['cmd_evol_tracks'][pd['evol_track']][1]))

        for syst in pd['all_syst_filters']:
            print("in the '{}' photometric system.\n".format(
                pd['cmd_systs'][syst[0]][0]))

        # Get isochrones and their extra parameters (mass, etc.).
        isoch_list, extra_pars = read_isochs.main(met_f_filter, age_values,
                                                  **pd)

        # Interpolate extra points into all the isochrones.
        isochs_interp = [[] for _ in isoch_list]
        for i, _ in enumerate(isoch_list):
            for isoch in _:
                isochs_interp[i].append(interp_isoch(isoch))

        # Take the synthetic data from the unique filters read, create the
        # necessary colors, and position the magnitudes and colors in the
        # same sense they are read from the cluster's data file.
        pd['isochs_theor'] = arrange_filters()

        # Pack parameters.
        pd['param_values'] = [met_values, age_values] + param_ranges[2:]

        # Obtain number of models in the solutions space.
        lens = [len(_) for _ in pd['param_values']]
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

        return pd
