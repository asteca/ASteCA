
import numpy as np
import read_isochs
import met_ages_values
from ..check import CMD_phot_systs_filts


def interp_isoch(isochrone):
    '''
    Interpolate extra color, magnitude and masses into the isochrone.
    '''
    N = 1500
    t, xp = np.linspace(0, 1, N), np.linspace(0, 1, len(isochrone[0]))
    # Store isochrone's interpolated values.
    isoch_inter = np.asarray([np.interp(t, xp, _) for _ in isochrone])

    return isoch_inter


def main(cmd_evol_tracks, evol_track, bf_params, iso_paths, par_ranges,
         all_syst_filters, **kwargs):
    '''
    Read isochrones and parameters if best fit function is set to run.
    '''

    ip_list = []
    # Only read files of best fit method is set to run.
    bf_flag = bf_params[0]
    if bf_flag:
        # Print info about tracks.
        print("Process {} theoretical isochrones".format(
            cmd_evol_tracks[evol_track][1]))

        for syst in all_syst_filters:
            print("in the '{}' photometric system.\n".format(
                CMD_phot_systs_filts.main()[syst[0]][0]))

        # Obtain allowed metallicities and ages. Use the first photometric
        # system defined.
        # *WE ASUME ALL PHOTOMETRIC SYSTEMS CONTAIN THE SAME NUMBER OF
        # METALLICITY FILES*
        param_ranges, met_f_filter, met_values, age_values =\
            met_ages_values.main(iso_paths, par_ranges)

        # Get isochrones and their parameter values.
        isoch_list = read_isochs.main(met_f_filter, age_values)

        # Interpolate extra points into all the isochrones.
        isochs_interp = [[] for _ in isoch_list]
        for i, _ in enumerate(isoch_list):
            for isoch in _:
                isochs_interp[i].append(interp_isoch(isoch))

        # Pack params.
        param_values = [met_values, age_values] + param_ranges[2:]
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

    return ip_list
