
from os.path import isdir, join


def check(mypath, pd):
    """
    Check all parameters related to the generation of synthetic clusters.

    Check that all the required photometric systems evolutionary tracks
    are present with checkEvolTracks()
    """

    # If best fit method is set to run.
    pd['fundam_params_all'] = []
    if pd['best_fit_algor'] != 'n':

        # Check the random seed
        if pd['synth_rand_seed'] in ('n', 'N', 'None', 'none', 'no'):
            pd['synth_rand_seed'] = None
        else:
            try:
                pd['synth_rand_seed'] = int(pd['synth_rand_seed'])
            except ValueError:
                raise ValueError("random seed is not a valid integer")

        pd = getParamVals(pd)

        pd = checkEvolTracks(mypath, pd)

        checkSynthClustParams(pd)

    return pd


def getParamVals(pd):
    """
    Generate the proper lists with the parameter ranges and their priors
    """
    fundam_params_all = {}
    for cl_pars in pd['par_ranges']:
        par_ranges = cl_pars[1:]

        tlst = []
        for idx, param in enumerate(par_ranges):
            # If only one value is defined.
            if len(param) == 1:
                if idx not in (0, 1):
                    tlst.append([float(param[0])])
                else:
                    tlst.append([param[0]])
            else:
                # Store range of values.
                if idx not in (0, 1):
                    tlst.append(list(map(float, param.split('/'))))
                else:
                    tlst.append(param.split('/'))

        # Check for min/max presence in z,log(age) parameters
        for idx in (0, 1):
            p_rng = tlst[idx]
            if len(p_rng) == 1:
                try:
                    p = float(p_rng[0])
                except ValueError:
                    p = p_rng[0]
                    if p not in ('min', 'max'):
                        raise ValueError(
                            "R1: unrecognized string '{}'".format(p))
                tlst[idx] = [p]
            else:
                try:
                    pmin = float(p_rng[0])
                except ValueError:
                    pmin = p_rng[0]
                    if pmin != 'min':
                        raise ValueError(
                            ("R1: unrecognized string '{}'.\nOnly 'min' "
                             + "string is accepted as the lower "
                             + "range.").format(pmin))
                try:
                    pmax = float(p_rng[1])
                except ValueError:
                    pmax = p_rng[1]
                    if pmax != 'max':
                        raise ValueError(
                            ("R1: unrecognized string '{}'.\nOnly 'max' "
                             + "string is accepted as the upper "
                             + "range.").format(pmax))
                tlst[idx] = [pmin, pmax]

        if len(tlst) != 7:
            raise ValueError("Missing parameters in line 'R1'")
        fundam_params_all[cl_pars[0]] = tlst

    pd['fundam_params_all'] = fundam_params_all
    return pd


def checkEvolTracks(mypath, pd):
    """
    Check photometric systems for all the evolutionary tracks that will be
    used. Store and return their IDs, filter names, and paths.
    """
    # Remove duplicate filters (if they exist), and combine them into one
    # tuple per photometric system.
    # The resulting list looks like this:
    # [('phot_syst1', 'T1', 'C'), ('phot_syst2', 'B', 'V'), ...]
    # where the first element of each tuple points to the photometric
    # system, and the remaining elements are the unique filters in that
    # system.
    all_syst_filters = list(set(pd['filters'] + pd['c_filters']))
    d = {}
    for k, v in all_syst_filters:
        d.setdefault(k, [k]).append(v)
    all_syst_filters = sorted(map(tuple, d.values()))

    # # Dictionary of photometric systems defined in the CMD service.
    # all_systs = pd['cmd_systs']

    # Generate name for the isochrones' folders.
    mpi = mypath + 'isochrones/'
    iso_paths = []
    for p_syst in all_syst_filters:
        # Set iso_path according to the above values. The lower-case is
        # important!
        iso_paths.append(join(mpi + p_syst[0].lower()))

    # Check if '/isochrones' folder exists.
    for iso_path in iso_paths:
        if not isdir(iso_path):
            raise ValueError(
                "'Best synthetic cluster fit' function is set to"
                " run but the folder:\n\n {}\n\ndoes not exists."
                .format(iso_path))

    # Extract lambdas from 'filterslambdas' files
    cmd_systs = {}
    for iso_path in iso_paths:
        with open(iso_path + "/filterslambdas.dat") as f:
            ls = f.read().split()[1:]
            phot_syst = iso_path.split('/')[-1]
            # Number of filters/lambdas/omegas
            Nf = int(len(ls) / 3.)
            # Filters
            filters = tuple(ls[:Nf])
            # Lambdas
            lambdas = tuple(map(float, ls[Nf:2 * Nf]))
            cmd_systs[phot_syst] = (filters, lambdas)

    # Check that all filters exist
    for syst in all_syst_filters:
        psyst, filters = syst[0], syst[1:]
        for filt in filters:
            if filt not in cmd_systs[psyst][0]:
                raise ValueError(("Filter '{}' not present in photometric "
                                  "system '{}'".format(filt, psyst)))

    # Read evolutionary tracks used and check that its is the same track used
    # for all systems.
    if len(iso_paths) > 1:
        all_evol_tracks = []
        for iso_path in iso_paths:
            with open(iso_path + "/filterslambdas.dat") as f:
                all_evol_tracks.append(f.read().split()[0])

        if len(set(all_evol_tracks)) > 1:
            raise ValueError(
                "The same evolutionary track must be used\n"
                + "for all the photometric systems. Tracks found:\n"
                + "{}".format(', '.join(all_evol_tracks)))
        evol_track = all_evol_tracks[0]
    else:
        # This is for printing to screen only.
        with open(iso_paths[0] + "/filterslambdas.dat") as f:
            evol_track = f.read().split()[0]

    pd['all_syst_filters'], pd['iso_paths'], pd['cmd_systs'],\
        pd['evol_track'] = all_syst_filters, iso_paths, cmd_systs, evol_track

    return pd


def checkSynthClustParams(pd):
    """
    """
    # Check IMF defined.
    if pd['IMF_name'] not in pd['imf_funcs']:
        raise ValueError("Name of IMF ({}) is incorrect.".format(
            pd['IMF_name']))

    if not 0. <= pd['min_bmass_ratio'] <= 1.:
        raise ValueError(
            "Binary mass ratio set ('{}') is out of\nboundaries. Please select"
            " a value in the range [0., 1.]".format(pd['min_bmass_ratio']))

    # Check maximum magnitude limit defined.
    if isinstance(pd['max_mag'], str):
        if pd['max_mag'] != 'max':
            raise ValueError("Maximum magnitude value selected ({}) is"
                             " not valid.".format(pd['max_mag']))

    for key, vals in pd['fundam_params_all'].items():
        # Check E_BV
        if vals[2][0] < 0.:
            raise ValueError("Minimum extinction must be a positive value")
        # Check d_m
        if vals[3][0] < 0.:
            raise ValueError(
                "Minimum distance modulus must be a positive value")
        # Check mass
        if vals[4][0] < 100.:
            raise ValueError("Minimum mass must be >= 100")
        # Check b_fr
        if vals[5][0] < 0. or vals[5][0] > 1.:
            raise ValueError("Binary fraction must be in the range [0, 1]")
        # Check R_V
        if vals[6][0] <= 0.:
            raise ValueError("Minimum R_V must be > 0")
