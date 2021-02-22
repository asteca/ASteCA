
from os.path import isdir, join


def check(mypath, pd):
    """
    Check all parameters related to the generation of synthetic clusters.

    Check that all the required photometric systems evolutionary tracks
    are present with checkEvolTracks()
    """

    # If best fit method is set to run.
    pd['fundam_params'] = []
    if pd['best_fit_algor'] != 'n':

        # Check the random seed
        if pd['synth_rand_seed'] in ('n', 'N', 'None', 'none', 'no'):
            pd['synth_rand_seed'] = None
        else:
            try:
                pd['synth_rand_seed'] = int(pd['synth_rand_seed'])
            except ValueError:
                raise ValueError("random seed is not a valid integer")

        checkParamRanges(pd)

        pd = checkEvolTracks(mypath, pd)

        checkSynthClustParams(pd)

        pd = getParamVals(pd)

    return pd


def checkParamRanges(pd):
    """
    """
    # Check that no parameter range is empty.
    p_names = [
        'metallicity', 'age', 'extinction', 'distance', 'mass', 'binarity']
    for i, p in enumerate(pd['par_ranges']):
        if not p:
            raise ValueError("Range defined for '{}' parameter is"
                             " empty.".format(p_names[i]))
    m_range = pd['par_ranges'][4]
    if m_range[0] == 0.:
        print("  WARNING: minimum total mass is zero in params_input"
              " file.\n  Changed minimum mass to 10.\n")
        m_range[0] = 10.


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

    # Check R_V defined.
    if pd['R_V'] <= 0.:
        raise ValueError(
            "Ratio of total to selective absorption\nR_V ({}) must be positive"
            " defined.".format(pd['R_V']))

    # Check maximum magnitude limit defined.
    if isinstance(pd['max_mag'], str):
        if pd['max_mag'] != 'max':
            raise ValueError("Maximum magnitude value selected ({}) is"
                             " not valid.".format(pd['max_mag']))


def getParamVals(pd):
    """
    Properly format parameter ranges to be used by the selected best fit
    method.
    """

    fundam_params = []
    for i, param in enumerate(pd['par_ranges']):
        # If only one value is defined.
        if len(param) == 1:
            fundam_params.append([param[0]])
        # If min == max store single value in array.
        elif param[0] == param[1]:
            fundam_params.append([param[0]])
        else:
            # Store range of values.
            fundam_params.append(param)

    # Check for min/max presence in z,log(age) parameters
    for idx in (0, 1):
        p_rng = fundam_params[idx]
        t = 'RZ' if idx == 0 else 'RA'
        if len(p_rng) == 1:
            try:
                p = float(p_rng[0])
            except ValueError:
                p = p_rng[0]
                if p not in ('min', 'max'):
                    raise ValueError("'{}': unrecognized string '{}'".format(
                        t, p))
            fundam_params[idx] = [p]
        else:
            try:
                pmin = float(p_rng[0])
            except ValueError:
                pmin = p_rng[0]
                if pmin != 'min':
                    raise ValueError(
                        ("'{}': unrecognized string '{}'.\nOnly 'min' "
                         + "string is accepted as the lower range.").format(
                            t, pmin))
            try:
                pmax = float(p_rng[1])
            except ValueError:
                pmax = p_rng[1]
                if pmax != 'max':
                    raise ValueError(
                        ("'{}': unrecognized string '{}'.\nOnly 'max' "
                         + "string is accepted as the upper range.").format(
                            t, pmax))
            fundam_params[idx] = [pmin, pmax]

    pd['fundam_params'] = fundam_params
    return pd
