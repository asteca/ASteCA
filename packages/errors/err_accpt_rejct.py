
import numpy as np
import warnings


def max_err_cut(cld, err_max):
    """
    Accept stars with photometric / kinematic errors < e_max in its
    magnitudes, colors, parallax, proper motions, and radial velocity.

    All 'nan' values are kept, i.e.: evaluated to True.
    Source: https://stackoverflow.com/a/48584644/1391441

    This means that there will be a different number of stars with valid data
    in each data dimension.
    """

    # Prepare values.
    em_float = []
    for err in err_max:
        if err == 'n':
            em_float.append(np.inf)
        else:
            em_float.append(float(err))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)

        # Photometric data.
        m_msk = np.logical_or(cld['em'] < em_float[0], np.isnan(cld['em']))
        c_msk = np.logical_or(cld['ec'] < em_float[0], np.isnan(cld['ec']))

        # Kinematic data.
        plx_msk = np.logical_or(
            cld['ek'][0] < em_float[1], np.isnan(cld['ek'][0]))
        pmx_msk = np.logical_or(
            cld['ek'][1] < em_float[2], np.isnan(cld['ek'][1]))
        pmy_msk = np.logical_or(
            cld['ek'][2] < em_float[2], np.isnan(cld['ek'][2]))
        rv_msk = np.logical_or(
            cld['ek'][3] < em_float[3], np.isnan(cld['ek'][3]))

    acpt_indx = np.flatnonzero(
        (m_msk.all(0) & c_msk.all(0) & plx_msk & pmx_msk & pmy_msk &
         rv_msk)).tolist()
    rjct_indx = np.setdiff1d(
        np.arange(len(cld['em'][0])), acpt_indx).tolist()

    return acpt_indx, rjct_indx, em_float


def main(i_c, cld, clp, err_max, **kwargs):
    """
    Accept and reject stars in and out of the cluster's boundaries according to
    a given `err_max` error value.
    """

    # Call function to reject stars with errors > e_max.
    acpt_indx, rjct_indx, em_float = max_err_cut(cld, err_max)
    if not acpt_indx:
        raise ValueError(
            "ERROR: No stars left after error rejection.\n"
            "Try increasing the maximum accepted error value.")

    # Filter elements.
    # In place for #243
    import sys
    if sys.version_info[0] == 2:
        acpt = {k: v[..., acpt_indx] for k, v in cld.iteritems()}
        rjct = {k: v[..., rjct_indx] for k, v in cld.iteritems()}
    else:
        acpt = {k: v[..., acpt_indx] for k, v in cld.items()}
        rjct = {k: v[..., rjct_indx] for k, v in cld.items()}

    # Store each star separately. This part is important since it is here
    # where we define the position of the data.
    acpt_stars = [
        list(_) for _ in zip(*[
            acpt['ids'], acpt['x'], acpt['y'], acpt['mags'].T, acpt['em'].T,
            acpt['cols'].T, acpt['ec'].T, acpt['kine'].T, acpt['ek'].T])]
    rjct_stars = [
        list(_) for _ in zip(*[
            rjct['ids'], rjct['x'], rjct['y'], rjct['mags'].T, rjct['em'].T,
            rjct['cols'].T, rjct['ec'].T, rjct['kine'].T, rjct['ek'].T])]

    print("  Stars rejected based on their errors ({}).".format(
        len(rjct_stars)))

    if i_c == 'comp':
        clp['em_float'] = em_float
        clp['acpt_stars_c'], clp['rjct_stars_c'] = acpt_stars, rjct_stars
    elif i_c == 'incomp':
        clp['acpt_stars_i'], clp['rjct_stars_i'] = acpt_stars, rjct_stars

    return clp
