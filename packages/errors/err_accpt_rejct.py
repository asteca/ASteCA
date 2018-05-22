
import err_accpt_rejct_max as e_a_r_mx


def main(i_c, cld, clp, err_max, **kwargs):
    """
    Accept and reject stars in and out of the cluster's boundaries according to
    a given `err_max` photometric error value.
    """

    # Call function to reject stars with errors > e_max.
    acpt_indx, rjct_indx = e_a_r_mx.main(cld, err_max)
    if not acpt_indx:
        raise ValueError(
            "ERROR: No stars left after error rejection.\n"
            "Try increasing the maximum accepted error value.")

    # Filter elements.
    acpt = {k: v[..., acpt_indx] for k, v in cld.iteritems()}
    rjct = {k: v[..., rjct_indx] for k, v in cld.iteritems()}

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

    if i_c == 'incomp':
        clp['acpt_stars_i'], clp['rjct_stars_i'] = acpt_stars, rjct_stars
    elif i_c == 'comp':
        clp['acpt_stars_c'], clp['rjct_stars_c'] = acpt_stars, rjct_stars

    return clp
