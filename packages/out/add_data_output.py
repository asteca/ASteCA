
import os
from astropy.table import Table
from ..inp.get_data import flatten


def main(
    npd, pd, cont_index, n_memb, frac_cl_area, kde_cent,
    clust_rad, KP_Bys_rc, KP_Bys_rt, KP_memb_num, isoch_fit_params,
        isoch_fit_errors, MassT_dist_vals, binar_dist_vals, **kwargs):
    """
    Add data obtained to the 'asteca_output.dat' file.
    """

    # Unpack data.
    out_file_name, write_name = npd['out_file_name'], npd['write_name']

    # Round structure parameters.
    frmt = "{:.6f}"
    # Center + radii and uncertainties
    cre_r = [
        frmt.format(_) for _ in [
            kde_cent[0], kde_cent[1], clust_rad,
            KP_Bys_rc['median'], KP_Bys_rc['16th'], KP_Bys_rc['84th'],
            KP_Bys_rt['median'], KP_Bys_rt['16th'], KP_Bys_rt['84th']]]

    # Cluster parameters and errors.
    cpe_r = [
        MassT_dist_vals['mean_sol'], MassT_dist_vals['median_sol'],
        MassT_dist_vals['mode_sol'], MassT_dist_vals['errors'][0],
        MassT_dist_vals['errors'][1], MassT_dist_vals['errors'][2]]
    # Add (z, a, beta, E_BV, DR, RV, dm)
    cpe_r += [
        item for t in zip(
            isoch_fit_params['mean_sol'], isoch_fit_params['median_sol'],
            isoch_fit_params['mode_sol'],
            list(zip(*isoch_fit_errors))[0], list(zip(*isoch_fit_errors))[1],
            list(zip(*isoch_fit_errors))[2])
        for item in t]
    # Insert b_fr before E_BV
    b_fr = [
        binar_dist_vals['mean_sol'], binar_dist_vals['median_sol'],
        binar_dist_vals['mode_sol'], binar_dist_vals['errors'][0],
        binar_dist_vals['errors'][1], binar_dist_vals['errors'][2]]
    idx = 3 * 6  # 3 pars * 6 vals
    cpe_r = cpe_r[:idx] + b_fr + cpe_r[idx:]
    cpe_r = [round(_, 5) for _ in cpe_r]

    # Store all parameter values in list.
    line = [
        write_name, cre_r, round(cont_index, 2), KP_memb_num,
        n_memb, round(frac_cl_area, 2), cpe_r]
    # Flatten list.
    line_f = list(flatten(line))

    with open(out_file_name, mode='a') as f:
        # Some platforms don't automatically seek to end when files opened
        # in append mode
        f.seek(0, os.SEEK_END)
        tt = Table(list(zip(*[line_f])))
        tt.write(f, format='ascii.no_header', delimiter=',')

    print("Analysis results added to output file")
