
# import numpy as np
# import matplotlib.pyplot as plt
# import display_errors  # DEPRECATED
import err_accpt_rejct_max as e_a_r_mx


# def err_sel_stars(acpt_indx, rjct_indx, cld):
#     '''
#     Select accepted/rejected stars.
#     '''
#     mags_z, em_z = np.array(zip(*cld['mags'])), np.array(zip(*cld['em']))
#     cols_z, ec_z = np.array(zip(*cld['cols'])), np.array(zip(*cld['ec']))

#     idx_a = np.array(cld['ids'])[acpt_indx].tolist()
#     x_a = np.array(cld['x'])[acpt_indx].tolist()
#     y_a = np.array(cld['y'])[acpt_indx].tolist()
#     mags_a = mags_z[acpt_indx].tolist()
#     em_a = em_z[acpt_indx].tolist()
#     cols_a = cols_z[acpt_indx].tolist()
#     ec_a = ec_z[acpt_indx].tolist()

#     idx_r = np.array(cld['ids'])[rjct_indx].tolist()
#     x_r = np.array(cld['x'])[rjct_indx].tolist()
#     y_r = np.array(cld['y'])[rjct_indx].tolist()
#     mags_r = mags_z[rjct_indx].tolist()
#     em_r = em_z[rjct_indx].tolist()
#     cols_r = cols_z[rjct_indx].tolist()
#     ec_r = ec_z[rjct_indx].tolist()

#     # Store everything as lists.
#     acpt_stars = [
#         list(_) for _ in zip(*[idx_a, x_a, y_a, mags_a, em_a, cols_a, ec_a])]
#     rjct_stars = [
#         list(_) for _ in zip(*[idx_r, x_r, y_r, mags_r, em_r, cols_r, ec_r])]

#     return acpt_stars, rjct_stars


# DEPRECATED
# def maxError(cld):
#     """
#     Find the maximum error value across magnitudes and colors.
#     """
#     e_max_all = 0.
#     for em in cld['em']:
#         e_max_all = max(e_max_all, max(em))
#     for ec in cld['ec']:
#         e_max_all = max(e_max_all, max(ec))

#     return e_max_all


def main(i_c, cld, clp, err_max, **kwargs):
    """
    Accept and reject stars in and out of the cluster's boundaries according to
    a given `err_max` photometric error value.
    """

    # if pd['run_mode'] in ('auto', 'semi'):
    # e_max_val = pd['err_max']
    # Call function to reject stars with errors > e_max.
    acpt_indx, rjct_indx = e_a_r_mx.main(cld, err_max)
    if not acpt_indx:
        raise ValueError(
            "ERROR: No stars left after error rejection.\n"
            "Try increasing the maximum accepted error value.")

    # Call function to store stars according to the returned indexes.
    # acpt_stars, rjct_stars = err_sel_stars(acpt_indx, rjct_indx, cld)

    # Filter elements.
    acpt = {k: v[..., acpt_indx] for k, v in cld.iteritems()}
    rjct = {k: v[..., rjct_indx] for k, v in cld.iteritems()}

    # Store each star separately.
    acpt_stars = [
        list(_) for _ in zip(*[
            acpt['ids'], acpt['x'], acpt['y'], acpt['mags'].T, acpt['em'].T,
            acpt['cols'].T, acpt['ec'].T])]
    rjct_stars = [
        list(_) for _ in zip(*[
            rjct['ids'], rjct['x'], rjct['y'], rjct['mags'].T, rjct['em'].T,
            rjct['cols'].T, rjct['ec'].T])]

    # DEPRECATED (at least for now, 08/05/18)
    # # If 'manual' mode is set, display errors distributions and ask the user
    # # to accept it or else use all stars except those with errors > e_max in
    # # either the magnitude or the color.
    # elif pd['run_mode'] == 'manual':
    #     e_max_all = maxError(cld)
    #     move_on = False
    #     while not move_on:
    #         while True:
    #             answer_rad = int(raw_input(
    #                 "Choose a method for error-based stars rejection:\n"
    #                 "  1 (emax), 2 (use all stars): "))

    #             if answer_rad == 1:
    #                 e_max_val = float(raw_input(
    #                     'Select maximum error value: '))
    #                 # Call function to reject stars with errors > e_max.
    #                 acpt_indx, rjct_indx = e_a_r_mx.main(cld, e_max_val)
    #                 break
    #             elif answer_rad == 2:
    #                 # Store all indexes.
    #                 acpt_indx, e_max_val = range(len(mmag)), e_max_all
    #                 break
    #             if answer_rad not in (1, 2):
    #                 print("Sorry, input is not valid. Try again.")
    #             else:
    #                 if len(acpt_indx) == 0:
    #                     print("No stars left after error rejection. Try"
    #                           " again.")
    #                 else:
    #                     break

    #         # Call function to store stars according to the returned indexes.
    #         acpt_stars, rjct_stars = err_sel_stars(acpt_indx, rjct_indx, cld)

    #         if answer_rad == 1:
    #             # Display automatic errors rejection.
    #             display_errors.main(
    #                 pd['filters'], pd['colors'], mmag, acpt_stars,
    #                 rjct_stars, e_max_val)
    #             plt.show()
    #             # Ask if keep or reject.
    #             while True:
    #                 try:
    #                     mv_on = raw_input('Accept error based rejection?'
    #                                       ' (y/n) ')
    #                     if mv_on == 'y':
    #                         print("Fit accepted.")
    #                         move_on = True
    #                         break
    #                     elif mv_on == 'n':
    #                         print("Fit rejected.")
    #                         break
    #                     else:
    #                         print("Sorry, input is not valid. Try again.")
    #                 except Exception:
    #                     print("Sorry, input is not valid. Try again.")
    #         else:
    #             move_on = True

    print("  Stars rejected based on their errors ({}).".format(
        len(rjct_stars)))

    if i_c == 'incomp':
        clp['acpt_stars_i'], clp['rjct_stars_i'] = acpt_stars, rjct_stars
    elif i_c == 'comp':
        clp['acpt_stars_c'], clp['rjct_stars_c'] = acpt_stars, rjct_stars

    return clp
