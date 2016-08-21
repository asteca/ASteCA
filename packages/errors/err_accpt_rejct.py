
import matplotlib.pyplot as plt
import display_errors
import err_accpt_rejct_lowexp as e_a_r_le
import err_accpt_rejct_eyefit as e_a_r_ef
import err_accpt_rejct_max as e_a_r_mx


def err_sel_stars(acpt_indx, rjct_indx, cld):
    '''
    Go through the lists of indexes to select accepted/rejected stars.
    '''
    # Unpack data.
    ids, x, y, mags, em, cols, ec = cld['ids'], cld['x'], cld['y'],\
        cld['mags'], cld['em'], cld['cols'], cld['ec']
    # Initialize empty lists.
    acpt_stars, rjct_stars = [], []

    for st_ind in acpt_indx:
        acpt_stars.append([ids[st_ind], x[st_ind], y[st_ind],
                          mags[st_ind], em[st_ind], cols[st_ind],
                          ec[st_ind]])

    for st_ind in rjct_indx:
        rjct_stars.append([ids[st_ind], x[st_ind], y[st_ind],
                          mags[st_ind], em[st_ind], cols[st_ind],
                          ec[st_ind]])

    return acpt_stars, rjct_stars


def main(cld, clp, pd):
    """
    Accept and reject stars in and out of the cluster's boundaries according to
    a given criteria based on their photometric errors.
    """

    # Unpack data.
    mode = pd['mode']
    mag, e_mag, e_col1 = cld['mags'], cld['em'], cld['ec']
    er_mode, e_max, be = pd['er_params'][:3]

    # Get value of brightest and dimmest stars.
    min_mag, max_mag = min(mag), max(mag)
    # Define max limit for the box that holds the brightest stars.
    be_m = (min_mag + be)
    # Create a segmented list in magnitude.
    # Magnitude range.
    delta_mag = max_mag - be_m
    # Width of the intervals in magnitude.
    interv_mag = 0.5
    # Number of intervals.
    n_interv = int(round(delta_mag / interv_mag))
    # Define list of points spanning the magnitude range starting from the
    # bright end. The '+ interv_mag' term is intentional so that the
    # err_medians function defines ranges around these values and they get
    # positioned in the middle of the magnitude interval.
    mag_value = [be_m + interv_mag * (q + interv_mag) for q in range(n_interv)]

    # Pack params to pass. These values are used by the 'eyefit' function and
    # more importantly the err_medians function which is called by the
    # synth_clust function.
    err_pck = [be_m, interv_mag, mag_value]

    # Flag indicates that the function had to fall back to the
    # 'e_max'-based rejection method since the selected one failed.
    err_max_fallback = False
    # Flag indicates that no error rejection was possible.
    err_all_fallback = False

    # Check selected mode.
    if mode in ('auto', 'semi'):
        # If 'semi' is set, check for the flag that indicates which method to
        # use and override the one in the input params file.
        if mode == 'semi':
            # Unpack semi flag
            err_flag_semi = pd['err_flag_semi']
            if err_flag_semi != 0:
                if err_flag_semi in {1, 2, 3}:
                    # Set error mode to use.
                    mode_map = {1: 'emax', 2: 'lowexp', 3: 'eyefit'}
                    er_mode = mode_map[int(err_flag_semi)]
                    # Save mode in 'semi_input' to variable.
                    pd['er_params'][0] = er_mode
                    print("Semi: using method selected: {}.".format(er_mode))
                else:
                    print("  WARNING: wrong error method in semi input file."
                          "\n  Falling back to emax.")
                    er_mode, pd['er_params'][0] = 'emax', 'emax'

        # Check which error rejection algorithm was selected in the input
        # file.
        if er_mode == 'emax':
            # Call function to reject stars with errors > e_max.
            acpt_indx, rjct_indx, err_plot = e_a_r_mx.main(cld, **pd)

        elif er_mode in ('lowexp', 'eyefit'):
            try:
                if er_mode == 'lowexp':
                    # Call N sigma exp function.
                    acpt_indx, rjct_indx, err_plot = e_a_r_le.main(
                        cld, be_m, **pd)

                elif er_mode == 'eyefit':
                    # Call 'eyefit' function.
                    acpt_indx, rjct_indx, err_plot = e_a_r_ef.main(
                        err_pck, cld, **pd)
            except RuntimeError:
                print("  WARNING: {} function could not be fitted.\n"
                      "  Falling back to e_max function.".format(er_mode))
                # Call function to reject stars with errors > e_max.
                acpt_indx, rjct_indx, err_plot = e_a_r_mx.main(cld, **pd)
                err_max_fallback = True

        # If list of accepted stars is empty, fall back to e_max limit.
        if not acpt_indx and er_mode != 'emax':
            # Call function to reject stars with errors > e_max.
            acpt_indx, rjct_indx, err_plot = e_a_r_mx.main(cld, **pd)
            err_max_fallback = True

            if acpt_indx:
                print("  WARNING: No stars accepted based on their errors.\n"
                      "  Using all stars with errors < {}".format(e_max))

            # If there's still no accepted stars, use all.
            else:
                print("  WARNING: No stars accepted based on their errors.\n"
                      "  Using all stars.")
                # Store all indexes.
                acpt_indx, rjct_indx = [i for i in range(len(mag))], []
                err_all_fallback = True

        # If the method used was e_max, use all stars.
        elif not acpt_indx and er_mode == 'emax':
            print("  WARNING: No stars accepted based on their errors.\n"
                  "  Using all stars.")
            # Store all indexes.
            acpt_indx, rjct_indx = [i for i in range(len(mag))], []
            err_plot = []
            err_all_fallback = True

        # Call function to store stars according to the returned indexes.
        acpt_stars, rjct_stars = err_sel_stars(acpt_indx, rjct_indx, cld)

    # If 'manual' mode is set, display errors distributions and ask the user
    # to accept it or else use all stars except those with errors > e_max in
    # either the magnitude or the color.
    elif mode == 'manual':
        move_on = False
        while not move_on:
            while True:
                try:
                    answer_rad = int(raw_input(
                        "Choose a method for error-based stars rejection:\n"
                        "  1 (emax), 2 (lowexp), 3 (eyefit),"
                        " 4 (use all stars): "))

                    if answer_rad == 1:
                        e_max_n = raw_input(
                            'Select maximum error value: ')
                        pd['er_params'][1] = float(e_max_n)
                        # Call function to reject stars with errors > e_max.
                        acpt_indx, rjct_indx, err_plot = e_a_r_mx.main(
                            e_mag, e_col1)
                        er_mode, pd['er_params'][0] = 'emax', 'emax'
                        break
                    elif answer_rad == 2:
                        N_sig = raw_input(
                            "Select number of sigmas to lift the curve: ")
                        pd['er_params'][4] = float(N_sig)
                        # Call N sigma exp function.
                        acpt_indx, rjct_indx, err_plot = e_a_r_le.main(
                            cld, be_m, **pd)
                        er_mode, pd['er_params'][0] = 'lowexp', 'lowexp'
                        break
                    elif answer_rad == 3:
                        # Call 'eyefit' function.
                        acpt_indx, rjct_indx, err_plot = e_a_r_ef.main(
                            err_pck, cld, **pd)
                        er_mode, pd['er_params'][0] = 'eyefit', 'eyefit'
                        break
                    elif answer_rad == 4:
                        # Store all indexes.
                        acpt_indx, rjct_indx = [i for i in range(len(mag))], []
                        err_plot, pd['er_params'][0] = [], ''
                        err_all_fallback = True
                        break

                    if answer_rad not in {1, 2, 3, 4}:
                        print("Sorry, input is not valid. Try again.")
                    else:
                        if len(acpt_indx) == 0:
                            print("No stars left after error rejection.\nTry"
                                  " again with a different method or use all "
                                  "stars (4).")
                        else:
                            break
                except:
                    print("Sorry, input is not valid. Try again.")

            # Call function to store stars according to the returned indexes.
            acpt_stars, rjct_stars = err_sel_stars(acpt_indx, rjct_indx,
                                                   cld)

            if answer_rad != 4:
                # Display automatic errors rejection.
                display_errors.main(er_mode, mag, err_plot, acpt_stars,
                                    rjct_stars, err_pck, **pd)
                plt.show()
                # Ask if keep or reject.
                while True:
                    try:
                        mv_on = raw_input('Accept error based rejection?'
                                          ' (y/n) ')
                        if mv_on == 'y':
                            print("Fit accepted.")
                            move_on = True
                            break
                        elif mv_on == 'n':
                            print("Fit rejected.")
                            break
                        else:
                            print("Sorry, input is not valid. Try again.")
                    except:
                        print("Sorry, input is not valid. Try again.")
            else:
                move_on = True

    print("Stars accepted/rejected based on their errors.")

    clp['err_flags'] = [err_all_fallback, err_max_fallback]
    clp['acpt_stars'], clp['rjct_stars'], clp['err_plot'],\
        clp['err_pck'] = acpt_stars, rjct_stars, err_plot, err_pck

    return clp, pd
