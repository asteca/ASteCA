
import matplotlib.pyplot as plt
import display_errors
import err_accpt_rejct_max as e_a_r_mx


def err_sel_stars(acpt_indx, cld):
    '''
    Go through the lists of indexes to select accepted/rejected stars.
    '''
    # Unpack data.
    ids, x, y, mags, em, cols, ec = cld['ids'], cld['x'], cld['y'],\
        cld['mags'], cld['em'], cld['cols'], cld['ec']
    # Initialize empty lists.
    acpt_stars, rjct_stars = [], []
    # For every star in the cluster's file.
    for i, st_idx in enumerate(ids):
        if i in acpt_indx:
            # Separate magnitudes and their errors.
            m_a, e_ma = [], []
            for j, m in enumerate(mags):
                m_a.append(m[i])
                e_ma.append(em[j][i])
            # Separate colors and their errors.
            c_a, e_ca = [], []
            for j, c in enumerate(cols):
                c_a.append(c[i])
                e_ca.append(ec[j][i])
            # Store all data on accepted star in list.
            acpt_stars.append([st_idx, x[i], y[i], m_a, e_ma, c_a, e_ca])
        else:
            # Separate magnitudes and their errors.
            m_r, e_mr = [], []
            for j, m in enumerate(mags):
                m_r.append(m[i])
                e_mr.append(em[j][i])
            # Separate colors and their errors.
            c_r, e_cr = [], []
            for j, c in enumerate(cols):
                c_r.append(c[i])
                e_cr.append(ec[j][i])
            # Store all data on accepted star in list.
            rjct_stars.append([st_idx, x[i], y[i], m_r, e_mr, c_r, e_cr])

    return acpt_stars, rjct_stars


def main(cld, clp, pd):
    """
    Accept and reject stars in and out of the cluster's boundaries according to
    a given criteria based on their photometric errors.
    """
    # Use main magnitude.
    mmag = cld['mags'][0]
    # Flag indicates that no error rejection was possible.
    err_all_fallback = False

    if pd['run_mode'] in ('auto', 'semi'):
        e_max_man = pd['err_max']
        # Call function to reject stars with errors > e_max.
        acpt_indx = e_a_r_mx.main(cld, e_max_man)
        if not acpt_indx:
            print("  WARNING: No stars accepted based on their errors.\n"
                  "  Using all stars.")
            # Store all indexes.
            acpt_indx, err_all_fallback = range(len(mmag)), True

        # Call function to store stars according to the returned indexes.
        acpt_stars, rjct_stars = err_sel_stars(acpt_indx, cld)

    # If 'manual' mode is set, display errors distributions and ask the user
    # to accept it or else use all stars except those with errors > e_max in
    # either the magnitude or the color.
    elif pd['run_mode'] == 'manual':
        move_on = False
        while not move_on:
            while True:
                answer_rad = int(raw_input(
                    "Choose a method for error-based stars rejection:\n"
                    "  1 (emax), 2 (use all stars): "))

                if answer_rad == 1:
                    e_max_man = float(raw_input(
                        'Select maximum error value: '))
                    # Call function to reject stars with errors > e_max.
                    acpt_indx = e_a_r_mx.main(cld, e_max_man)
                    break
                elif answer_rad == 2:
                    # Store all indexes.
                    acpt_indx, e_max_man = range(len(mmag)), 'nan'
                    break
                if answer_rad not in (1, 2):
                    print("Sorry, input is not valid. Try again.")
                else:
                    if len(acpt_indx) == 0:
                        print("No stars left after error rejection. Try"
                              " again.")
                    else:
                        break

            # Call function to store stars according to the returned indexes.
            acpt_stars, rjct_stars = err_sel_stars(acpt_indx, cld)

            if answer_rad == 1:
                # Display automatic errors rejection.
                display_errors.main(
                    pd['filters'], pd['colors'], mmag, acpt_stars,
                    rjct_stars, e_max_man)
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
                    except Exception:
                        print("Sorry, input is not valid. Try again.")
            else:
                move_on = True

    print("Stars rejected based on their errors ({}).".format(len(rjct_stars)))

    clp['err_all_fallback'], clp['acpt_stars'], clp['rjct_stars'],\
        clp['err_max'] = err_all_fallback, acpt_stars, rjct_stars, e_max_man

    return clp, pd
