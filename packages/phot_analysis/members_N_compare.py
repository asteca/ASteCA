

def main(n_memb, bayes_da_return):
    '''
    Compare the estimated number of true members obtained via the stars density
    analysis done in `get_members_number` with the number of stars in the
    cluster region that are assigned a MP of 0.5 or more. These stars are the
    ones with a greater probability of being cluster members than field region
    stars.
    '''

    memb_prob_avrg_sort, flag_decont_skip = bayes_da_return

    memb_par, n_memb_da, flag_memb_par = float("inf"), -1., False
    # Obtain parameter if the DA was applied.
    if not flag_decont_skip:

        n_memb_da = 0
        # Number of stars assigned a MP>=0.5.
        for star in memb_prob_avrg_sort:
            if star[7] >= 0.5:
                n_memb_da += 1

        # Obtain only if n_memb was calculated.
        if n_memb != -1.:

            if n_memb != 0 or n_memb_da != 0:
                # Obtain parameter.
                memb_par = (float(n_memb) - float(n_memb_da)) / \
                    (float(n_memb) + float(n_memb_da))

            # Set flag.
            if abs(memb_par) > 0.33:
                flag_memb_par = True
                print("  WARNING: the number of estimated true cluster"
                      " members\n  (structural vs. photometric) differ"
                      " greatly.")

    return memb_par, n_memb_da, flag_memb_par
