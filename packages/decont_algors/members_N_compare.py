

def main(clp):
    '''
    Compare the estimated number of true members obtained via the stars density
    analysis done in `members_number` with the number of stars in the
    cluster region that are assigned a MP of 0.5 or more. These stars are the
    ones with a greater probability of being cluster members than field region
    stars.
    '''
    memb_par, n_memb_da, flag_memb_par = float("inf"), -1., False
    # Obtain parameter if the DA was applied.
    if not clp['flag_decont_skip']:

        n_memb_da = 0
        # Number of stars assigned a MP>=0.5.
        for star in clp['memb_prob_avrg_sort']:
            if star[7] >= 0.5:
                n_memb_da += 1

        # Obtain only if n_memb was calculated.
        if not clp['flag_num_memb_low']:

            if clp['n_memb'] != 0 or n_memb_da != 0:
                # Obtain parameter.
                memb_par = (float(clp['n_memb']) - float(n_memb_da)) / \
                    (float(clp['n_memb']) + float(n_memb_da))
            else:
                print("  WARNING: could not obtain 'memb_par' parameter.")

            # Set flag.
            if abs(memb_par) > 0.33:
                flag_memb_par = True
                print("  WARNING: structural vs. photometric true cluster\n"
                      "  members estimated differ by a factor > 2.")

    clp['memb_par'], clp['n_memb_da'], clp['flag_memb_par'] =\
        memb_par, n_memb_da, flag_memb_par
    return clp
