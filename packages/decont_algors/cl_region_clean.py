
import local_cell_clean


def nmemb_sel(n_memb, memb_prob_avrg_sort):
    '''
    Algorithm to select which stars to use by the best fit function.
    '''
    cl_reg_fit, cl_reg_no_fit, min_prob = memb_prob_avrg_sort, [], 0.

    # Check approximate number of true members obtained by the structural
    # analysis.
    if n_memb > 10:

        # Total number of stars in the cluster region.
        n_tot = len(memb_prob_avrg_sort)

        # If there are less stars in the cluster region than than n_memb stars,
        # use all stars in the cluster region.
        if n_memb >= n_tot:
            txt = "  WARNING: fewer stars in cluster region ({}) than\n" +\
                "  the number of approximate cluster members ({}).\n" +\
                "  No removal applied."
            print(txt.format(n_tot, n_memb))
            # Use all stars in the cluster region.
            indx, min_prob = n_tot, 0.
        else:
            # Use the first n_memb stars, ie: those stars with the highest
            # membership probability values.
            indx, min_prob = n_memb, zip(*memb_prob_avrg_sort)[-1][n_memb]

        cl_reg_fit, cl_reg_no_fit = memb_prob_avrg_sort[:indx], \
            memb_prob_avrg_sort[indx:]

    else:
        print("  WARNING: less than 10 stars identified as\n"
              "  cluster members. No removal applied.")

    return cl_reg_fit, cl_reg_no_fit, min_prob


def top_h(memb_prob_avrg_sort):
    '''
    Reject stars in the lower half of the membership probabilities list.
    '''

    cl_reg_fit, cl_reg_no_fit, min_prob = memb_prob_avrg_sort, [], 0.

    middle_indx = int(len(memb_prob_avrg_sort) / 2)
    rem_fit = memb_prob_avrg_sort[:middle_indx]
    # Check number of stars left.
    if len(rem_fit) > 10:
        cl_reg_fit, cl_reg_no_fit, min_prob = rem_fit,\
            memb_prob_avrg_sort[middle_indx:],\
            memb_prob_avrg_sort[middle_indx][-1]
    else:
        print("  WARNING: less than 10 stars left after reducing\n"
              "  by top half membership probability. No removal applied.")

    return cl_reg_fit, cl_reg_no_fit, min_prob


def manual(memb_prob_avrg_sort, fld_clean_prob, min_prob_i=None):
    '''
    Find index of star with membership probability < min_prob_i.
    '''
    cl_reg_fit, cl_reg_no_fit, min_prob = memb_prob_avrg_sort, [], 0.

    if min_prob_i is None:
        # Manual mode.
        min_prob_man = fld_clean_prob
    else:
        # MP>=0.5 mode.
        min_prob_man = min_prob_i

    indx = 0
    for star in memb_prob_avrg_sort:
        if star[-1] < min_prob_man:
            break
        else:
            indx += 1

    if len(memb_prob_avrg_sort[:indx]) > 10:
        cl_reg_fit, cl_reg_no_fit, min_prob = \
            memb_prob_avrg_sort[:indx], memb_prob_avrg_sort[indx:],\
            min_prob_man
    else:
        if min_prob_i is None:
            print("  WARNING: fewer than 10 stars left after removing\n"
                  "  stars with MP<={:.2f}. No removal applied.".format(
                      min_prob_man))
        else:
            print("  WARNING: fewer than 10 stars left after removing\n"
                  "  stars with MP<=0.5. No removal applied.")

    return cl_reg_fit, cl_reg_no_fit, min_prob


def main(clp, fld_clean_mode, fld_clean_bin, fld_clean_prob, **kwargs):
    '''
    Remove stars from cluster region, according to a given membership
    probability lower limit, minimum magnitude limit or local density-based
    removal.
    '''

    # Default assignment.
    cl_reg_fit, cl_reg_no_fit, min_prob, bin_edges =\
        clp['memb_prob_avrg_sort'], [], 0., 0.

    if fld_clean_mode == 'all':
        # Skip reduction process.
        print('No field star removal applied on cluster region.')

    # If the DA was skipped and any method but 'local' or 'mag' is selected,
    # don't run.
    elif clp['flag_decont_skip'] and fld_clean_mode not in ('local', 'mag'):
        print("  WARNING: decontamination algorithm was skipped.\n"
              "  Can't apply '{}' field stars removal method.".format(
                  fld_clean_mode))

    # If no field regions were defined, this mode won't work.
    elif clp['flag_no_fl_regs_c'] and fld_clean_mode == 'local':
        print("  WARNING: no field regions were defined. Can't apply\n"
              "  '{}' field stars removal method.".format(fld_clean_mode))

    else:
        # This mode works if the DA did not run but it needs field regions
        # defined.
        if fld_clean_mode == 'local':
            cl_reg_fit, cl_reg_no_fit, min_prob, bin_edges =\
                local_cell_clean.main(
                    clp['field_regions_c'], clp['memb_prob_avrg_sort'],
                    clp['flag_decont_skip'], fld_clean_bin)

        if fld_clean_mode == 'n_memb':
            cl_reg_fit, cl_reg_no_fit, min_prob = nmemb_sel(
                clp['n_memb'], clp['memb_prob_avrg_sort'])

        if fld_clean_mode == 'mp_05':
            cl_reg_fit, cl_reg_no_fit, min_prob = manual(
                clp['memb_prob_avrg_sort'], fld_clean_prob, 0.5)

        elif fld_clean_mode == 'top_h':
            cl_reg_fit, cl_reg_no_fit, min_prob = top_h(
                clp['memb_prob_avrg_sort'])

        elif fld_clean_mode == 'man':
            cl_reg_fit, cl_reg_no_fit, min_prob = manual(
                clp['memb_prob_avrg_sort'], fld_clean_prob)

        print("Stars selected in cluster region, '{}' method ({}).".format(
            fld_clean_mode, len(cl_reg_fit)))

    clp['cl_reg_fit'], clp['cl_reg_no_fit'], clp['cl_reg_clean_plot'] =\
        cl_reg_fit, cl_reg_no_fit, [min_prob, bin_edges]
    return clp
