
import local_cell_clean


def nmemb_sel(n_memb, memb_prob_avrg_sort):
    '''
    Algorithm to select which stars to use by the best fit function.
    '''
    cl_reg_fit, cl_reg_no_fit, cl_reg_clean_plot = memb_prob_avrg_sort, [], \
        [0.]

    # Check approximate number of true members obtained by the structural
    # analysis.
    if n_memb > 10:

        # Total number of stars in the cluster region.
        n_tot = len(memb_prob_avrg_sort)

        # If there are less stars in the cluster region than than n_memb stars,
        # use all stars in the cluster region.
        if n_memb >= n_tot:
            # Use all stars in the cluster region.
            indx, cl_reg_clean_plot = n_tot, [0.]
        else:
            # Use the first n_memb stars, ie: those stars with the highest
            # membership probability values.
            indx, cl_reg_clean_plot = n_memb, \
                [zip(*memb_prob_avrg_sort)[-1][n_memb]]

        cl_reg_fit, cl_reg_no_fit = memb_prob_avrg_sort[:indx], \
            memb_prob_avrg_sort[indx:]

    else:
        print("  WARNING: less than 10 stars identified as true\n"
              "  cluster members. Using full list.")

    return cl_reg_fit, cl_reg_no_fit, cl_reg_clean_plot


def top_h(memb_prob_avrg_sort):
    '''
    Reject stars in the lower half of the membership probabilities list.
    '''

    cl_reg_fit, cl_reg_no_fit, cl_reg_clean_plot = memb_prob_avrg_sort, [], \
        [0.]

    middle_indx = int(len(memb_prob_avrg_sort) / 2)
    rem_fit = memb_prob_avrg_sort[:middle_indx]
    # Check number of stars left.
    if len(rem_fit) > 10:
        cl_reg_fit, cl_reg_no_fit, cl_reg_clean_plot = rem_fit, \
            memb_prob_avrg_sort[middle_indx:], \
            [memb_prob_avrg_sort[middle_indx][-1]]
    else:
        print("  WARNING: less than 10 stars left after reducing\n"
              "  by top half membership probability. Using full list.")

    return cl_reg_fit, cl_reg_no_fit, cl_reg_clean_plot


def manual(memb_prob_avrg_sort, rm_params, min_prob=None):
    '''
    Find index of star with membership probability < min_prob.
    '''
    cl_reg_fit, cl_reg_no_fit, cl_reg_clean_plot = memb_prob_avrg_sort, [], \
        [0.]

    if min_prob is None:
        # Manual mode.
        min_prob_man = rm_params[2]
    else:
        # MP>=0.5 mode.
        min_prob_man = min_prob

    indx = 0
    for star in memb_prob_avrg_sort:
        if star[-1] < min_prob_man:
            break
        else:
            indx += 1

    if len(memb_prob_avrg_sort[:indx]) > 10:
        cl_reg_fit, cl_reg_no_fit, cl_reg_clean_plot = \
            memb_prob_avrg_sort[:indx], memb_prob_avrg_sort[indx:],\
            [min_prob_man]
    else:
        if min_prob is None:
            print("  WARNING: less than 10 stars left after reducing\n"
                  "  by manual membership probability. Using full list.")
        else:
            print("  WARNING: less than 10 stars left after reducing\n"
                  "  by MP>=0.5 selection. Using full list.")

    return cl_reg_fit, cl_reg_no_fit, cl_reg_clean_plot


def man_mag(memb_prob_avrg_sort, rm_params):
    '''
    Reject stars beyond the given magnitude limit.
    '''

    min_prob_man = rm_params[2]
    cl_reg_fit, cl_reg_no_fit, cl_reg_clean_plot = memb_prob_avrg_sort, [],\
        [0.]

    rem_fit, rem_not_fit = [], []
    for star in memb_prob_avrg_sort:
        if star[3] <= min_prob_man:
            rem_fit.append(star)
        else:
            rem_not_fit.append(star)

    # Check number of stars left.
    if len(rem_fit) > 10:
        cl_reg_fit, cl_reg_no_fit, cl_reg_clean_plot = rem_fit, rem_not_fit,\
            [min_prob_man]
    else:
        print("  WARNING: less than 10 stars left after reducing\n"
              "  by magnitude limit. Using full list.")

    return cl_reg_fit, cl_reg_no_fit, cl_reg_clean_plot


def main(clp, rm_params, **kwargs):
    '''
    Remove stars from cluster region, according to a given membership
    probability lower limit, minimum magnitude limit or local density-based
    removal.
    '''

    n_memb, flag_no_fl_regs, field_regions, memb_prob_avrg_sort,\
        flag_decont_skip = [
            clp[_] for _ in ['n_memb', 'flag_no_fl_regs', 'field_regions',
                             'memb_prob_avrg_sort', 'flag_decont_skip']]
    mode_rem_memb = rm_params[0]

    # Default assignment.
    cl_reg_fit, cl_reg_no_fit, cl_reg_clean_plot = memb_prob_avrg_sort, [], \
        [0.]

    if mode_rem_memb == 'skip':
        # Skip reduction process.
        print('Membership-based removal function skipped.')

    # If the DA was skipped and any method but 'local' or 'mag' is selected,
    # don't run.
    elif flag_decont_skip and mode_rem_memb not in ['local', 'mag']:
        print("  WARNING: decontamination algorithm was skipped.\n"
              "  Can't apply '{}' MP removal method.\n"
              "  Using full list.").format(mode_rem_memb)

    # If no field regions were defined, this mode won't work.
    elif flag_no_fl_regs and mode_rem_memb == 'local':
        print("  WARNING: no field regions were defined. Can't apply\n"
              "  '{}' MP removal method. Using full list.").format(
                  mode_rem_memb)

    else:
        # This mode works if the DA did not run but it needs field regions
        # defined.
        if mode_rem_memb == 'local':
            cl_reg_fit, cl_reg_no_fit, cl_reg_clean_plot = \
                local_cell_clean.main(field_regions, memb_prob_avrg_sort,
                                      flag_decont_skip, rm_params)

        if mode_rem_memb == 'n_memb':
            cl_reg_fit, cl_reg_no_fit, cl_reg_clean_plot = nmemb_sel(
                n_memb, memb_prob_avrg_sort)

        if mode_rem_memb == 'mp_05':
            cl_reg_fit, cl_reg_no_fit, cl_reg_clean_plot = \
                manual(memb_prob_avrg_sort, rm_params, 0.5)

        elif mode_rem_memb == 'top_h':
            cl_reg_fit, cl_reg_no_fit, cl_reg_clean_plot = \
                top_h(memb_prob_avrg_sort)

        elif mode_rem_memb == 'man':
            cl_reg_fit, cl_reg_no_fit, cl_reg_clean_plot = \
                manual(memb_prob_avrg_sort, rm_params)

        elif mode_rem_memb == 'mag':
            cl_reg_fit, cl_reg_no_fit, cl_reg_clean_plot = \
                man_mag(memb_prob_avrg_sort)

        print('Membership-based removal function applied.')

    clp['cl_reg_fit'], clp['cl_reg_no_fit'], clp['cl_reg_clean_plot'] =\
        cl_reg_fit, cl_reg_no_fit, cl_reg_clean_plot
    return clp
