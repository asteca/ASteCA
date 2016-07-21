
from ..inp import input_params as g
import local_cell_clean


def nmemb_sel(n_memb, memb_prob_avrg_sort):
    '''
    Algorithm to select which stars to use by the best fit function.
    Will set the minimum probability value such that an equal number of
    stars are used in the best fit process, as the approximate number of
    members found when comparing the density of the cluster region with that
    of the field regions defined.
    '''
    rem_memb_fit, rem_memb_no_fit, rem_plot_pars = memb_prob_avrg_sort, [], \
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
            indx, rem_plot_pars = n_tot, [0.]
        else:
            # Use the first n_memb stars, ie: those stars with the highest
            # membership probability values.
            indx, rem_plot_pars = n_memb, \
                [zip(*memb_prob_avrg_sort)[-1][n_memb]]

        rem_memb_fit, rem_memb_no_fit = memb_prob_avrg_sort[:indx], \
            memb_prob_avrg_sort[indx:]

    else:
        print("  WARNING: less than 10 stars identified as true\n"
              "  cluster members. Using full list.")

    return rem_memb_fit, rem_memb_no_fit, rem_plot_pars


def top_h(memb_prob_avrg_sort):
    '''
    Reject stars in the lower half of the membership probabilities list.
    '''

    rem_memb_fit, rem_memb_no_fit, rem_plot_pars = memb_prob_avrg_sort, [], \
        [0.]

    middle_indx = int(len(memb_prob_avrg_sort) / 2)
    rem_fit = memb_prob_avrg_sort[:middle_indx]
    # Check number of stars left.
    if len(rem_fit) > 10:
        rem_memb_fit, rem_memb_no_fit, rem_plot_pars = rem_fit, \
            memb_prob_avrg_sort[middle_indx:], \
            [memb_prob_avrg_sort[middle_indx][-1]]
    else:
        print("  WARNING: less than 10 stars left after reducing\n"
              "  by top half membership probability. Using full list.")

    return rem_memb_fit, rem_memb_no_fit, rem_plot_pars


def manual(memb_prob_avrg_sort, min_prob=None):
    '''
    Find index of star with membership probability < min_prob.
    '''
    rem_memb_fit, rem_memb_no_fit, rem_plot_pars = memb_prob_avrg_sort, [], \
        [0.]

    if min_prob is None:
        # Manual mode.
        min_prob_man = g.rm_params[2]
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
        rem_memb_fit, rem_memb_no_fit, rem_plot_pars = \
            memb_prob_avrg_sort[:indx], memb_prob_avrg_sort[indx:],\
            [min_prob_man]
    else:
        if min_prob is None:
            print("  WARNING: less than 10 stars left after reducing\n"
                  "  by manual membership probability. Using full list.")
        else:
            print("  WARNING: less than 10 stars left after reducing\n"
                  "  by MP>=0.5 selection. Using full list.")

    return rem_memb_fit, rem_memb_no_fit, rem_plot_pars


def man_mag(memb_prob_avrg_sort):
    '''
    Reject stars beyond the given magnitude limit.
    '''

    min_prob_man = g.rm_params[2]
    rem_memb_fit, rem_memb_no_fit, rem_plot_pars = memb_prob_avrg_sort, [],\
        [0.]

    rem_fit, rem_not_fit = [], []
    for star in memb_prob_avrg_sort:
        if star[3] <= min_prob_man:
            rem_fit.append(star)
        else:
            rem_not_fit.append(star)

    # Check number of stars left.
    if len(rem_fit) > 10:
        rem_memb_fit, rem_memb_no_fit, rem_plot_pars = rem_fit, rem_not_fit,\
            [min_prob_man]
    else:
        print("  WARNING: less than 10 stars left after reducing\n"
              "  by magnitude limit. Using full list.")

    return rem_memb_fit, rem_memb_no_fit, rem_plot_pars


def main(n_memb, flag_no_fl_regs, bayes_da_return, field_region):
    '''
    Remove stars from cluster region, according to a given membership
    probability lower limit, minimum magnitude limit or local density-based
    removal.
    '''

    memb_prob_avrg_sort, flag_decont_skip = bayes_da_return
    mode_rem_memb = g.rm_params[0]

    # Default assignment.
    rem_memb_fit, rem_memb_no_fit, rem_plot_pars = memb_prob_avrg_sort, [], \
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
            rem_memb_fit, rem_memb_no_fit, rem_plot_pars = \
                local_cell_clean.main(bayes_da_return, field_region)

        if mode_rem_memb == 'n_memb':
            rem_memb_fit, rem_memb_no_fit, rem_plot_pars = nmemb_sel(
                n_memb, memb_prob_avrg_sort)

        if mode_rem_memb == 'mp_05':
            rem_memb_fit, rem_memb_no_fit, rem_plot_pars = \
                manual(memb_prob_avrg_sort, 0.5)

        elif mode_rem_memb == 'top_h':
            rem_memb_fit, rem_memb_no_fit, rem_plot_pars = \
                top_h(memb_prob_avrg_sort)

        elif mode_rem_memb == 'man':
            rem_memb_fit, rem_memb_no_fit, rem_plot_pars = \
                manual(memb_prob_avrg_sort)

        elif mode_rem_memb == 'mag':
            rem_memb_fit, rem_memb_no_fit, rem_plot_pars = \
                man_mag(memb_prob_avrg_sort)

        print('Membership-based removal function applied.')

    return rem_memb_fit, rem_memb_no_fit, rem_plot_pars
