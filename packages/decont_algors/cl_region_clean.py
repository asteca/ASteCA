
import numpy as np
from . import local_cell_clean
from .decont_algors import sort_members


def main(clp, fld_clean_mode, fld_clean_bin, fld_clean_prob, **kwargs):
    """
    Remove stars from cluster region, according to a given membership
    probability lower limit, minimum magnitude limit or local density-based
    removal.
    """

    # Default assignment.
    cl_reg_fit, cl_reg_no_fit, local_rm_edges = clp['memb_prob_avrg_sort'],\
        [], None

    if fld_clean_mode == 'all':
        # Skip reduction process.
        print("No field star removal applied on cluster region")
        # Remove stars with a color far away from the x median.
        cl_reg_fit, cl_reg_no_fit = rmColorOutliers(cl_reg_fit, cl_reg_no_fit)

    # If the DA was skipped and any method but 'local' or 'mag' is selected,
    # don't run.
    elif clp['flag_decont_skip'] and fld_clean_mode not in ('local', 'mag'):
        print("  WARNING: decontamination algorithm was skipped\n"
              "  Can't apply '{}' field stars removal method".format(
                  fld_clean_mode))

    # If no field regions were defined, this mode won't work.
    elif clp['flag_no_fl_regs'] and fld_clean_mode == 'local':
        print("  WARNING: no field regions were defined. Can't apply\n"
              "  '{}' field stars removal method".format(fld_clean_mode))

    else:
        # This mode works if the DA did not run but it needs field regions
        # defined.
        if fld_clean_mode == 'local':
            cl_reg_fit, cl_reg_no_fit, local_rm_edges =\
                local_cell_clean.main(
                    clp['n_memb'], clp['field_regions'],
                    clp['memb_prob_avrg_sort'], clp['flag_decont_skip'],
                    fld_clean_bin)

        if fld_clean_mode == 'n_memb':
            cl_reg_fit, cl_reg_no_fit = nmemb_sel(
                clp['n_memb'], clp['memb_prob_avrg_sort'])

        if fld_clean_mode == 'mp_05':
            cl_reg_fit, cl_reg_no_fit = manual(
                clp['memb_prob_avrg_sort'], fld_clean_prob, 0.5)

        elif fld_clean_mode == 'top_h':
            cl_reg_fit, cl_reg_no_fit = top_h(
                clp['memb_prob_avrg_sort'])

        elif fld_clean_mode == 'man':
            cl_reg_fit, cl_reg_no_fit = manual(
                clp['memb_prob_avrg_sort'], fld_clean_prob)

        # Remove stars with a color far away from the x median.
        cl_reg_fit, cl_reg_no_fit = rmColorOutliers(cl_reg_fit, cl_reg_no_fit)

        # Sort according to the largest MPs.
        cl_reg_fit = sort_members(cl_reg_fit)
        cl_reg_no_fit = sort_members(cl_reg_no_fit)

        print("Stars selected in cluster region, '{}' method ({})".format(
            fld_clean_mode, len(cl_reg_fit)))

    clp['cl_reg_fit'], clp['cl_reg_no_fit'], clp['local_rm_edges'] =\
        cl_reg_fit, cl_reg_no_fit, local_rm_edges
    return clp


def nmemb_sel(n_memb, memb_prob_avrg_sort):
    """
    Algorithm to select which stars to use by the best fit function.
    """
    cl_reg_fit, cl_reg_no_fit = memb_prob_avrg_sort, []

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
            indx = n_tot
        else:
            # Use the first n_memb stars, ie: those stars with the highest
            # membership probability values.
            indx = n_memb

        cl_reg_fit, cl_reg_no_fit = memb_prob_avrg_sort[:indx], \
            memb_prob_avrg_sort[indx:]

    else:
        print("  WARNING: less than 10 stars identified as\n"
              "  cluster members. No removal applied.")

    return cl_reg_fit, cl_reg_no_fit


def top_h(memb_prob_avrg_sort):
    """
    Reject stars in the lower half of the membership probabilities list.
    """

    cl_reg_fit, cl_reg_no_fit = memb_prob_avrg_sort, []

    middle_indx = int(len(memb_prob_avrg_sort) / 2.)
    rem_fit = memb_prob_avrg_sort[:middle_indx]
    # Check number of stars left.
    if len(rem_fit) > 10:
        cl_reg_fit, cl_reg_no_fit = rem_fit, memb_prob_avrg_sort[middle_indx:]
    else:
        print("  WARNING: less than 10 stars left after reducing\n"
              "  by top half membership probability. No removal applied.")

    return cl_reg_fit, cl_reg_no_fit


def manual(memb_prob_avrg_sort, fld_clean_prob, min_prob_i=None):
    """
    Find index of star with membership probability < min_prob_i.
    """
    cl_reg_fit, cl_reg_no_fit = memb_prob_avrg_sort, [],

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
        cl_reg_fit, cl_reg_no_fit = \
            memb_prob_avrg_sort[:indx], memb_prob_avrg_sort[indx:]
    else:
        if min_prob_i is None:
            print("  WARNING: fewer than 10 stars left after removing\n"
                  "  stars with MP<={:.2f}. No removal applied.".format(
                      min_prob_man))
        else:
            print("  WARNING: fewer than 10 stars left after removing\n"
                  "  stars with MP<=0.5. No removal applied.")

    return cl_reg_fit, cl_reg_no_fit


def rmColorOutliers(cl_reg_fit, cl_reg_no_fit, Nstd=5.):
    """
    Reject stars that are far away from the sequence.
    """
    mag = np.array(list(zip(*cl_reg_fit))[3]).flatten()
    colors = np.array(list(zip(*cl_reg_fit))[5])

    # Median, std for colors
    col_med, col_std = np.median(colors), np.std(colors)
    # Median, std for magnitude.
    mag_med, mag_std = np.median(mag), np.std(mag)

    cl_reg_fit2, cl_reg_no_fit2 = [], []
    for i, st in enumerate(cl_reg_fit):
        # For the brightest stars
        if st[3] <= mag_med - .5 * mag_std:
            if np.any(st[5] < col_med - Nstd * col_std) or np.any(
                    st[5] > col_med + Nstd * col_std):
                cl_reg_no_fit2.append(st)
            else:
                cl_reg_fit2.append(st)
        else:
            # Make the left side a bit tighter than the right
            if np.any(st[5] < col_med - (Nstd - 2) * col_std) or np.any(
                    st[5] > col_med + (Nstd - 1) * col_std):
                cl_reg_no_fit2.append(st)
            else:
                cl_reg_fit2.append(st)

    # Add rejected stars to the old list.
    cl_reg_no_fit2 = cl_reg_no_fit + cl_reg_no_fit2

    if len(cl_reg_fit) > len(cl_reg_fit2):
        print("Removed {} CMD outlier stars".format(
            len(cl_reg_fit) - len(cl_reg_fit2)))

    return cl_reg_fit2, cl_reg_no_fit2
