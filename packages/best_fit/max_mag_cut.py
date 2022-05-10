
import numpy as np
import copy


def main(cl_reg_fit, Max_mag):
    """
    Reject stars beyond the maximum magnitude given.
    """

    # Maximum observed (main) magnitude.
    max_mag_obs = np.max(list(zip(*list(zip(*cl_reg_fit))[1:][2]))[0])

    if Max_mag == 'max':
        # No magnitude cut applied.
        cl_max_mag, max_mag_syn = copy.deepcopy(cl_reg_fit), max_mag_obs
    else:
        Max_mag = float(Max_mag)
        star_lst = []
        for star in cl_reg_fit:
            # Check main magnitude value.
            if star[3][0] <= Max_mag:
                # Keep stars brighter that the magnitude limit.
                star_lst.append(star)

        # Check number of stars left.
        if len(star_lst) > 10:
            # For the synthetic clusters, use the minimum value between the
            # selected 'Max_mag' value and the maximum observed magnitude.
            # This prevents large 'Max_mag' values from generating synthetic
            # clusters with low mass stars in the not-observed region.
            cl_max_mag, max_mag_syn = star_lst, min(Max_mag, max_mag_obs)
            print("Maximum magnitude cut applied ({:.1f} mag)".format(
                max_mag_syn))
        else:
            cl_max_mag, max_mag_syn = copy.deepcopy(cl_reg_fit), max_mag_obs
            print("  WARNING: less than 10 stars left after removing\n"
                  "  stars by magnitude limit. No removal applied.")

    return cl_max_mag, max_mag_syn
