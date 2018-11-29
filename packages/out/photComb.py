
from random import shuffle


def main(pd, clp):
    """
    Combine the photometry of the accepted/rejected stars in the cluster
    and field regions defined. Used to define the limits for the CMDs used
    in the following plots.
    """

    # Stars in complete photometry, error accepted, cluster region.
    cl_ac_col_0 = list(list(zip(*list(zip(*clp['cl_region_c']))[5]))[0])
    cl_ac_mag_0 = list(list(zip(*list(zip(*clp['cl_region_c']))[3]))[0])
    cl_ac_col_1 = []
    if len(pd['colors']) > 1:
        cl_ac_col_1 = list(list(zip(*list(zip(*clp['cl_region_c']))[5]))[1])

    # Stars in complete photometry, error rejected, cluster region.
    cl_rj_col_0, cl_rj_mag_0, cl_rj_col_1 = [], [], []
    if len(clp['cl_region_rjct_c']) > 0:
        cl_rj_col_0 = list(list(zip(*list(zip(
            *clp['cl_region_rjct_c']))[5]))[0])
        cl_rj_mag_0 = list(list(zip(*list(zip(
            *clp['cl_region_rjct_c']))[3]))[0])
        if len(pd['colors']) > 1:
            cl_rj_col_1 = list(list(zip(*list(zip(
                *clp['cl_region_rjct_c']))[5]))[1])

    # Only use 25$ of the rejected stars. This way they have a lesser impact
    # on the CMD limits.
    Nh = int(len(cl_rj_mag_0) / 4.)
    shuffle(cl_rj_mag_0)
    shuffle(cl_rj_col_0)
    shuffle(cl_rj_col_1)

    # Combine all data into a single list for each dimension.
    mag_0_comb = cl_ac_mag_0 + cl_rj_mag_0[:Nh]
    col_0_comb = cl_ac_col_0 + cl_rj_col_0[:Nh]
    col_1_comb = cl_ac_col_1 + cl_rj_col_1[:Nh]

    stars_f_rjct, stars_f_acpt = field_region_stars(
        clp['field_regions_c'], clp['field_regions_rjct_c'], pd['colors'])

    clp.update({
        'stars_f_rjct': stars_f_rjct, 'stars_f_acpt': stars_f_acpt,
        'mag_0_comb': mag_0_comb, 'col_0_comb': col_0_comb,
        'col_1_comb': col_1_comb})
    return clp


def field_region_stars(field_regions, field_regions_rjct, colors):
    """
    Generate list with accepted/rejected stars within all the defined field
    regions.
    """
    stars_f_acpt = [[], [], []]
    if field_regions:
        # Extract color(s) and main magnitude defined.
        stars_f_acpt[0] = [
            star[3][0] for flrg in field_regions for star in flrg]
        stars_f_acpt[1] = [
            star[5][0] for flrg in field_regions for star in flrg]
        if len(colors) > 1:
            stars_f_acpt[2] = [
                star[5][1] for flrg in field_regions for star in flrg]

    stars_f_rjct = [[], [], []]
    if field_regions_rjct:
        # Extract color(s) and main magnitude defined.
        stars_f_rjct[0] = [
            star[3][0] for flrg in field_regions_rjct for star in flrg]
        stars_f_rjct[1] = [
            star[5][0] for flrg in field_regions_rjct for star in flrg]
        if len(colors) > 1:
            stars_f_rjct[2] = [
                star[5][1] for flrg in field_regions_rjct for star in flrg]

    return stars_f_rjct, stars_f_acpt
