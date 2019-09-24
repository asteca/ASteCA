
import numpy as np


def main(st_dist_mass, M_total):
    """
    http://www.astro.ru.nl/~slarsen/teaching/Galaxies/cmd.pdf
    http://python4mpia.github.io/fitting_data/MC-sampling-from-Salpeter.html

    Returns a mass distribution according to a given IMF and a total cluster
    mass.

    Generate N_stars for each interval (m, m+dm) with masses randomly
    distributed within this interval.

    """

    # This is not in use since May 2019 (see 'imf.py'), because all the
    # mass distributions need to be equal. It could be useful when #239 is
    # implemented.
    #
    # if st_dist_mass[M_total][-1] is True:
    #     base, scale, N_stars = st_dist_mass[M_total][:-1]
    #     mass_dist = np.random.random(N_stars) * scale + base
    # else:
    # mass_dist = st_dist_mass[M_total]
    mass_dist = st_dist_mass[0][:np.searchsorted(st_dist_mass[1], M_total)]

    return mass_dist
