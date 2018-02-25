
import numpy as np


def main(st_dist, M_total):
    """
    http://www.astro.ru.nl/~slarsen/teaching/Galaxies/cmd.pdf
    http://python4mpia.github.io/fitting_data/MC-sampling-from-Salpeter.html

    Returns a mass distribution according to a given IMF and a total cluster
    mass.

    Generate  N_stars for each interval (m, m+dm) with masses
    randomly distributed within the limiting masses that delimit the
    interval, ie: (m, m+dm).

    """
    if st_dist[M_total][-1] is True:
        base, scale, N_stars = st_dist[M_total][:-1]
        mass_dist = np.random.random(N_stars) * scale + base
    else:
        mass_dist = st_dist[M_total]

    return mass_dist
