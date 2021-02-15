
import numpy as np


def main(st_dist_mass, mean_bin_mr, bin_frac, M_total):
    """
    Returns the sampled IMF up to a total mass.
    """

    # # Select a fraction of stars to be binaries
    # bin_indxs = binar_probs[:isoch_cut.shape[-1]] <= bin_frac
    # # Secondary masses of the binary systems
    # binar_masses = isochrone[-1][bin_indxs]
    # # Excess of mass that will be introduced by the secondary masses
    # M_excess = binar_masses.sum()
    # # Mass of primary systems
    # single_masses = isochrone[m_ini_idx][bin_indxs].sum()
    # # Excess mass as a fraction of the single masses
    # M_excess_frac = M_excess / single_masses
    # # Remove that excess from the total mass to compensate
    # M_total = max(10, M_total - M_excess_frac * bin_frac * M_total)
    # #
    # This line is a faster approximation to the (exact) block above
    M_total = max(10, M_total - mean_bin_mr * bin_frac * M_total)

    return st_dist_mass[0][:np.searchsorted(st_dist_mass[1], M_total)]
