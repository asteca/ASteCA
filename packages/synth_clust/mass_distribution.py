
import numpy as np


def main(st_dist_mass, mean_bin_mr, bin_frac, M_total):
    """
    Returns the sampled IMF up to a total mass.

    The binary fraction correction is derived as follows:

    1. The initial mass of all single systems is: M_T = sum(m_s)_1^N_T
    2. A fraction 'bin_frac' is converted to binary systems:
       N_B = bin_frac * N_T
    3. The new (larger) mass is thus: M_F = M_T + M_B = sum(m_s)_1^N_T +
       sum(m_b)_1^N_B
    4. The excess mass introduced is : M_B = sum(m_b)_1^N_B ~ N_B * m_b' =
       bin_frac * N_T * m_b'; where m_b' is the mean secondary mass of the
       binary systems.
    5. The mean secondary mass can be obtained as: m_b' = m_s' * r; where
       m_s' is the mean single systems' mass, and r is the mean mass ratio
       of the 'q' distribution.
    6. Hence: M_T = N_T * m_s' <--> m_s' = M_T/N_T
    7. Finally: M_B = M_T * bin_frac * r

    """
    # Correct the total mass
    M_total = max(10, M_total - mean_bin_mr * bin_frac * M_total)

    return st_dist_mass[0][:np.searchsorted(st_dist_mass[1], M_total)]
