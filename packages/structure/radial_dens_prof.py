
import numpy as np
from ..aux_funcs import circFrac


def main(clp, RDP_rings, rings_rm=.1, Nmin=10, **kwargs):
    """
    Obtain the RDP using the concentric rings method. For plotting only.

    HARDCODED

    rings_rm: remove the more conflicting last X% of radii values.
    Nmin: minimum number of stars that a ring should contain. Else, expand it.
    """

    # Frame limits
    x0, x1 = min(clp['xy_filtered'].T[0]), max(clp['xy_filtered'].T[0])
    y0, y1 = min(clp['xy_filtered'].T[1]), max(clp['xy_filtered'].T[1])

    # Handle the case where int()==0
    max_i = max(1, int(rings_rm * RDP_rings))
    # The +1 adds a ring accounting for the initial 0. in the array
    radii = np.linspace(
        0., clp['xy_cent_dist'].max(), RDP_rings + 1 + max_i)[:-max_i]

    # Areas and #stars for all rad values.
    rdp_radii, rdp_points, rdp_stddev = [], [], []
    l_prev, N_in_prev = np.inf, 0.
    for l, h in zip(*[radii[:-1], radii[1:]]):
        # Stars within this ring.
        N_in = (
            (clp['xy_cent_dist'] >= l) & (clp['xy_cent_dist'] < h)).sum() +\
            N_in_prev

        l_now = min(l, l_prev)

        # Require that at least 'Nmin' stars are within the ring.
        if N_in > Nmin:
            # Area of ring.
            fr_area_l = circFrac(
                (clp['kde_cent']), l_now, x0, x1, y0, y1, clp['N_MC'],
                clp['rand_01_MC'], clp['cos_t'], clp['sin_t'])
            fr_area_h = circFrac(
                (clp['kde_cent']), h, x0, x1, y0, y1, clp['N_MC'],
                clp['rand_01_MC'], clp['cos_t'], clp['sin_t'])
            ring_area = (np.pi * h**2 * fr_area_h) -\
                (np.pi * l_now**2 * fr_area_l)

            # Store RDP parameters.
            rad_med = h if l_now == 0. else .5 * (l_now + h)
            rdp_radii.append(rad_med)
            rdp_points.append(N_in / ring_area)
            rdp_stddev.append(np.sqrt(N_in) / ring_area)

            # Reset
            l_prev, N_in_prev = np.inf, 0.

        else:
            l_prev = l_now
            N_in_prev += N_in

    if not rdp_radii:
        raise ValueError("ERROR: RDP is empty. Check the center coordinates")

    if RDP_rings != len(rdp_radii):
        print("RDP: N={} rings with <10 stars inside were merged".format(
            RDP_rings - len(rdp_radii)))

    clp['rdp_radii'], clp['rdp_points'], clp['rdp_stddev'] = rdp_radii,\
        rdp_points, rdp_stddev
    return clp


# DELETE 23/01/20 Tried this approach briefly but it did not work.
# def kNNRDP(fr_dist, fr_dens):
#     """
#     Use the kNN's per-star densities. Average these values for several circular
#     rings to obtain the RDP.
#     """

#     # HARDCODED: remove the more conflicting last ~10% of radii values.
#     radii = np.linspace(0., fr_dist.max(), 55)[:-5]

#     rdp_radii, rdp_NN, rdp_stddev = [], [], []

#     # This method uses the median of the 'fr_dens' values for
#     # stars within a ring. Not sure why (09/01/20) it is not working
#     # properly.
#     for l, h in zip(*[radii[:-1], radii[1:]]):
#         # Stars within this ring.
#         msk_in = (fr_dist >= l) & (fr_dist < h)
#         if sum(msk_in) > 0:
#             rdp_radii.append(.5 * (l + h))
#             rdp_NN.append(np.median(fr_dens[msk_in]))
#             rdp_stddev.append(np.std(fr_dens[msk_in]))

#     if not rdp_radii:
#         raise ValueError("ERROR: RDP is empty. Check the center coordinates")

#     return rdp_radii, rdp_NN, rdp_stddev
