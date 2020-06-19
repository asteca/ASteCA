
import numpy as np
from scipy.spatial.distance import cdist
from .contamination_index import CIfunc
from ..aux_funcs import circFrac
from ..out import prep_plots
from .. import update_progress


def main(cld_i, clp, coords, rad_manual, nsteps_rad, NN_dd, **kwargs):
    """
    Estimate the radius through the optimization of the #cluster-members vs
    #field-stars values. Assign the uncertainty through a bootstrap process
    based on the field density's uncertainty.
    """

    # Parameters used internally only.
    rad_radii, rad_areas, st_dists_cent = rdpAreasDists(cld_i, **clp)

    # Obtain the optimal radius and arrays for plotting.
    clp['clust_rad'], clp['rad_rads'], clp['rad_N_membs'], clp['rad_N_field'],\
        clp['rad_CI'] = optimalRadius(
            rad_radii, rad_areas, st_dists_cent, **clp)

    coord = prep_plots.coord_syst(coords)[0]
    if rad_manual == 'n':
        print("Estimating radius")

        if nsteps_rad > 2:
            # Use bootstrap to estimate the uncertainty.
            clp['e_rad'] = radError(
                rad_radii, rad_areas, st_dists_cent, nsteps_rad, **clp)
            print("Radius found: {:g} {}".format(clp['clust_rad'], coord))
        else:
            clp['e_rad'] = np.array([np.nan, np.nan])
            print("Radius found (no bootstrap): {:g} {}".format(
                clp['clust_rad'], coord))

    elif rad_manual != 'n':
        # Update radius, assign zero error.
        clp['clust_rad'], clp['e_rad'] = rad_manual, np.array([np.nan, np.nan])
        print("Manual radius set: {:g} {}".format(clp['clust_rad'], coord))

    return clp


def rdpAreasDists(
    cld_i, kde_cent, N_MC, rand_01_MC, cos_t, sin_t, Nstars=20, Nrads=300,
        **kwargs):
    """
    The areas for each radius value in 'rad_radii' are obtained here once.
    We also calculate here the distance of each star to the defined center.

    HARDCODED
    Nstars: number of stars used to define the minimum and maximum radii values
    used to define the 'rad_radii' array.
    Nrads: number of values used to generate the 'rad_radii' array.
    """

    # Array of coordinates.
    xy = np.array([cld_i['x'], cld_i['y']]).T

    # Distances of all stars to the center of the cluster.
    st_dists_cent = cdist([kde_cent], xy)[0]
    sort_idx = np.argsort(st_dists_cent)

    # Use min,max radius values as those that leave 20 stars before and after,
    # respectively.
    ri_min, ri_max = st_dists_cent[
        sort_idx[Nstars]], st_dists_cent[sort_idx[-Nstars]]
    # This array gives the radius finding function a reasonable resolution.
    rad_radii = np.linspace(ri_min, ri_max, Nrads)

    # Frame limits
    x0, x1 = min(xy.T[0]), max(xy.T[0])
    y0, y1 = min(xy.T[1]), max(xy.T[1])

    dx0, dx1 = abs(kde_cent[0] - x0), abs(kde_cent[0] - x1)
    dy0, dy1 = abs(kde_cent[1] - y0), abs(kde_cent[1] - y1)
    dxy = min(dx0, dx1, dy0, dy1)

    # Areas associated to the radii defined in 'rad_radii'.
    rad_areas = np.pi * np.array(rad_radii)**2
    for i, rad in enumerate(rad_radii):
        fr_area = 1.
        if rad > dxy:
            fr_area = circFrac(
                (kde_cent), rad, x0, x1, y0, y1, N_MC, rand_01_MC, cos_t,
                sin_t)
        rad_areas[i] *= fr_area

    return rad_radii, rad_areas, st_dists_cent


def optimalRadius(
    rad_radii, rad_areas, st_dists_cent, field_dens, bttstrp_flag=False,
        **kwargs):
    """
    Estimate the optimal radius as the value that maximizes the (normalized)
    number of members stars versus the number of field stars. The rationale
    is that we want the maximum possible of member stars, ans the minimum
    possible of field stars within the region.
    """

    data, break_counter = [], 0
    for i, rad in enumerate(rad_radii):

        # Stars within radius.
        n_in_cl_reg = (st_dists_cent <= rad).sum()
        if n_in_cl_reg == 0:
            continue

        # Tried (21/01/20) the variable field density method, but it failed.
        # It is too fragile and ends up selecting very large radius values.
        # Field density
        # (dist > rad).sum(): #stars outside the cluster region
        # (area_frame - rdp_areas[i]): are outside the cluster region.
        # field_dens_i = (N_tot - n_in_cl_reg) / (area_tot - rad_areas[i])
        # fd_list.append(field_dens_i)
        # field_dens = np.median(fd_list)

        CI, n_memb, n_fl = CIfunc(n_in_cl_reg, field_dens, rad_areas[i])

        # This keeps the N_memb trending always upwards. Nicer graphs but
        # hides information.
        # n_memb = max(n_memb_old, n_memb)
        # n_memb_old = n_memb

        # Break check
        if (n_memb <= 0. or CI > 1.1) and i > .5 * len(rad_radii):
            break_counter += 1
            if break_counter > 3:
                break
        data.append([rad, n_memb, n_fl, CI])

    # rads, N_membs, N_field, CI
    data = np.clip(data, a_min=0., a_max=None).T

    # Normalizing separately is important. Otherwise it selects low radii
    # values.
    N_membs = data[1] / data[1].max()
    N_field = data[2] / data[2].max()
    idx = np.argmax(N_membs - N_field)
    clust_rad = data[0][idx]

    if bttstrp_flag:
        return clust_rad

    return clust_rad, data[0], N_membs, N_field, data[3]


def radError(
        rad_radii, rad_areas, st_dists_cent, N_btstrp, field_dens, **kwargs):
    """
    Bootstrap the distances to estimate the radius uncertainty.
    """

    steps = int(.1 * N_btstrp)
    all_rads = np.empty(N_btstrp)
    for i in range(N_btstrp):

        # Sample distances with replacement.
        st_dists_btstrp = np.random.choice(st_dists_cent, st_dists_cent.size)

        # Obtain the optimal radius and arrays for plotting.
        all_rads[i] = optimalRadius(
            rad_radii, rad_areas, st_dists_btstrp, field_dens, True)

        if i + 1 == N_btstrp:
            pass
        elif (i + 1) % steps:
            continue
        update_progress.updt(N_btstrp, i + 1)

    return np.percentile(all_rads, (16, 84))
