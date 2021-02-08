
import numpy as np
from . import extin_coefs
from . import imf
from . import binarity
from . import add_errors
from .. import update_progress


def main(pd, clp):
    """
    Return list structured as:

    theor_tracks = [m1, m2, .., mN]
    mX = [age1, age2, ..., ageM]
    ageX = [f1,.., c1, c2,.., f1b,.., c1b, c2b,.., bp, mb, m_ini,.., m_bol]
    where:
    fX:  individual filters (mags). Currently a single main mag is allowed.
    cX:  colors
    Mini,...: extra parameters.
    fXb: filters with binary data  ------
    cXb: colors with the binary data    | Only stored if binarity is set as a
    mb:  binary masses                  | free param, or fixed to a value > 0.
    pb:  binary probabilities------------

    theor_tracks.shape = (Nz, Na, Nd, Ni)
    Nz: number of metallicities
    Na: number of log(age)s
    Nd: number of data columns
    Ni: number of interpolated values
    """
    # Set the random seed in numpy
    from . import set_rand_seed
    set_rand_seed.main(pd['synth_rand_seed'])

    # Obtain extinction coefficients.
    ext_coefs = extin_coefs.main(pd['cmd_systs'], pd['filters'], pd['colors'])

    # Obtain mass distribution using the selected IMF.
    st_dist_mass = imf.main(pd['IMF_name'], pd['fundam_params'][4])

    # Set the binary flag
    binar_fracs = pd['fundam_params'][5]
    pd['binar_flag'] = False
    if len(binar_fracs) > 1 or binar_fracs[0] > 0.:
        pd['binar_flag'] = True

    # Store the number of defined filters and colors.
    N_fc = [len(pd['filters']), len(pd['colors'])]
    # Index of 'M_ini' (theoretical initial mass), stored in the interpolated
    # isochrones: right after the magnitude and color(s)
    pd['m_ini_idx'] = N_fc[0] + N_fc[1]

    # Interpolate the sampled masses in the isochrones
    print("Populating isochrones")
    interp_tracks, mags_cols_intp = interpIsochs(
        pd['data_tracks'], st_dist_mass, pd['binar_flag'], N_fc)

    # Combine all data into a single array of shape:
    # (N_z, N_age, N_data, N_IMF_interp), where 'N_data' depends on the number
    # of colors defines, and whether the binarity is on/off
    if pd['binar_flag']:
        all_met_vals, all_age_vals = pd['fundam_params'][:2]
        interp_tracks = binarity.binarGen(
            pd['data_tracks'], pd['m_ini_idx'], interp_tracks, mags_cols_intp,
            all_met_vals, all_age_vals, pd['bin_mr'], N_fc)
    else:
        interp_tracks = interp_tracks[:, :, :pd['m_ini_idx'] + 1, :]

    pd['theor_tracks'] = interp_tracks

    print("(Size of array: {:.0f} Mbs)".format(
        pd['theor_tracks'].nbytes / 1024.**2))

    # Error parameters
    err_rand = add_errors.randIdxs(
        pd['lkl_method'], pd['theor_tracks'].shape[-1])

    # For the move_isoch() function. In place for #174
    rand_unif = np.random.uniform(0., 1., pd['theor_tracks'].shape[-1])

    # Pack
    mv_err_pars = clp['err_lst'], clp['em_float'], err_rand, rand_unif

    return pd, ext_coefs, st_dist_mass, N_fc, mv_err_pars


def interpIsochs(data_tracks, st_dist_mass, binar_flag, N_fc):
    """
    Interpolate sampled masses into all the filters, colors, filters of colors,
    and extra parameters (masses, etc).
    """

    # Sampled masses
    t = st_dist_mass[0]

    # Define empty array to hold the final data
    # Magnitude and colors defined. The '1' is there for the initial mass
    Nd = N_fc[0] + N_fc[1] + 1
    if binar_flag:
        # The '2' is there for the binarity probability and mass
        Nd += N_fc[0] + N_fc[1] + 2
    Nz, Na, Np = len(data_tracks[0][0]), len(data_tracks[0][0][0]), len(t)
    interp_tracks = np.empty([Nz, Na, Nd, Np])

    # 'data_tracks' is defined in read_met_files() as:
    # pd['data_tracks'] = [mags_theor, cols_theor, extra_pars, mags_cols_theor]
    # mags_theor, cols_theor, extra_pars, mags_cols_theor = data_tracks
    extra_pars = data_tracks[-2]

    # Do not add the 'mags_cols_theor' array (at the end of 'data_tracks') to
    # 'interp_tracks'
    kk = -1
    for nn, data in enumerate(data_tracks[:-1]):
        kk += 1
        # For each filter/color/extra parameter in list.
        for idx, fce in enumerate(data):
            kk += idx
            # For each metallicity value list.
            for i, met in enumerate(fce):
                # For each age value list.
                for j, age in enumerate(met):
                    # Assume 'M_ini' is in the '0th' position
                    xp = extra_pars[0][i][j]
                    yp = fce[i][j]
                    interp_tracks[i][j][kk] = np.interp(t, xp, yp)

        update_progress.updt(len(data_tracks[:-1]), nn + 1)

    # The magnitudes for each color ('mags_cols_intp') are defined here
    # and discarded after the colors (and magnitudes) with binarity assignment
    # are obtained.
    Nz, Na, _, Np = interp_tracks.shape
    mags_cols_intp = np.empty([Nz, Na, len(data_tracks[-1]), Np])
    if binar_flag:
        for k, fce in enumerate(data_tracks[-1]):
            for i, met in enumerate(fce):
                for j, age in enumerate(met):
                    xp = extra_pars[0][i][j]
                    yp = fce[i][j]
                    mags_cols_intp[i][j][k] = np.interp(t, xp, yp)

    # # In place for #415
    # filename = 'temp.memmap'
    # sp = np.array(theor_tracks).shape
    # theor_tracks_mp = np.memmap(
    #     filename, dtype='float64', mode='w+', shape=sp)
    # theor_tracks_mp[:] = theor_tracks[:]

    return interp_tracks, mags_cols_intp
