
import numpy as np
from .. import update_progress


def main(data_tracks, st_dist_mass, interp_tracks, binar_flag):
    """
    Interpolate sampled masses into all the filters, colors, filters of colors,
    and extra parameters (masses, etc).
    """

    # 'data_tracks' is defined in read_met_files() as:
    # pd['data_tracks'] = [mags_theor, cols_theor, extra_pars, mags_cols_theor]
    extra_pars = data_tracks[-2]

    # Sampled masses
    t = st_dist_mass[0]

    # Do not add the 'mags_cols_theor' array (at the end of 'data_tracks') to
    # 'interp_tracks'
    idx = 0
    for nn, data in enumerate(data_tracks[:-1]):
        # For each metallicity value list.
        for i, met in enumerate(data):
            # For each age value list.
            for j, age in enumerate(met):
                # Assume 'M_ini' is in the '0th' position
                xp = extra_pars[i][j][0]
                # For each filter/color/extra parameter in list.
                for k, fce in enumerate(age):
                    if i == 0 and j == 0:
                        idx += k
                    interp_tracks[i][j][nn + idx] = np.interp(t, xp, fce)

        update_progress.updt(len(data_tracks[:-1]), nn + 1)

    # The magnitudes for each defined color ('mags_cols_intp') are used here
    # and discarded after the colors (and magnitudes) with binarity assignment
    # are obtained.
    mags_cols_intp = []
    if binar_flag:
        for i, met in enumerate(data_tracks[-1]):
            m = []
            for j, age in enumerate(met):
                a = []
                for fce in age:
                    xp = extra_pars[i][j][0]
                    a.append(np.interp(t, xp, fce))
                m.append(a)
            mags_cols_intp.append(m)

    # # In place for #415
    # filename = 'temp.memmap'
    # sp = np.array(theor_tracks).shape
    # theor_tracks_mp = np.memmap(
    #     filename, dtype='float64', mode='w+', shape=sp)
    # theor_tracks_mp[:] = theor_tracks[:]

    return interp_tracks, mags_cols_intp
