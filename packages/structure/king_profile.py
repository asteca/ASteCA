
import numpy as np
from scipy.optimize import curve_fit
import warnings
from scipy.optimize import OptimizeWarning
import king_prof_funcs as kpf
from ..inp import input_params as g
from ..out import prep_plots


def fit_3P_King_prof(fd, radii_k, ring_dens_k, guess3):
    '''
    Fit central density, core radius and tidal radius, using a fixed
    value *only* for the field density (fd).
    '''

    with warnings.catch_warnings():
        warnings.simplefilter("error", OptimizeWarning)

        try:
            popt, pcov = curve_fit(lambda x, cd, rc, rt: kpf.three_params(
                                   x, rt, cd, rc, fd), radii_k, ring_dens_k,
                                   guess3)
        except:
            raise ValueError("  WARNING: (1) no convergence, core & tidal"
                             " radii not found.")

    # Unpack obtained core and tidal radius and calculate their errors.
    cd, rc, rt = popt
    # Obtain error in core and tidal radius.
    if np.isfinite(pcov).all():
        e_rc = np.sqrt(pcov[1][1]) if pcov[1][1] > 0 else -1.
        e_rt = np.sqrt(pcov[2][2]) if pcov[2][2] > 0 else -1.
    else:
        e_rc, e_rt = -1., -1.

    # If fit converged to tidal radius that extends beyond 500 times
    # the core radius, or either radius is equal or less than zero;
    # discard the fit.
    if rt > rc * 500.:
        raise ValueError("  WARNING: (1) tidal radius is too large.")
    elif rt <= 0. or rc <= 0.:
        raise ValueError("  WARNING: (1) core and/or tidal radius is <=0.")
    elif e_rc > (5. * rc) or e_rt > (5. * rt):
        raise ValueError("  WARNING: (1) core and/or tidal radius error is"
                         " too large.")

    return cd, rc, e_rc, rt, e_rt


def fit_2P_King_prof(fd, radii_k, ring_dens_k, guess2):
    '''
    Fit a 2P King profile. The maximum central density and core radius are
    left as free parameters while the field density is fixed (fd).
    '''

    with warnings.catch_warnings():
        warnings.simplefilter("error", OptimizeWarning)
        try:
            popt, pcov = curve_fit(lambda x, cd, rc: kpf.two_params(
                x, cd, rc, fd), radii_k, ring_dens_k, guess2)
        except:
            raise ValueError("  WARNING: (2) no convergence, core radius not"
                             " found.")

    # Unpack max density and core radius.
    cd, rc = popt
    # Obtain error in core radius.
    if np.isfinite(pcov).all():
        e_rc = np.sqrt(pcov[1][1]) if pcov[1][1] > 0 else -1.
    else:
        e_rc = -1.

    # If fit converged, check values.
    if rc <= 0.:
        raise ValueError("  WARNING: (2) core radius is <=0.")
    elif e_rc > (5. * rc):
        raise ValueError("  WARNING: (2) core radius error is too large.")

    return cd, rc, e_rc


def fit_3P_2P_King_prof(fd, rc, radii_k, ring_dens_k, guess3):
    '''
    Fit central density and tidal radius, using a fixed value for the field
    density (fd) and the core radius (rc).
    '''

    with warnings.catch_warnings():
        warnings.simplefilter("error", OptimizeWarning)
        try:
            popt, pcov = curve_fit(lambda x, cd, rt: kpf.three_params(
                                   x, rt, cd, rc, fd), radii_k, ring_dens_k,
                                   guess3)
        except:
            raise ValueError("  WARNING: (3) no convergence, tidal radius not"
                             " found.")

    # Unpack tidal radius and its error.
    rt = popt[1]
    # Obtain error in tidal radius.
    if np.isfinite(pcov).all():
        e_rt = np.sqrt(pcov[1][1]) if pcov[1][1] > 0 else -1.
    else:
        e_rt = -1.

    # If the concentration parameter is larger than ~2.7 (too
    # large even for a globular cluster) or either radius is
    # equal or less than zero; discard the fit.
    if (rt / rc) > 500.:
        raise ValueError("  WARNING: (3) tidal radius is too large.")
    elif rt <= 0.:
        raise ValueError("  WARNING: (3) tidal radius is <=0.")
    elif e_rt > (5. * rt):
        raise ValueError("  WARNING: (3) tidal radius error is too large.")

    return rt, e_rt


def num_memb_conc_param(flag_3pk_conver, cd, rt, rc):
    '''
    If 3-P King profile converged, ie: the tidal radius was found,
    calculate approximate number of cluster members with Eq (3) from
    Froebrich et al. (2007); 374, 399-408 and the concentration
    parameter.
    '''

    n_c_k, kcp = -1., -1.
    if flag_3pk_conver:
        # Obtain approximate number of members.
        x = 1 + (rt / rc) ** 2
        n_c_k = int(round((np.pi * cd * rc ** 2) * (np.log(x) -
                    4 + (4 * np.sqrt(x) + (x - 1)) / x)))
        # Obtain concentration parameter.
        kcp = np.log10(rt / rc)

    return n_c_k, kcp


def main(clust_rad, field_dens, radii, rdp_points):
    '''
    Function to fit a King profile to a given radial density.
    The field density value is fixed and the core radius, tidal radius and
    maximum central density are fitted.
    '''

    # Flags that indicate either no convergence or that the fits were not
    # attempted.
    flag_2pk_conver, flag_3pk_conver = False, False

    # Initial dummy values, used if function is skipped.
    rc, e_rc, rt, e_rt, n_c_k, kcp, cd = -1., -1., -1., -1., -1., -1., -1.

    # Check flag to run or skip.
    if g.kp_flag:

        # Field density value is fixed.
        fd = field_dens
        # Initial guesses for fit: max_dens, rt, rc
        max_dens, rt_guess, rc_guess = max(rdp_points), clust_rad, \
            clust_rad / 2.

        # Skip first radius value if it is smaller than the second value. This
        # makes it easier for the KP to converge.
        if rdp_points[0] > rdp_points[1]:
            radii_k, ring_dens_k = radii, rdp_points
        else:
            radii_k, ring_dens_k = radii[1:], rdp_points[1:]

        # USE AT SOME POINT?
        # # Find maximum density value and assume this is the central density.
        # # Do not use previous values.
        # max_dens_ind = np.argmax(rdp_points)
        # radii_k, ring_dens_k = radii[max_dens_ind:],rdp_points[max_dens_ind:]

        # Attempt to fit a 3P-KP, fixing *only* the field density.
        try:
            guess3 = (max_dens, rc_guess, rt_guess)
            cd, rc, e_rc, rt, e_rt = fit_3P_King_prof(fd, radii_k,
                                                      ring_dens_k,
                                                      guess3)
            # Fit converged.
            flag_3pk_conver = True

        except Exception as e:
            print(e)
            flag_3pk_conver = False

            # Attempt to fit a 2P-KP to obtain the core radius.
            try:
                guess2 = (max_dens, rc_guess)
                cd, rc, e_rc = fit_2P_King_prof(fd, radii_k, ring_dens_k,
                                                guess2)
                # Fit converged.
                flag_2pk_conver = True

                # Attempt to fit a 3P-KP, fixing the obtained r_core value, and
                # the field density.
                try:
                    guess3 = (max_dens, rt_guess)
                    rt, e_rt = fit_3P_2P_King_prof(fd, rc, radii_k,
                                                   ring_dens_k, guess3)
                    # Fit converged.
                    flag_3pk_conver = True

                except Exception as e:
                    print(e)

            except Exception as e:
                print(e)
                # 2P-KP did not converge.
                flag_2pk_conver = False

        # Obtain number of members and concentration parameter.
        n_c_k, kcp = num_memb_conc_param(flag_3pk_conver, cd, rt, rc)

        # Print results.
        coord = prep_plots.coord_syst()[0]
        if flag_3pk_conver:
            # Set precision of printed values.
            text2 = '{:.1f}, {:.1f}' if coord == 'px' else '{:g}, {:g}'
            text = 'Core & tidal radii obtained: ' + text2 + ' {}.'
            print text.format(rc, rt, coord)

        elif flag_2pk_conver:
            # Set precision of printed values.
            text2 = '{:.1f}' if coord == 'px' else '{:g}'
            text = 'Only core radius obtained: ' + text2 + ' {}.'
            print text.format(rc, coord)

        else:
            print("Core & tidal radii not found.")

    return rc, e_rc, rt, e_rt, n_c_k, kcp, cd, flag_2pk_conver, flag_3pk_conver
