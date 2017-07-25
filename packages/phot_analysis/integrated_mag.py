
import numpy as np
from scipy.interpolate import spline


def calc_integ_mag(st_reg):
    '''
    Calculate integrated magnitude up to a certain maximum magnitude value.
    '''

    # Define magnitude range.
    mag_range = np.linspace(min(st_reg), max(st_reg), 100)

    # Define final output lists.
    reg_mag = [[], []]

    # Sort magnitude list.
    sort_lis = sorted(st_reg)

    int_mag_val = 0.
    for mag_limit in mag_range:

        # Calculate integrated magnitude up to this mag value.
        trim_lis = [i for i in sort_lis if i <= mag_limit]

        if len(trim_lis) == 0:
            int_mag_val = mag_limit
        else:
            energ_sum = 0.
            for mag_i in trim_lis:
                energ_sum = energ_sum + 10 ** (mag_i / -2.5)

            int_mag_val = -2.5 * np.log10(energ_sum)

        reg_mag[0].append(mag_limit)
        reg_mag[1].append(int_mag_val)

    return reg_mag


def field_reg_integ_mag_curve(fl_reg_m):
    '''
    Return smooth averaged curve for field star regions integrated magnitudes.
    '''
    # Average all field regions for the integrated magnitude.
    fl_reg_avr = np.average(fl_reg_m, axis=0)

    # Smooth curve for the averaged integrated magnitude.
    fl_reg_mag = []
    # Magnitude values.
    fl_reg_mag.append(np.linspace(min(fl_reg_avr[0]), max(fl_reg_avr[0]), 100))
    # Integ magnitude values.
    fl_reg_mag.append(spline(fl_reg_avr[0], fl_reg_avr[1], fl_reg_mag[0]))

    return fl_reg_mag


def main(clp, im_flag, **kwargs):
    '''
    Obtain integrated magnitude using all stars inside the cluster's radius for
    several limits in magnitude.
    '''
    cl_region, field_regions, flag_no_fl_regs = [
        clp[_] for _ in ['cl_region', 'field_regions', 'flag_no_fl_regs']]
    if im_flag:
        # Only use stars inside cluster's radius.
        # For each magnitude defined.
        cl_reg_imag = []
        for mag in zip(*zip(*cl_region)[3]):
            cl_reg_imag.append(calc_integ_mag(mag))

        if flag_no_fl_regs is False:

            # Run for every field region defined.
            fl_regs_int_m = []
            for f_reg in field_regions:
                fl_reg_m = []
                for mag in zip(*zip(*f_reg)[3]):
                    fl_reg_m.append(calc_integ_mag(mag))
                fl_regs_int_m.append(fl_reg_m)

            # For each defined magnitude.
            fl_reg_imag = []
            for f_mags in zip(*fl_regs_int_m):
                fl_reg_imag.append(field_reg_integ_mag_curve(f_mags))

            # Obtain integrated magnitude of clean cluster region, ie:
            # subtracting the field contribution.
            integ_mag = []
            for fl_m, cl_m in zip(fl_reg_imag, cl_reg_imag):
                if min(fl_m[1]) >= min(cl_m[1]):
                    integ_mag.append(-2.5 * np.log10(
                        1 - 10**((min(fl_m[1]) - min(cl_m[1])) / -2.5)) +
                        min(cl_m[1]))
                else:
                    print("  WARNING: integrated magnitude of field regions\n"
                          "  is larger than that of cluster+field region.\n"
                          "  Storing cluster+field region integrated mag.")
                    # If the field is brighter than the cluster.
                    integ_mag.append(min(cl_m[1]))
        else:
            print("  WARNING: no field regions defined. Integrated magnitude\n"
                  "  is not cleaned from field star contamination.")
            # Pass dummy lists.
            fl_reg_imag, integ_mag = [], [np.nan]
            for cl_m in cl_reg_imag:
                integ_mag.append(min(cl_m[1]))

        print('Integrated magnitude distribution obtained ({:.2f}).'.format(
            integ_mag[0]))
    else:
        print('Skipping integrated magnitudes function.')
        cl_reg_imag, fl_reg_imag, integ_mag = [], [], [np.nan]

    clp['cl_reg_imag'], clp['fl_reg_imag'], clp['integ_mag'] =\
        cl_reg_imag, fl_reg_imag, integ_mag
    return clp
