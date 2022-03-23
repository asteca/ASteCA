

def check(pd):
    """
    Check that the parameters are properly written.
    """

    manual_struct = {}
    for line in pd['manual_struct']:
        if line[0] in manual_struct.keys():
            raise ValueError((
                "duplicated entries found in 'Structure' block: '{}'".format(
                    line[0])))

        if line[1] == 'a' or line[2] == 'a':
            if not (line[1] == 'a' and line[2] == 'a'):
                raise ValueError(
                    "Center coordinates must both be either floats or 'a'")
            cx, cy = 'a', 'a'
        else:
            cx, cy = float(line[1]), float(line[2])

        float_flag = False
        try:
            fdens = float(line[3])
            float_flag = True
        except ValueError:
            fdens = line[3]
        if float_flag:
            if fdens < 0:
                raise ValueError("The field density value must be >= 0")
        else:
            if fdens != 'a':
                raise ValueError("Unrecognized field density mode '{}'".format(
                    line[3]))

        float_flag = False
        try:
            rad = float(line[4])
            float_flag = True
        except ValueError:
            rad = line[4]
        if float_flag:
            if rad <= 0:
                raise ValueError("The radius value must be > 0")
        else:
            if rad not in pd['rad_modes_accpt']:
                raise ValueError("Radius mode '{}' not recognized".format(
                    rad))

        float_flag = False
        try:
            fregs = int(line[5])
            float_flag = True
        except ValueError:
            fregs = line[5]
        if float_flag:
            if fregs < 0:
                raise ValueError("The number of field regions must be >= 0")
        else:
            if fregs != 'a':
                raise ValueError("Unrecognized field regions mode '{}'".format(
                    fregs))

        manual_struct[line[0]] = (cx, cy, fdens, rad, fregs)

    # Re-write entry in 'pd'
    pd['manual_struct'] = manual_struct

    if pd['kp_ndim'] not in (0, 2, 4):
        raise ValueError(
            "Unrecognized value for King profile 'ndim' parameter")
    elif pd['kp_ndim'] in (2, 4):
        # # DEPRECATED 05/2021
        # if 'emcee' not in pd['inst_packgs_lst']:
        #     raise ValueError("King profile is selected to run, but 'emcee'"
        #                      " is not installed")
        if pd['kp_nchains'] < 10:
            raise ValueError(
                "set a minimum of 10 chains for KP Bayesian analysis")
        if pd['kp_nburn'] <= 0. or pd['kp_nburn'] >= 1.:
            raise ValueError("KP 'nburn' should be in the range (0., 1.)")

    return pd
