
import copy


def main(pd, clust_name, **kwargs):
    """
    Range values for trimming the frame, if given.

    Read center, field density, radius and number of field regions for this
    cluster. If this cluster has no values set for these parameters, use the
    general values from the 'CLUSTER' entry.
    """

    pd['flag_tf'], pd['tf_range'] = False, [0., 0., 0., 0.]
    for line in pd['trim_frame_range']:
        if line[0] == clust_name:
            pd['flag_tf'], pd['tf_range'] = True, line[1]

    clFlag = False
    try:
        cx, cy, fdens, rad, fregs = copy.deepcopy(
            pd['manual_struct'][clust_name])
        clFlag = True
    except KeyError:
        pass
    # Default values
    if clFlag is False:
        if 'CLUSTER' in pd['manual_struct'].keys():
            cx, cy, fdens, rad, fregs = copy.deepcopy(
                pd['manual_struct']['CLUSTER'])
        else:
            raise ValueError("Line 'S0' in 'asteca.ini' missing 'CLUSTER'")

    pd['cent_method'] = (cx, cy)
    pd['fdens_method'] = fdens
    pd['rad_method'] = rad
    pd['fregs_method'] = fregs

    return pd
