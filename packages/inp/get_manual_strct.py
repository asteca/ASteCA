

def main(pd, clust_name, **kwargs):
    """
    Read the range values for trimming the frame, if given.
    Read center, field density, radius and number of field regions for this
    cluster.
    """

    pd['flag_tf'], pd['tf_range'] = False, [0., 0., 0., 0.]
    for line in pd['trim_frame_range']:
        if line[0] == clust_name:
            pd['flag_tf'], pd['tf_range'] = True, line[1]

    # Default values
    pd['cent_method'], pd['fdens_method'], pd['rad_method'],\
        pd['fregs_method'] = ('a', 'a'), 'a', 'a', 'a'

    try:
        cx, cy, fdens, rad, fregs = pd['manual_struct'][clust_name]
        pd['cent_method'] = (cx, cy)
        pd['fdens_method'] = fdens
        pd['rad_method'] = rad
        pd['fregs_method'] = fregs
    except KeyError:
        pass

    return pd
