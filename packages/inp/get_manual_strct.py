

def main(pd, clust_name, **kwargs):
    """
    Read the range values for trimming the frame, if given.
    Read manual center, radius and number of field regions for this cluster.
    """

    pd['flag_tf'], pd['tf_range'] = False, [0., 0., 0., 0.]
    for line in pd['trim_frame_range']:
        if line[0] == clust_name:
            pd['flag_tf'], pd['tf_range'] = True, line[1]

    # Fill with dummy values since these variables are required later on.
    pd['cent_manual'], pd['rad_manual'], pd['f_regs_manual'] = 'n', 'n', 'n'

    for line in pd['manual_struct']:
        # Cluster name found in file.
        if line[0] == clust_name:
            try:
                pd['cent_manual'] = [float(line[1]), float(line[2])]
            except ValueError:
                pd['cent_manual'] = 'n'
            try:
                pd['rad_manual'] = float(line[3])
            except ValueError:
                pd['rad_manual'] = 'n'
            try:
                pd['f_regs_manual'] = int(line[4])
            except ValueError:
                pd['f_regs_manual'] = 'n'
            break

    return pd
