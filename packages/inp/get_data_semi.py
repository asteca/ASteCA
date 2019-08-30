

def main(pd, clust_name, **kwargs):
    '''
    Get center, radius and flags for semi automatic mode.
    '''
    if pd['run_mode'] == 'semi':
        for line in pd['semi_input']:
            # Cluster name found in file.
            if line[0] == clust_name:
                pd['cl_cent_semi'] = [float(line[1]), float(line[2])]
                pd['cl_rad_semi'] = float(line[3])
                pd['cl_f_regs_semi'] = int(line[4])
                pd['cent_flag_semi'] = int(line[5])
                pd['rad_flag_semi'] = int(line[6])
                pd['freg_flag_semi'] = int(line[7])
    else:
        # Fill with dummy values since these variables are required later on.
        pd['cl_cent_semi'], pd['cl_rad_semi'], pd['cl_f_regs_semi'],\
            pd['cent_flag_semi'], pd['rad_flag_semi'],\
            pd['freg_flag_semi'] = [], -1., 0, 0, 0, 0

    return pd
