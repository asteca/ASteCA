

def main(pd, clust_name, **kwargs):
    '''
    Get center, radius and flags for semi automatic mode.
    '''
    semi_file = 'semi_input.dat'  # HARDCODED
    if pd['run_mode'] == 'semi':
        with open(semi_file, "r") as f_cl_dt:
            for line in f_cl_dt:
                li = line.strip()
                # Skip comments.
                if not li.startswith("#"):
                    reader = li.split()
                    # Prevent empty lines with spaces detected as a cluster
                    # line from crashing the code.
                    if reader:
                        # Cluster name found in file.
                        if reader[0] == clust_name:
                            pd['cl_cent_semi'] = [float(reader[1]),
                                                  float(reader[2])]
                            pd['cl_rad_semi'] = float(reader[3])
                            pd['cl_f_regs_semi'] = int(reader[4])
                            pd['cent_flag_semi'] = int(reader[5])
                            pd['rad_flag_semi'] = int(reader[6])
                            pd['freg_flag_semi'] = int(reader[7])
                            pd['err_flag_semi'] = int(reader[8])
    else:
        # Fill with dummy values since these variables are required later on.
        pd['cl_cent_semi'], pd['cl_rad_semi'], pd['cl_f_regs_semi'],\
            pd['cent_flag_semi'], pd['rad_flag_semi'],\
            pd['freg_flag_semi'], pd['err_flag_semi'] = [], -1., 0, 0,\
            0, 0, 0

    return pd
