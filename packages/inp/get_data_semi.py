
import input_params as g


def main(clust_name):
    '''
    Get center, radius and flags for semi automatic mode.
    '''

    # Define as global so as to be able to change it if the cluster is not
    # found in semi_input.dat file.
    global mode

    semi_return = []
    # Mode is semi.
    if g.mode == 'semi':

        semi_file = 'semi_input.dat'
        # Flag to indicate if cluster was found in file.
        flag_clust_found = False
        with open(semi_file, "r") as f_cl_dt:
            for line in f_cl_dt:
                li = line.strip()
                # Skip comments.
                if not li.startswith("#"):
                    reader = li.split()

                    # Prevent empty lines with spaces detected as a cluster
                    # line from crashing the code.
                    if reader:
                        # If cluster is found in file.
                        if reader[0] == clust_name:
                            cl_cent_semi = [float(reader[1]), float(reader[2])]
                            cl_rad_semi = float(reader[3])
                            cl_f_regs_semi = int(reader[4])
                            cent_flag_semi, rad_flag_semi, freg_flag_semi, \
                                err_flag_semi = int(reader[5]), \
                                int(reader[6]), int(reader[7]), int(reader[8])
                            # Set flag to True if the cluster was found.
                            flag_clust_found = True

        # If cluster was found.
        if flag_clust_found:
            semi_return = [cl_cent_semi, cl_rad_semi, cl_f_regs_semi,
                           cent_flag_semi, rad_flag_semi, freg_flag_semi,
                           err_flag_semi]
        else:
            # If the cluster was not found in the file, default to 'manual'.
            print ("  WARNING: cluster's name not found in 'semi_input.dat'\n"
                   "  file. Defaulting to 'auto' mode.")
            # Re-define global variable.
            g.mode = 'auto'

    return semi_return
