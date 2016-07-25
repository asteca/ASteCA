
from ..inp import input_params as g


def main(synth_clst, synth_file_out):
    '''
    Create output data file with stars in the best fit synthetic cluster found
    by the 'Best Fit' function.
    '''

    bf_flag = g.bf_params[0]

    # Check if function should run.
    if bf_flag:

        # If cluster is not empty.
        if synth_clst.any():

            # Save best fit synthetic cluster found to file.
            with open(synth_file_out, "w") as f_out:
                f_out.write("#mag    e_mag    col1   e_col1    m_ini\n")
                for line in zip(*synth_clst):
                    f_out.write(
                        "{:<8.4f} {:>8.4f} {:>8.4f} {:>8.4f}"
                        " {:>8.4f}\n".format(line[2], line[3], line[0],
                                             line[1], line[4]))

            print('Best fit synthetic cluster saved to file.')

        else:
            print("  WARNING: empty synthetic cluster could not be saved\n"
                  "  to file.")
