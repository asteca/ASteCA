

def main(clp, memb_file_out, **kwargs):
    '''
    Create output data file with stars inside the cluster radius along with
    their membership probabilities.
    '''

    cl_reg_fit, cl_reg_no_fit = clp['cl_reg_fit'], clp['cl_reg_no_fit']

    with open(memb_file_out, 'w') as out_data_file:
        out_data_file.write(
            "#ID                  x             y       "
            "mag     e_mag       col     e_col memb_prob     sel\n")

    frmt = '{:<8}  {:>12.6f}  {:>12.6f}  {:>8.3f}  {:>8.3f}  {:>8.3f}  ' +\
        '{:>8.3f}  {:>8.2f}  {:>6}\n'

    idx = ['1', '0']
    # Save cluster region with MPs obtained by the decontamination algorithm.
    with open(memb_file_out, "a") as f_out:
        for i, reg in enumerate([cl_reg_fit, cl_reg_no_fit]):
            for line in reg:
                # Identify stars selected by the removal function.
                txt = line[:3] + [line[3][0]] + [line[4][0]] + [line[5][0]] +\
                    [line[6][0]] + [line[7]] + [idx[i]]
                f_out.write(frmt.format(*txt))

    print('Cluster region saved to file.')
