

def main(clp, memb_file_out, **kwargs):
    '''
    Create output data file with stars inside the cluster radius along with
    their membership probabilities.
    '''

    cl_reg_fit, cl_reg_no_fit = clp['cl_reg_fit'], clp['cl_reg_no_fit']

    with open(memb_file_out, 'w') as out_data_file:
        out_data_file.write("#ID x y mag e_mag col e_col memb_prob sel\n")

    txt = '{:<8}  {:>12.6f}  {:>12.6f}  {:>8.3f}  {:>8.3f}  {:>8.3f}  ' +\
        '{:>8.3f}  {:>8.2f}  {:>6}\n'

    # Save cluster region with MPs obtained by the decontamination algorithm.
    with open(memb_file_out, "a") as f_out:
        for line in cl_reg_fit:
            # Identify stars selected by the removal function.
            l = line + ['1']
            f_out.write(txt.format(*l))
        # Stars rejected by the removal function.
        for line in cl_reg_no_fit:
            l = line + ['0']
            f_out.write(txt.format(*l))

    print('Membership probabilities saved to file.')
