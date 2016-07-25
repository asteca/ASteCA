

def main(memb_file_out, red_return):
    '''
    Create output data file with stars inside the cluster radius along with
    their membership probabilities.
    '''

    red_memb_fit, red_memb_no_fit = red_return[:2]

    with open(memb_file_out, 'w') as out_data_file:
        out_data_file.write("#ID x y mag e_mag col1 e_col1 memb_prob sel\n")

    # Save cluster region with MPs obtained by the decont algor.
    with open(memb_file_out, "a") as f_out:
        for line in red_memb_fit:
            # Identify stars selected by the red_mem function.
            l = line + ['1']
            f_out.write('{:<8} {:>8.6f} {:>8.6f} {:>8.4f} {:>8.4f} {:>8.4f} \
{:>8.4f} {:>8.4f} {:>6}'.format(*l))
            f_out.write('\n')
        # ^ Stars selected by the reduced membership function.
        for line in red_memb_no_fit:
            l = line + ['0']
            f_out.write('{:<8} {:>8.6f} {:>8.6f} {:>8.4f} {:>8.4f} {:>8.4f} \
{:>8.4f} {:>8.4f} {:>6}'.format(*l))
            f_out.write('\n')

    print('Membership probabilities saved to file.')
