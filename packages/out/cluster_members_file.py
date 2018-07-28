
from astropy.io import ascii


def main(clp, memb_file_out, **kwargs):
    '''
    Create output data file with stars inside the cluster radius along with
    their membership probabilities.

    TODO write all data out
    '''

    cl_reg_fit, cl_reg_no_fit = clp['cl_reg_fit'], clp['cl_reg_no_fit']

    data = []
    idx = ['1', '0']
    for i, reg in enumerate([cl_reg_fit, cl_reg_no_fit]):
        for line in reg:
            # Identify stars selected by the removal function.
            data.append(
                line[:3] + [line[3][0]] + [line[4][0]] + [line[5][0]] +
                [line[6][0]] + [line[9]] + [idx[i]])

    ascii.write(
        list(zip(*data)), memb_file_out, format='csv', overwrite=True,
        names=['ID', 'x', 'y', 'mag1', 'e_mag', 'col1', 'e_col1', 'MP', 'sel'])

    print('Cluster region saved to file.')
