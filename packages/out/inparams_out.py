
from os.path import isfile
from time import strftime


def main(npd, file_end, **kwargs):
    """
    Save input 'params_input.dat' file as output file.
    """

    fname = './params_input.dat'
    if isfile('./params_input' + file_end + '.dat'):
        fname = './params_input' + file_end + '.dat'

    with open(fname, 'r') as f:
        # Read file into data var.
        data = f.readlines()

    # Add current date & time
    now_time = strftime("%Y-%m-%d %H:%M:%S")
    data[2] = "#                      Created: [{}]\n".format(now_time)

    # Write
    with open(npd['params_out'], 'w') as f:
        f.writelines(data)
