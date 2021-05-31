
from os.path import isfile
from time import strftime


def main(npd, file_end, **kwargs):
    """
    Save input 'asteca.ini' file as output file.
    """

    fname = './asteca.ini'
    if isfile('./asteca' + file_end + '.ini'):
        fname = './asteca' + file_end + '.ini'

    with open(fname, 'r') as f:
        # Read file into data var.
        data = f.readlines()

    # Add current date & time
    now_time = strftime("%Y-%m-%d %H:%M:%S")
    data[2] = "#                      Created: [{}]\n".format(now_time)

    # Write
    with open(npd['params_out'], 'w') as f:
        f.writelines(data)
