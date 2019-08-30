
from time import strftime


def main(npd):
    """
    Save input 'params_input.dat' file as output file.
    """
    # Read
    with open('./params_input.dat', 'r') as f:
        # Read file into data var.
        data = f.readlines()

    # Add current date & time
    now_time = strftime("%Y-%m-%d %H:%M:%S")
    data[2] = "#                      Created: [{}]\n".format(now_time)

    # Write
    with open(npd['params_out'], 'w') as f:
        f.writelines(data)
