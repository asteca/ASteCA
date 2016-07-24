
from ..inp import input_params as g


def girardi_age_format():
    """
    Define reg expression to isolate the age of an isochrone.
    """
    age_format = r"Age = \t(.+?) yr"
    return age_format


def main(met_f):
    '''
    Read line start format and columns indexes for the selected set of
    Girardi isochrones and chosen CMD.
    '''

    cmd_select, iso_select = g.ps_params[1], g.ps_params[2]

    # Assign values according to the system and set of isochrones selected.
    if iso_select in ['PAR10', 'PAR11', 'PAR12', 'PAR12C']:
        # String that identifies the beginning of a new isochrone.
        line_start = "#\tIsochrone  Z = "
    elif iso_select in ['GIR02', 'MAR08', 'MAR08B', 'MAR08A']:
        # String that identifies the beginning of a new isochrone.
        line_start = "#\tIsochrone\tZ = "

    # Open the metallicity file.
    with open(met_f, mode="r") as f_iso:
        # Iterate through each line in the file.
        for line in f_iso:
            # When this line is found, extract columns names from the
            # following line.
            if line.startswith(line_start):
                # Line found, split into column names.
                l = f_iso.next().split()
                # Remove first element from list, the comment character.
                if l[0] == '#':
                    del l[0]
                else:
                    print("  WARNING: theoretical isochrones are not"
                          " formatted\nas expected.")
                break

    # Get indexes of initial mass, magnitudes, etc.
    mass = l.index('M_ini')
    if cmd_select == 1:
        # V, B
        mag1, mag2 = l.index('V'), l.index('B')
    elif cmd_select == 2:
        # V, I
        mag1, mag2 = l.index('V'), l.index('I')
    if cmd_select == 3:
        # V, U
        mag1, mag2 = l.index('V'), l.index('U')
    elif cmd_select == 4:
        # T1, C
        mag1, mag2 = l.index('T1'), l.index('C')
    elif cmd_select == 5:
        # J, H
        mag1, mag2 = l.index('J'), l.index('H')
    elif cmd_select == 6:
        # H, J
        mag1, mag2 = l.index('H'), l.index('J')
    elif cmd_select == 7:
        # K_s, H
        mag1, mag2 = l.index('K_s'), l.index('H')
    elif cmd_select == 8:
        # g, u
        mag1, mag2 = l.index('g'), l.index('u')
    elif cmd_select == 9:
        # g, r
        mag1, mag2 = l.index('g'), l.index('r')
    elif cmd_select == 10:
        # y, b
        mag1, mag2 = l.index('y'), l.index('b')
    elif cmd_select == 11:
        # y, v
        mag1, mag2 = l.index('y'), l.index('v')
    elif cmd_select == 12:
        # y, u
        mag1, mag2 = l.index('y'), l.index('u')
    elif cmd_select == 13:
        # F606W, F814W
        mag1, mag2 = l.index('F606W'), l.index('F814W')

    # Get age format.
    age_format = girardi_age_format()

    return line_start, age_format, mass, mag1, mag2

# # Mass column.
# mass = 2
# if cmd_select == 1:
#     # V, B
#     mag1, mag2 = 10, 9
# elif cmd_select == 2:
#     # V, I
#     mag1, mag2 = 10, 12
# if cmd_select == 3:
#     # V, U
#     mag1, mag2 = 10, 8
# elif cmd_select == 4:
#     # T1, C
#     mag1, mag2 = 10, 8
# elif cmd_select == 5:
#     # J, H
#     mag1, mag2 = 8, 9
# elif cmd_select == 6:
#     # H, J
#     mag1, mag2 = 9, 8
# elif cmd_select == 7:
#     # K_s, H
#     mag1, mag2 = 10, 9
# elif cmd_select == 8:
#     # g, u
#     mag1, mag2 = 9, 8
# elif cmd_select == 9:
#     # g, r
#     mag1, mag2 = 9, 10
# elif cmd_select == 10:
#     # y, b
#     mag1, mag2 = 12, 11
# elif cmd_select == 11:
#     # y, v
#     mag1, mag2 = 12, 10
# elif cmd_select == 12:
#     # y, u
#     mag1, mag2 = 12, 9
# elif cmd_select == 13:
#     # F606W, F814W
#     mag1, mag2 = 13, 18