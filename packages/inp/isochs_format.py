

def cmd_age_format():
    """
    Define reg expression to isolate the age of an isochrone.
    """
    age_format = r"Age = \t(.+?) yr"
    return age_format


def cmd_line_start_format(evol_track):
    """
    Return the format of the line where column names are read from, for each
    set of evolutionary tracks.
    """
    if evol_track in ['PAR10', 'PAR11', 'PAR12', 'PAR12C']:
        # String that identifies the beginning of a new isochrone.
        line_start = "#\tIsochrone  Z = "
    elif evol_track in ['GIR02', 'MAR08', 'MAR08B', 'MAR08A']:
        # String that identifies the beginning of a new isochrone.
        line_start = "#\tIsochrone\tZ = "

    return line_start


def read_line_start(met_f, line_start):
    """
    Read the line where the filter names are positioned, for each set of
    evolutionary tracks.
    """
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
                    print("  WARNING: theoretical isochrones are not\n"
                          "  formatted as expected!")
                break

    return l


def cmd_common_ids(l):
    """
    These parameters are equivalent across photometric systems in a given
    set of CMD evolutionary tracks. Their column numbers should also be equal
    (as far as I checked!)
    """
    # Get indexes of initial mass, magnitudes, etc.
    M_ini, M_act, logLLo, logTe, logG, mbol = l.index('M_ini'),\
        l.index('M_act'), l.index('logL/Lo'), l.index('logTe'),\
        l.index('logG'), l.index('mbol')

    if [M_ini, M_act, logLLo, logTe, logG, mbol] != [2, 3, 4, 5, 6, 7]:
        print("  WARNING: extra parameters in isochrones are not\n"
              "  positioned as expected!")

    return [M_ini, M_act, logLLo, logTe, logG, mbol]


def girardi_filters_ids(l, uniq_fltrs):
    '''
    Read id numbers for each filter defined in this photometric system.
    '''
    filters_ids = []
    for filt in uniq_fltrs:
        filters_ids.append(l.index(filt))

    return filters_ids
