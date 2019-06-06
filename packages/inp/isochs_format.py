

def cmd_age_format():
    """
    Define reg expression to isolate the age of a CMD isochrone from its
    commented title line.
    """
    age_format = r"Age = (.+?) yr"
    return age_format


def cmd_line_start_format(cmd_evol_tracks, evol_track):
    """
    Return the format of the line where column names are read from, for each
    set of evolutionary tracks.
    """
    if evol_track in cmd_evol_tracks.keys():
        # String that identifies the beginning of a new isochrone.
        line_start = "# Zini"
    else:
        print("ERROR: evolutionary track not recognized.")

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
                # Split into column names.
                ls = line.split()
                # Remove comment character.
                if ls[0] == '#':
                    del ls[0]
                else:
                    print("  WARNING: theoretical isochrones are not\n"
                          "  formatted as expected!\n")
                break

    return ls


def cmd_common_ids(evol_track, l):
    """
    These parameters are equivalent across photometric systems in a given
    set of CMD evolutionary tracks.
    """
    # Get indexes of initial mass, magnitudes, etc.
    Mini, int_IMF, Mass, logL, logTe, logg, label, mbolmag = l.index('Mini'),\
        l.index('int_IMF'), l.index('Mass'), l.index('logL'),\
        l.index('logTe'), l.index('logg'), l.index('label'),\
        l.index('mbolmag')

    # Simple check
    if evol_track in ['PAR10', 'PAR11', 'PAR12']:
        if [Mini, int_IMF, Mass, logL, logTe, logg, label, mbolmag] !=\
                [2, 3, 4, 5, 6, 7, 8, 9]:
            print("  WARNING: extra parameters in isochrones are not\n"
                  "  positioned as expected!\n")

    # TODO these extra params consume a lot of memory and are not used for
    # now. Don't read until they are used or I find a more efficient method
    # of storing them.
    return [Mini]  #, int_IMF, Mass, logL, logTe, logg, label, mbolmag]


def girardi_filters_ids(l, uniq_fltrs):
    '''
    Read id numbers for each filter defined in this photometric system.
    '''
    filters_ids = []
    for filt in uniq_fltrs:
        filters_ids.append(l.index(filt))

    return filters_ids
