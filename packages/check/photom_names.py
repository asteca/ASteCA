
import sys
from os.path import join
import CMD_phot_systs_filts


def phot_syst_filt_check(entry, phot_syst, filter_name):
    """
    Check photometric system and filter names.
    """
    all_systs = CMD_phot_systs_filts.main()
    if phot_syst not in all_systs.keys():
        sys.exit("\nERROR: unknown photometric system '{}', given in"
                 "'{}'.".format(phot_syst, entry))
    if filter_name not in all_systs[phot_syst][1]:
        sys.exit("\nERROR: filter '{}' given in '{}' is not present\n"
                 "in '{}' photometric system,\n.".format(
                     filter_name, entry, all_systs[phot_syst][0]))


def check(mypath, pd):
    """
    Check that the magnitudes and colors are properly defined. Store data to
    read from cluster's file, to load the correct set of isochrones, and to
    properly generate the synthetic clusters (if the best match function is
    set to run).
    """

    # Extract magnitudes (filters) data.
    mag_columns, filters = [], []
    for mag in pd['id_mags'][0::2]:
        try:
            colum_indx, phot_syst, filter_name = mag.split(',')
        except:
            sys.exit("ERROR: bad formatting for filter '{}'".format(mag))
        # Used to read data from cluster file in 'get_data.
        mag_columns.append(int(colum_indx))
        # Check.
        phot_syst_filt_check(mag, phot_syst, filter_name)
        # Name of photometric system and filter, used to extract its
        # synthetic data from the correct theoretical isochrone.
        filters.append((int(phot_syst), filter_name))
    # Store error column number.
    e_mag_columns = []
    for e_mag_idx in pd['id_mags'][1::2]:
            e_mag_columns.append(int(e_mag_idx))

    # Extract colors data.
    col_columns, c_filters, colors = [], [], []
    for col in pd['id_cols'][0::2]:
        try:
            colum_indx, phot_syst, color = col.split(',')
        except:
            sys.exit("ERROR: bad formatting for color '{}'".format(col))
        col_columns.append(int(colum_indx))
        # Separate filters from color name.
        try:
            filter_name1, filter_name2 = color.split('-')
        except:
            sys.exit("ERROR: bad color name '{}'".format(color))
        # Check.
        phot_syst_filt_check(col, phot_syst, filter_name1)
        phot_syst_filt_check(col, phot_syst, filter_name2)
        c_filters.append((int(phot_syst), filter_name1))
        c_filters.append((int(phot_syst), filter_name2))
        colors.append((int(phot_syst), color))
    e_col_columns = []
    for e_col_idx in pd['id_cols'][1::2]:
        e_col_columns.append(int(e_col_idx))

    # Remove duplicate filters (if they exist), and combine them into one
    # tuple per photometric system.
    # The resulting list looks like this:
    # [(0, 'B', 'V'), (1, 'J'), (2, 'T1', 'C')]
    # where the first element of each tuple points to the photometric system,
    # and the remaining elements are the unique filters in that system.
    all_filters = list(set(filters + c_filters))
    d = {}
    for k, v in all_filters:
        d.setdefault(k, [k]).append(v)
    all_filters = sorted(map(tuple, d.values()))

    print mag_columns, e_mag_columns, filters
    print col_columns, e_col_columns, colors
    print all_filters
    import pdb; pdb.set_trace()  # breakpoint 0f6792af //


    # Fix isochrones location according to the CMD and set selected.
    text1, text2 = 'none', 'none'
    # Map isochrones set selection to proper name.
    iso_sys = {'PAR12C': 'parsec12C', 'PAR12': 'parsec12',
               'PAR10': 'parsec10', 'PAR11': 'parsec11',
               'MAR08A': 'marigo08A', 'MAR08B': 'marigo08B',
               'MAR08': 'marigo08', 'GIR02': 'girardi02'}
    text1 = iso_sys.get(iso_select)
    # Generate correct name for the isochrones path.
    if cmd_select in {1, 2, 3}:
        text2 = 'ubvrijhk'
    elif cmd_select in {4}:
        text2 = 'washington'
    elif cmd_select in {5, 6, 7}:
        text2 = '2mass'
    elif cmd_select in {8, 9}:
        text2 = 'sloan'
    elif cmd_select in {10, 11, 12}:
        text2 = 'stroemgren'
    elif cmd_select in {13}:
        text2 = 'acs_wfc'
    # Set iso_path according to the above values.
    iso_path = join(mypath + '/isochrones/' + text1 + '_' + text2)
    pd['iso_path'] = iso_path

    return pd
