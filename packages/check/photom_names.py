
import sys
from os.path import join


def phot_syst_filt_check(all_systs, entry, phot_syst, filter_name):
    """
    Check photometric system and filter names.
    """
    if phot_syst not in all_systs.keys():
        sys.exit("\nERROR: unknown photometric system '{}', given in"
                 "'{}'.".format(phot_syst, entry))
    if filter_name not in all_systs[phot_syst][1]:
        sys.exit("\nERROR: filter '{}' given in '{}' is not present\n"
                 "in '{}' photometric system.".format(
                     filter_name, entry, all_systs[phot_syst][0]))


def check(mypath, pd):
    """
    Check that the magnitudes and colors are properly defined. Store data to
    read from cluster's file, to load the correct set of isochrones, and to
    properly generate the synthetic clusters (if the best match function is
    set to run).
    """

    # Read column indexes for the IDs and the coordinates, and the coordinate
    # format used (pixels or degrees)
    try:
        id_indx, x_indx, y_indx, coords = pd['id_coords']
        id_indx, x_indx, y_indx = map(int, [id_indx, x_indx, y_indx])
        coords = str(coords)
    except ValueError:
        sys.exit("ERROR: bad format for the IDs and coordinates\n"
                 "column indexes in 'params_input.dat'.")

    # Check px/deg.
    if coords not in ('px', 'deg'):
        sys.exit("ERROR: coordinate units '{}' given in the input parameters\n"
                 "file are incorrect.".format(coords))

    # Dictionary of photometric systems defined in the CMD service.
    all_systs = pd['cmd_systs']

    # Extract magnitudes (filters) data.
    mag_indx, e_mag_indx, filters = [], [], []
    for mag in pd['id_mags'][0::2]:
        try:
            colum_indx, phot_syst, filter_name = mag.split(',')
        except ValueError:
            sys.exit("ERROR: bad formatting for filter '{}'".format(mag))
        # Used to read data from cluster file in 'get_data.
        mag_indx.append(int(colum_indx))
        if pd['bf_flag']:
            # Check.
            phot_syst_filt_check(all_systs, mag, phot_syst, filter_name)
        # Name of photometric system and filter, used to extract its
        # synthetic data from the correct theoretical isochrone.
        filters.append((phot_syst, filter_name))
    # Extract magnitude error columns.
    if len(pd['id_mags'][1::2]) == len(pd['id_mags'][0::2]):
        for colum_indx in pd['id_mags'][1::2]:
            e_mag_indx.append(int(colum_indx))
    elif len(pd['id_mags'][1::2]) < len(pd['id_mags'][0::2]):
        sys.exit("ERROR: missing error column index for filter"
                 " in 'params_input dat'.")

    # Extract colors data.
    col_indx, e_col_indx, c_filters, colors = [], [], [], []
    for col in pd['id_cols'][0::2]:
        try:
            colum_indx, phot_syst, filter_name1, filter_name2 = col.split(',')
        except ValueError:
            sys.exit("ERROR: bad formatting for color '{}'".format(col))
        col_indx.append(int(colum_indx))
        if pd['bf_flag']:
            # Check.
            phot_syst_filt_check(all_systs, col, phot_syst, filter_name1)
            phot_syst_filt_check(all_systs, col, phot_syst, filter_name2)
        c_filters.append((phot_syst, filter_name1))
        c_filters.append((phot_syst, filter_name2))
        colors.append((phot_syst, filter_name1 + ',' + filter_name2))
    # Extract colors error columns.
    if len(pd['id_cols'][1::2]) == len(pd['id_cols'][0::2]):
        for colum_indx in pd['id_cols'][1::2]:
            e_col_indx.append(int(colum_indx))
    elif len(pd['id_cols'][1::2]) < len(pd['id_cols'][0::2]):
        sys.exit("ERROR: missing error column index for color"
                 " in 'params_input dat'.")

    all_syst_filters, iso_paths = [], []
    if pd['bf_flag']:
        # Remove duplicate filters (if they exist), and combine them into one
        # tuple per photometric system.
        # The resulting list looks like this:
        # [('2', 'T1', 'C'), ('4', 'B', 'V'), ('65', 'J')]
        # where the first element of each tuple points to the photometric
        # system, and the remaining elements are the unique filters in that
        # system.
        all_syst_filters = list(set(filters + c_filters))
        d = {}
        for k, v in all_syst_filters:
            d.setdefault(k, [k]).append(v)
        all_syst_filters = sorted(map(tuple, d.values()))

        # Fix isochrones location according to the CMD and set selected.
        text1 = pd['cmd_evol_tracks'][pd['evol_track']][0]
        # Generate correct name for the isochrones path.
        iso_paths = []
        for p_syst in all_syst_filters:
            text2 = all_systs[p_syst[0]][0]
            # Set iso_path according to the above values.
            iso_paths.append(
                join(mypath + 'isochrones/' + text1 + '_' + text2))

        # REMOVE when support for multiple photometric systems is in place.
        if len(all_syst_filters) > 1:
            sys.exit("ERROR: more than one photometric system defined.")

    # Add data to parameters dictionary.
    pd['id_indx'], pd['x_indx'], pd['y_indx'], pd['coords'],\
        pd['mag_indx'], pd['e_mag_indx'], pd['filters'], pd['col_indx'],\
        pd['e_col_indx'], pd['colors'], pd['all_syst_filters'],\
        pd['iso_paths'] = id_indx, x_indx, y_indx, coords, mag_indx,\
        e_mag_indx, filters, col_indx, e_col_indx, colors,\
        all_syst_filters, iso_paths

    return pd
