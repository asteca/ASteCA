

def check(mypath, pd):
    """
    Check that the magnitudes and colors are properly defined.
    Store data to read from cluster's file, to load the correct set of
    isochrones, and to properly generate the synthetic clusters (if the best
    match function is set to run).
    """

    # Check read mode.
    if pd['read_mode'] not in pd['read_mode_accpt']:
        raise ValueError("read mode '{}' given in the input parameters\n"
                         "file is incorrect.".format(pd['read_mode']))

    # Check px/deg.
    if pd['coords'] not in pd['coord_accpt']:
        raise ValueError(
            "coordinate units '{}' given in the input parameters\n"
            "file are incorrect.".format(pd['coords']))

    # Read column indexes for the IDs and the coordinates.
    if pd['read_mode'] == 'num':
        # Name of columns when no header is present. The '+ 1' is because
        # astropy Tables' first column is named 'col1', not 'col0'.
        id_col, x_col, y_col = [
            'col' + str(int(i) + 1) for i in [
                pd['id_ids'], pd['id_xdata'], pd['id_ydata']]]
    else:
        id_col, x_col, y_col = pd['id_ids'], pd['id_xdata'], pd['id_ydata']

    # Check error values
    N_colors = int(len(pd['id_cols']) / 2.)
    for err in pd['err_max']:
        try:
            if float(err) < 0.:
                raise ValueError("max error value must be in the range >0.")
        except ValueError:
            # This assumes that there is no maximum number of colors that can
            # be defined
            if len(pd['err_max']) - 4 != N_colors:
                raise ValueError(
                    "there are {} 'e_*_max' values defined,\n"
                    "there should be {} ({} color(s) defined).".format(
                        len(pd['err_max']), 4 + N_colors, N_colors))

    # Dictionary of photometric systems defined in the CMD service.
    all_systs = pd['cmd_systs']

    # Extract magnitudes (filters) data.
    mag_col, e_mag_col, filters = [], [], []
    for mag in pd['id_mags'][0::2]:
        try:
            column_id, phot_syst, filter_name = mag.split(',')
        except ValueError:
            raise ValueError("bad formatting for filter '{}'".format(mag))
        # Used to read data from cluster file in 'get_data.
        mag_col += [
            'col' + str(int(column_id) + 1) if pd['read_mode'] == 'num' else
            column_id]
        if pd['bf_flag']:
            # Check.
            phot_syst_filt_check(all_systs, mag, phot_syst, filter_name)
        # Name of photometric system and filter, used to extract its
        # synthetic data from the correct theoretical isochrone.
        filters.append((phot_syst, filter_name))
    # Extract magnitude error columns.
    if len(pd['id_mags'][1::2]) == len(pd['id_mags'][0::2]):
        for column_id in pd['id_mags'][1::2]:
            e_mag_col += [
                'col' + str(int(column_id) + 1) if pd['read_mode'] == 'num'
                else column_id]
    elif len(pd['id_mags'][1::2]) < len(pd['id_mags'][0::2]):
        raise ValueError("missing error column name/index for filter"
                         " in 'params_input dat'.")

    # Extract colors data.
    col_col, e_col_col, c_filters, colors = [], [], [], []
    for col in pd['id_cols'][0::2]:
        try:
            column_id, phot_syst, filter_name1, filter_name2 = col.split(',')
        except ValueError:
            raise ValueError("bad formatting for color '{}'".format(col))
        col_col += [
            'col' + str(int(column_id) + 1) if pd['read_mode'] == 'num'
            else column_id]
        if pd['bf_flag']:
            # Check.
            phot_syst_filt_check(all_systs, col, phot_syst, filter_name1)
            phot_syst_filt_check(all_systs, col, phot_syst, filter_name2)
        c_filters.append((phot_syst, filter_name1))
        c_filters.append((phot_syst, filter_name2))
        colors.append((phot_syst, filter_name1 + ',' + filter_name2))
    # Extract colors error columns.
    if len(pd['id_cols'][1::2]) == len(pd['id_cols'][0::2]):
        for column_id in pd['id_cols'][1::2]:
            e_col_col += [
                'col' + str(int(column_id) + 1) if pd['read_mode'] == 'num'
                else column_id]
    elif len(pd['id_cols'][1::2]) < len(pd['id_cols'][0::2]):
        raise ValueError("missing error column name/index for color"
                         " in 'params_input dat'.")

    # Add data to parameters dictionary.
    pd['id_col'], pd['x_col'], pd['y_col'],\
        pd['mag_col'], pd['e_mag_col'], pd['filters'], pd['col_col'],\
        pd['e_col_col'], pd['colors'], pd['c_filters'] = id_col, x_col, y_col,\
        mag_col, e_mag_col, filters, col_col, e_col_col, colors, c_filters

    # Check max error values.
    for i, e in enumerate(pd['err_max']):
        try:
            float(e)
        except ValueError:
            if e not in ('n', 'N'):
                raise ValueError("bad value ('{}') for '{}' max error"
                                 " in 'params_input.dat'".format(e, i))

    return pd


def phot_syst_filt_check(all_systs, entry, phot_syst, filter_name):
    """
    Check photometric system and filter names.
    """
    if phot_syst not in all_systs.keys():
        raise ValueError("unknown photometric system '{}', given in"
                         "'{}'.".format(phot_syst, entry))
    if filter_name not in all_systs[phot_syst][1]:
        raise ValueError("filter '{}' given in '{}' is not present\n"
                         "in '{}' photometric system.".format(
                             filter_name, entry, all_systs[phot_syst][0]))
