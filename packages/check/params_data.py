

def check(mypath, pd):
    """
    Check that the magnitudes and colors are properly defined.
    Store data to read from cluster's file, to load the correct set of
    isochrones, and to properly generate the synthetic clusters (if the best
    match function is set to run).
    """
    # Read column indexes for the IDs and the coordinates.
    id_col, x_col, y_col = pd['id_ids'], pd['id_xdata'], pd['id_ydata']

    # Extract magnitudes (filters) data.
    try:
        # column_id, phot_syst, filter_name = mag.split(',')
        phot_syst, filter_name, mag_col, e_mag_col = pd['id_mags']
        mag_col, e_mag_col = [mag_col], [e_mag_col]
    except ValueError:
        raise ValueError("Bad formatting for the magnitude")
    # Name of photometric system and filter, used to extract its
    # synthetic data from the correct theoretical isochrone.
    filters = [(phot_syst, filter_name)]

    # Extract colors data.
    col_col, e_col_col, c_filters, colors = [], [], [], []
    for color in pd['id_cols']:
        try:
            phot_syst, filter1, filter2, color_id, e_color_id = color
        except ValueError:
            raise ValueError("Bad formatting for (one of) the color(s)")
        col_col += [color_id]
        e_col_col += [e_color_id]

        c_filters.append((phot_syst, filter1))
        c_filters.append((phot_syst, filter2))
        colors.append((phot_syst, filter1 + ',' + filter2))

    if not colors:
        raise ValueError("At least one color must be defined")

    # Add data to parameters dictionary.
    pd['id_col'], pd['x_col'], pd['y_col'],\
        pd['mag_col'], pd['e_mag_col'], pd['filters'], pd['col_col'],\
        pd['e_col_col'], pd['colors'], pd['c_filters'] = id_col, x_col, y_col,\
        mag_col, e_mag_col, filters, col_col, e_col_col, colors, c_filters

    # DEPRECATED 23/03/22
    # # Check error values
    # for err in pd['err_max']:
    #     try:
    #         ferr = float(err)
    #     except ValueError:
    #         if err != 'max':
    #             raise ValueError("bad value ('{}') for max error"
    #                              " in 'asteca.ini'".format(err))
    #         ferr = 1.
    #     if ferr < 0.:
    #         raise ValueError("max error value must be positive")

    # # This assumes that there is no maximum number of colors that can
    # # be defined
    # N_colors = len(pd['id_cols'])
    # if len(pd['err_max']) - 4 != N_colors:
    #     raise ValueError(
    #         "there are {} 'e_*_max' values defined, there should\n"
    #         "be {}, because {} color(s) is/are defined.".format(
    #             len(pd['err_max']), 4 + N_colors, N_colors))

    return pd
