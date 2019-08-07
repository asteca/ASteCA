
from astropy.io import ascii
from astropy.table import Table


def main(
    clp, npd, id_col, x_col, y_col, mag_col, e_mag_col, col_col, e_col_col,
    plx_col, e_plx_col, pmx_col, e_pmx_col, pmy_col, e_pmy_col, rv_col,
        e_rv_col, coords, read_mode, filters, colors, **kwargs):
    '''
    Create output data file with stars inside the cluster radius along with
    their membership probabilities and 'clean region' selection identifier.
    '''
    # Generate column names in case 'read_mode=num' was used.
    if read_mode == 'num':
        clp['col_names_keep'][0] = 'ID'
        if coords == 'deg':
            clp['col_names_keep'][1] = 'ra'
            clp['col_names_keep'][2] = 'dec'
        else:
            clp['col_names_keep'][1] = 'x'
            clp['col_names_keep'][2] = 'y'
        clp['col_names_keep'][3] = filters[0][1]
        for i, c in enumerate(col_col):
            j = clp['col_names_keep'].index(c)
            clp['col_names_keep'][j] = colors[i][1].replace(',', '')
        for i, ec in enumerate(e_col_col):
            j = clp['col_names_keep'].index(ec)
            clp['col_names_keep'][j] = 'e' + colors[i][1].replace(',', '')
        clp['col_names_keep'][4] = 'e' + filters[0][1]
        cnames = ('Plx', 'ePlx', 'PMx', 'ePMx', 'PMy', 'ePMy', 'RV', 'eRV')
        for i, col in enumerate((
                plx_col, e_plx_col, pmx_col, e_pmx_col, pmy_col, e_pmy_col,
                rv_col, e_rv_col)):
            if col is not False:
                j = clp['col_names_keep'].index(col)
                clp['col_names_keep'][j] = cnames[i]

    # Add ID associated to the use of the each star in the fundamental
    # parameters estimation process (ie: after cleaning the cluster region).
    kmsk = []
    for _ in (plx_col, pmx_col, pmy_col, rv_col):
        kmsk.append(True if _ is not False else False)
    data, idx = [], ['1', '0']
    for i, reg in enumerate([clp['cl_reg_fit'], clp['cl_reg_no_fit']]):
        for st in reg:
            mag, emag = st[3].tolist(), st[4].tolist()
            cols, ecols = st[5].tolist(), st[6].tolist()
            # Filter out 'nan's
            kinem, ekinem = st[7][kmsk].tolist(), st[8][kmsk].tolist()
            # Identify stars selected by the removal function.
            data.append(
                st[:3] + mag + emag + cols + ecols + kinem + ekinem +
                [st[9], idx[i]])

    # Add "incomplete" data in cluster region to file.
    ids_data = list(zip(*data))[0]
    for i, st in enumerate(clp['cl_region_i']):
        if st[0] not in ids_data:
            mag, emag = st[3].tolist(), st[4].tolist()
            cols, ecols = st[5].tolist(), st[6].tolist()
            # Filter out 'nan's
            kinem, ekinem = st[7][kmsk].tolist(), st[8][kmsk].tolist()
            # Identify stars selected by the removal function.
            data.append(
                st[:3] + mag + emag + cols + ecols + kinem + ekinem +
                [round(clp['memb_probs_cl_region_i'][i], 2), '-1'])

    t = Table(list(zip(*data)), names=clp['col_names_keep'] + ['MP', 'sel'])
    ascii.write(
        t, npd['memb_file_out'], overwrite=True, format='csv',
        fast_writer=False  # <-- TODO remove when the bug is fixed
    )

    print('Cluster region and MPs saved to file.')
