
import numpy as np
from astropy.table import Table
from ..inp import data_IO


def main(clp, npd, pd, td):
    """
    1. Create output data file with stars in the best fit synthetic cluster
       found by the 'Best Fit' function.
    2. Create output data file for the observed cluster with individual
       estimated masses
    """

    if pd['best_fit_algor'] != 'n':
        # If cluster is not empty.
        if not clp['synth_cl_phot'].any():
            print("  ERROR: empty synthetic cluster could not\n"
                  "  be saved to file")
            return

        writeFileOut(
            npd, pd['filters'], pd['colors'], td['m_ini_idx'],
            clp['synth_cl_phot'], clp['synth_cl_sigma'], clp['cl_max_mag'],
            clp['st_mass_mean'], clp['st_mass_std'], clp['st_mass_mean_binar'],
            clp['st_mass_std_binar'], clp['prob_binar'])
        print("Synthetic and observed clusters saved to file")


def writeFileOut(
    npd, filters, colors, m_ini_idx, synth_clust, sigma, cl_max_mag,
        st_mean, st_std, st_mean_binar, st_std_binar, prob_binar):
    """
    Save best fit synthetic cluster found, and per star masses to file.
    """

    # Write synthetic cluster file
    e_mags_cols = np.array(sigma).T
    # Create IDs identifying binary systems
    binar_ID = []
    for i, bi in enumerate(synth_clust[-1]):
        if bi == -99.:
            binar_ID.append('1' + str(i))
        else:
            binar_ID.append('2' + str(i))
    phot_col = [f[1] for f in filters] +\
        ['(' + c[1].replace(',', '-') + ')' for c in colors]
    ephot_col = ['e_' + f[1] for f in filters] +\
        ['e_(' + c[1].replace(',', '-') + ')' for c in colors]
    col_names = ['ID'] + phot_col + ephot_col + ['M_ini1', 'M_ini2']
    data = [binar_ID] + [_ for _ in synth_clust[:m_ini_idx]] +\
        [_ for _ in e_mags_cols.T] + [synth_clust[m_ini_idx], synth_clust[-1]]
    synth_table = Table(data, names=col_names)
    synth_table.meta['comments'] = ["Binary systems ID's begin with a '2'"]
    data_IO.dataSave(synth_table, npd['synth_file_out'], 'w')

    # Write observed cluster with masses to file
    if st_mean.any():
        st_ID = np.array(list(zip(*cl_max_mag))[0])
        main_mag = np.array(list(zip(*cl_max_mag))[3]).T[0]
        first_col = np.array(list(zip(*cl_max_mag))[5]).T[0]
        mass_table = Table(
            [st_ID, main_mag, first_col, st_mean, st_std, st_mean_binar,
             st_std_binar, prob_binar],
            names=['ID', 'Mag', 'Col', 'M1', 'M1_std', 'M2', 'M2_std',
                   'P_binar'])
        mass_table.meta['comments'] = [
            '', 'Subset of stars selected as members', 'M1: primary mass',
            'M2: secondary mass (if binary system)',
            'P_binar: binary probability', '']
        data_IO.dataSave(
            mass_table, npd['mass_file_out'], 'w',
            {'Mag': '%12.5f', 'Col': '%12.5f', 'M1': '%12.5f',
             'M1_std': '%12.5f', 'M2': '%12.5f', 'M2_std': '%12.5f',
             'P_binar': '%12.2f'})
    else:
        print("  WARNING: could not save masses/binary probabilities to file")
