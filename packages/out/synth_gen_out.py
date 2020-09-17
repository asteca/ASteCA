
import os
import numpy as np
from astropy.io import ascii
from astropy.table import Table


def createFile(
    filters, colors, extra_pars, model, synth_clust, sigma, x_cl, y_cl, x_fl,
        y_fl, synth_field, sigma_field, CI, rc, rt, name_gen_file):
    """
    """
    cl_ids = ["1" + str(_) for _ in np.arange(1, x_cl.size + 1)]
    fl_ids = ["2" + str(_) for _ in np.arange(1, x_fl.size + 1)]
    IDs = np.array([cl_ids + fl_ids]).astype(float)
    coords = np.array([
        x_cl.tolist() + x_fl.tolist(), y_cl.tolist() + y_fl.tolist()])
    if synth_field.any():
        photom = np.concatenate((synth_clust, synth_field), 1)
        e_photom = np.concatenate((sigma, sigma_field), 1)
    else:
        photom, e_photom = synth_clust, sigma

    extra_pars_fl = np.zeros((extra_pars.shape[0], x_fl.size))
    extra_pars_all = np.concatenate((extra_pars, extra_pars_fl), 1)

    data = np.concatenate((
        IDs, coords, photom, e_photom, extra_pars_all), 0)

    # Columns names
    mag_cols = [f[1].replace('mag', '') for f in filters]
    mag_cols += [c[1].replace(',', '-').replace('mag', '') for c in colors]
    mag_cols += ['e_' + f[1].replace('mag', '') for f in filters]
    mag_cols += [
        'e_' + c[1].replace(',', '-').replace('mag', '') for c in colors]
    names = ['ID', 'x', 'y'] + mag_cols + ['bfr_prob', 'bfr_mass', 'm_ini']

    # Same format for all columns
    dict_formats = {_: '%10.4f' for _ in names[1:]}
    dict_formats['ID'] = '%10.0f'

    t = Table(data.T, names=names)

    header = "#\n# CI, rc, rt : {:.2f}, {:.3f}, {:.3f}".format(
        CI, rc, rt) + "\n#\n# model (z, log(age), E_BV, dm, M, b_fr) : " +\
        "{:.5f} {:.4f} {:.3f} {:.3f} {:.0f} {:.2f}\n".format(*model) +\
        "#\n"
    with open(name_gen_file, mode='w') as f:
        f.write(header)

    with open(name_gen_file, mode='a') as f:
        f.seek(0, os.SEEK_END)
        ascii.write(
            t, f, format='csv', fast_writer=True, formats=dict_formats)
