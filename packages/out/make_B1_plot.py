
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
from . import mp_errors
from . import add_version_plot
from . import prep_plots
from . prep_plots import figsize_x, figsize_y, grid_x, grid_y


def main(
    npd, cld, pd, err_lst, cl_region, stars_out, col_0_comb,
        mag_0_comb, **kwargs):
    """
    Make B1 block plots.
    """
    fig = plt.figure(figsize=(figsize_x, figsize_y))
    gs = gridspec.GridSpec(grid_y, grid_x)
    add_version_plot.main(y_fix=1.)

    # Obtain plotting parameters and data.
    err_bar_all = prep_plots.error_bars(
        cld['mags'][0], np.nan, err_lst, 'all')

    # Magnitude vs uncertainties diagrams.
    arglist = [
        [gs, pd['colors'], pd['filters'], pd['id_kinem'], cld['mags'],
         cl_region, stars_out, err_bar_all]
    ]
    for n, args in enumerate(arglist):
        mp_errors.plot(n, *args)

    plt.suptitle(
        (r"$N={}$ (cluster + field regions)").format(
            len(cl_region) + len(stars_out)),
        x=.26, y=1.005)

    # Generate output file.
    fig.tight_layout()
    plt.savefig(
        join(npd['output_subdir'], str(npd['clust_name']) + '_B1'
             + npd['ext']))
    # Close to release memory.
    plt.clf()
    plt.close("all")
