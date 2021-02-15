
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
from . import mp_errors
from . import add_version_plot
from . import prep_plots
from . prep_plots import figsize_x, figsize_y, grid_x, grid_y


def main(
    npd, cld_c, pd, em_float, err_lst, cl_region_c, cl_region_rjct_c,
    stars_out_c, stars_out_rjct_c, N_st_err_rjct, col_0_comb, mag_0_comb,
        **kwargs):
    """
    Make B1 block plots.
    """
    fig = plt.figure(figsize=(figsize_x, figsize_y))
    gs = gridspec.GridSpec(grid_y, grid_x)
    add_version_plot.main(y_fix=1.)

    # Obtain plotting parameters and data.
    err_bar_all = prep_plots.error_bars(
        cld_c['mags'][0], np.nan, err_lst, 'all')

    # Magnitude vs uncertainties diagrams.
    arglist = [
        [gs, pd['colors'], pd['filters'], pd['id_kinem'], cld_c['mags'],
         em_float, cl_region_c, cl_region_rjct_c, stars_out_c,
         stars_out_rjct_c, N_st_err_rjct, err_bar_all]
    ]
    for n, args in enumerate(arglist):
        mp_errors.plot(n, *args)

    plt.suptitle(
        (r"$N_{{accpt}}={}$ , $N_{{rjct}}={}$ (cluster + "
         "field regions, compl frame)").format(
            len(cl_region_c) + len(stars_out_c),
            len(stars_out_rjct_c) + len(cl_region_rjct_c)),
        x=.26, y=1.005)

    # Generate output file.
    fig.tight_layout()
    plt.savefig(
        join(npd['output_subdir'], str(npd['clust_name']) + '_B1'),
        bbox_inches='tight')
    # Close to release memory.
    plt.clf()
    plt.close("all")
