
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from os.path import join
from . import mp_AD_test
from . import add_version_plot
from . prep_plots import figsize_x, figsize_y, grid_x, grid_y


def main(
    npd, pd, flag_ad_test, ad_cl, ad_fr, ad_cl_fr_p, ad_cl_fr_pk,
        **kwargs):
    """
    Make B3 block plots.
    """
    if 'B3' in pd['flag_make_plot'] and flag_ad_test:
        fig = plt.figure(figsize=(figsize_x, figsize_y))
        gs = gridspec.GridSpec(grid_y, grid_x)
        add_version_plot.main(y_fix=.999)

        arglist = [
            # pl_ad_test
            [gs, pd['plot_style'], flag_ad_test, ad_cl, ad_fr, pd['id_kinem']],
            # pl_ad_pvals_phot
            [gs, pd['plot_style'], flag_ad_test, ad_cl_fr_p],
            # pl_ad_pvals_pk
            [gs, pd['plot_style'], flag_ad_test, ad_cl_fr_pk, pd['id_kinem']]
        ]
        for n, args in enumerate(arglist):
            mp_AD_test.plot(n, *args)

        # Generate output file.
        fig.tight_layout()
        plt.savefig(
            join(npd['output_subdir'], str(npd['clust_name']) + '_B3'),
            bbox_inches='tight')
        # Close to release memory.
        plt.clf()
        plt.close("all")

        print("<<Plots for B3 block created>>")
    else:
        print("<<Skip B3 block plot>>")
