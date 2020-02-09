
import numpy as np
from os.path import join
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from . import add_version_plot
from . import prep_plots
from ..structure.king_profile import centDens, KingProf


def main(npd, pd, max_mag_syn, synth_gen_pars, **kwargs):
    """
    In place for #239
    Plot synthetic clusters.
    """

    if 'D0' in pd['flag_make_plot'] and pd['best_fit_algor'] == 'synth_gen':

        fig = plt.figure(figsize=(30, 25))
        gs = gridspec.GridSpec(10, 12)
        add_version_plot.main(y_fix=1.005)

        extra_pars, model, isoch_moved, isoch_binar, isoch_compl,\
            synth_clust, sigma, cl_dists, x_cl, y_cl, x_fl, y_fl, synth_field,\
            sigma_field, field_dens, CI, rc, rt, xmax, ymax, cx, cy =\
            synth_gen_pars

        coord, x_name, y_name = prep_plots.coord_syst(pd['coords'])
        x_ax0, y_ax = prep_plots.ax_names(
            pd['colors'][0], pd['filters'][0], 'mag')

        # # Plot field + cluster.
        # # figsize(x1, y1), GridSpec(y2, x2) --> To have square plots: x1/x2 =
        # # y1/y2 = 2.5
        # fig = plt.figure(figsize=(15, 10))  # create the top-level container
        # gs = gridspec.GridSpec(4, 6)  # create a GridSpec object

        # Cluster's stars CMD.
        ax = plt.subplot(gs[0:2, 0:2])
        # Set plot limits
        # col_min = min(min(synth_clust[1]), min(synth_field[1])) - 0.2
        # col_max = max(max(synth_clust[1]), max(synth_field[1])) + 0.2
        # mag_min, mag_max = max(max_mag_syn + .5, max(isoch_moved[0]) + .5), \
        #     min(isoch_moved[0]) - 0.2
        # plt.xlim(col_min, col_max)
        # plt.ylim(mag_min, mag_max)
        # # Set axis labels
        plt.xlabel('$' + x_ax0 + '$', fontsize=12)
        plt.ylabel('$' + y_ax + '$', fontsize=12)
        # Set minor ticks
        ax.minorticks_on()
        # Set grid
        ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5)
        t1 = r'$N_{{cluster}} = ${}'.format(len(isoch_moved[1]))
        t2 = r'${}_{{min}} = ${:.1f}'.format(y_ax, max_mag_syn)
        t3 = r'$z = ${:.4f}'.format(model[0])
        t4 = r'$\log(age) = ${:.1f}'.format(model[1])
        t5 = r'$A_{{V}} = ${:.1f} mag'.format(model[2] * 3.1)
        t6 = r'$\mu = ${:.2f} mag'.format(model[3])
        t7 = r'$M = ${:.0f} '.format(model[4]) + r'$M_{\odot}$'
        t8 = r'$b_{{fr}} = ${:.2f} '.format(model[5])
        t9 = r'$\sum m_i = ${:.0f}'.format(
            extra_pars[2].sum()) + r'$M_{\odot}$'
        text = t1 + '\n' + t2 + '\n' + t3 + '\n' + t4 + '\n' + t5 + '\n' +\
            t6 + '\n' + t7 + '\n' + t8 + '\n' + t9
        # Create an empty plot with the required text.
        plt.plot([], label=text)
        # Remove the handle from the legend box.
        plt.legend(handlelength=0, fontsize=11)
        # Plot stars.
        plt.scatter(
            isoch_moved[1], isoch_moved[0], marker='o', c='r', s=15., lw=0.3)
        plt.axhline(y=max_mag_syn, linestyle='-', color='green', lw=1.5)
        ax.invert_yaxis()

        # Completeness
        ax = plt.subplot(gs[0:2, 2:4])
        # Set plot limits
        # Set axis labels
        plt.xlabel(r'${}$'.format(y_ax), fontsize=12)
        plt.ylabel(r'$N$', fontsize=12)
        # Set minor ticks
        ax.minorticks_on()
        # Backg color.
        ax.set_facecolor('#D8D8D8')
        # Set grid
        ax.grid(b=True, which='major', color='w', linestyle='-', lw=.5)
        # Plot stars.
        _, edges, _ = plt.hist(
            isoch_binar[0], bins=50, color='blue', histtype='step', zorder=4,
            label='Before removal (N={})'.format(len(isoch_binar[0])))
        plt.hist(
            isoch_compl[0], bins=edges, color='red', histtype='step',
            label='After removal    (N={})'.format(len(isoch_compl[0])),
            ls='dashed', hatch="/", zorder=4)
        # Legends.
        leg = plt.legend(
            fancybox=True, loc='upper left', numpoints=1, fontsize=11)
        # Set the alpha value of the legend.
        leg.get_frame().set_alpha(0.7)

        # Magnitude error
        ax = plt.subplot(gs[0, 4:6])
        ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5)
        # Set plot limits
        plt.xlim(min(synth_clust[0]) - 0.5, max(synth_clust[0]) + 0.5)
        plt.ylim(-0.001, max(sigma[0]) + .25 * max(sigma[0]))
        # Set axis labels
        plt.ylabel(r'$\sigma_{}$'.format(y_ax), fontsize=12)
        plt.xlabel(r'${}$'.format(y_ax), fontsize=12)
        # Set minor ticks
        ax.minorticks_on()
        # plt.text(0.25, 0.85, '$N_{cluster} = %d$' % len(synth_clust[0]),
        #          transform=ax.transAxes,
        #          bbox=dict(facecolor='white', alpha=0.75), fontsize=13)
        # Plot stars errors.
        plt.scatter(synth_field[0], sigma_field[0], marker='o', c='k', s=5)
        plt.scatter(
            synth_clust[0], sigma[0], marker='o', c='r', s=5, zorder=3)
        plt.axvline(x=max_mag_syn, linestyle='-', color='green', lw=1.5)

        # Color color error
        ax = plt.subplot(gs[1, 4:6])
        ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5)
        # Set plot limits
        plt.xlim(min(synth_clust[0]) - 0.5, max(synth_clust[0]) + 0.5)
        plt.ylim(-0.001, max(sigma[1]) + .25 * max(sigma[1]))
        # Set axis labels
        plt.ylabel(r'$\sigma_{}$'.format(x_ax0), fontsize=12)
        plt.xlabel(r'${}$'.format(y_ax), fontsize=12)
        # Set minor ticks
        ax.minorticks_on()
        # Plot stars errors.
        plt.scatter(synth_field[0], sigma_field[1], marker='o', c='k', s=5)
        plt.scatter(synth_clust[0], sigma[1], marker='o', c='r', s=5, zorder=3)
        plt.axvline(x=max_mag_syn, linestyle='-', color='green', lw=1.5)

        # King's profile
        ax = plt.subplot(gs[2:4, 0:2])
        ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5)
        plt.title(r"$r_{{c}}=${:.0f}, $r_{{t}}=${:.0f}".format(rc, rt))
        r = np.linspace(0., rt, 1000)
        plt.plot(
            r, centDens(synth_clust.shape[1], rc, rt, r) *
            KingProf(r, rc, rt), c='r',
            label=r"$K_{{cp}}={:.2f}$".format(np.log10(rt / rc)), zorder=5)
        vals, edges = np.histogram(cl_dists, 20)
        areas = np.pi * edges**2  # TODO handle areas outside of frame?
        rings_a = areas[1:] - areas[:-1]
        dens = vals / rings_a
        plt.bar(edges[:-1], dens, width=5)
        plt.xlabel("Distance to center", fontsize=12)
        plt.ylabel(r"$\rho$ [st/area] (normalized)", fontsize=12)
        plt.legend(loc='upper right', fontsize=11)

        # x,y finding chart of full frame
        ax = plt.subplot(gs[2:4, 2:4])
        plt.xlim(0., xmax)
        plt.ylim(0., ymax)
        plt.xlabel('x (px)', fontsize=12)
        plt.ylabel('y (px)', fontsize=12)
        ax.set_title(r"$N_{{cl}}=${}, $N_{{fr}}=${}, CI={:.2f}".format(
            len(x_cl), len(x_fl), CI))
        # Set minor ticks
        ax.minorticks_on()
        circle = plt.Circle((cx, cy), rt, color='r', fill=False, lw=2.)
        fig.gca().add_artist(circle)
        mags = synth_clust[0].tolist() + synth_field[0].tolist()
        sizes = prep_plots.star_size(mags)
        plt.scatter(
            x_fl, y_fl, marker='o', c='black', s=sizes[synth_clust[0].size:])
        plt.scatter(
            x_cl, y_cl, marker='o', c='r', edgecolor='r',
            s=sizes[:synth_clust[0].size])

        # Full region CMD.
        ax = plt.subplot(gs[2:4, 4:6])
        plt.xlabel('$' + x_ax0 + '$', fontsize=12)
        plt.ylabel('$' + y_ax + '$', fontsize=12)
        # Set minor ticks
        ax.minorticks_on()
        # Set grid
        ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5,
                zorder=1)
        # Plot stars.
        plt.scatter(
            synth_field[1], synth_field[0], marker='o', c='k', s=5, lw=.3,
            edgecolor='w', zorder=1)
        # Identify binaries
        b_msk = extra_pars[0] == 2.
        plt.scatter(
            synth_clust[1][b_msk], synth_clust[0][b_msk], marker='o',
            c='b', s=10, lw=.2, edgecolor='k', zorder=4, label="Binary")
        plt.scatter(
            synth_clust[1][~b_msk], synth_clust[0][~b_msk], marker='o',
            c='r', s=10, lw=.2, edgecolor='k', zorder=2, label="Single")
        ax.invert_yaxis()
        plt.legend(numpoints=1, fontsize=11)

        fig.tight_layout()
        # Generate output file.
        plt.savefig(
            join(npd['synth_gen_fold'], str(npd['clust_name']) +
                 '_D0.' + pd['plot_frmt']), dpi=pd['plot_dpi'],
            bbox_inches='tight')
        # Close to release memory.
        plt.clf()
        plt.close("all")

        print("<<Plots for D0 block created>>")
    else:
        print("<<Skip D0 block plot>>")
