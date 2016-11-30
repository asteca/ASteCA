
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def synth_clust_plot(
        N_fc, mass_dist, isochrone, synth_cl_params, isoch_moved, isoch_cut,
        isoch_mass, isoch_binar, binar_idx, isoch_compl, compl_idx,
        synth_clust, path):
    '''
    Plot several diagrams related with the synthetic clusters.
    '''

    m, a, e, d, M_total, bin_frac = synth_cl_params
    # Index of m_ini, stored in the theoretical isochrones.
    m_ini = N_fc[0] + N_fc[1] + 2 * N_fc[1]

    # figsize(x1, y1), GridSpec(y2, x2) --> To have square plots: x1/x2 =
    # y1/y2 = 2.5
    fig = plt.figure(figsize=(10, 20))  # create the top-level container
    gs = gridspec.GridSpec(8, 4)  # create a GridSpec object

    ax1 = plt.subplot(gs[0:2, 0:2])
    ax1.set_title('Isochrone (interpolated)')
    ax1.invert_yaxis()
    ax1.set_xlabel('$color_o$', fontsize=15)
    ax1.set_ylabel('$M_{o}$', fontsize=15)
    ax1.minorticks_on()
    ax1.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    text1 = '$z = {}$\n'.format(m)
    text2 = '$log(age/yr) = %0.2f$' % a
    text = text1 + text2
    plt.text(0.1, 0.1, text, transform=ax1.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=15)
    plt.text(0.05, 0.92, 'a', transform=ax1.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
    plt.text(0.75, 0.92, 'N=%d' % len(isochrone[0]), transform=ax1.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=14)
    ax1.scatter(isochrone[N_fc[0]], isochrone[0], s=30, c='steelblue', lw=0.1)

    ax2 = plt.subplot(gs[0:2, 2:4])
    ax2.set_title('Shifted')
    ax2.invert_yaxis()
    ax2.minorticks_on()
    ax2.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    text1 = '$E_{(B-V)} = %0.2f$' '\n' % e
    text2 = '$(m-M)_o = %0.2f$' % d
    text = text1 + text2
    plt.text(0.1, 0.1, text, transform=ax2.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=15)
    plt.text(0.05, 0.92, 'b', transform=ax2.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
    plt.text(0.75, 0.92, 'N=%d' % len(isoch_moved[0]), transform=ax2.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=14)
    ax2.scatter(isoch_moved[N_fc[0]], isoch_moved[0], s=30, c='steelblue',
                lw=0.1)

    ax3 = plt.subplot(gs[2:4, 0:2])
    ax3.set_title('Max magnitude cut')
    ax3.invert_yaxis()
    ax3.minorticks_on()
    ax3.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    plt.text(0.05, 0.92, 'c', transform=ax3.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
    plt.text(0.75, 0.92, 'N=%d' % len(isoch_cut[0]), transform=ax3.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=14)
    ax3.scatter(isoch_cut[N_fc[0]], isoch_cut[0], s=30, c='steelblue', lw=0.3)

    ax4 = plt.subplot(gs[2:4, 2:4])
    ax4.set_title('Mass distribution')
    ax4.set_xlabel('$M_{\odot}$', fontsize=15)
    ax4.set_ylabel('$N$', fontsize=15)
    ax4.minorticks_on()
    plt.xlim(0, 10)
    ax4.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    plt.arrow(0.83, 0.085, 0.1, 0., transform=ax4.transAxes, fc="k", ec="k",
              lw=1.5, head_width=0.01)
    # This value should be entered manually to not disrupt the 'synth_cluster'
    # module.
    m_high = 300.
    m_high_str = '(' + str(m_high) + ')'
    plt.text(0.83, 0.033, m_high_str, transform=ax4.transAxes, fontsize=12)
    plt.text(0.05, 0.92, 'd', transform=ax4.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
    text1 = 'N=%d\n' % len(mass_dist)
    text2 = '$M_T = {:.0f}\,M_{{\odot}}$\n'.format(M_total)
    text3 = '$M = {:.0f}\,M_{{\odot}}$'.format(sum(mass_dist))
    text = text1 + text2 + text3
    plt.text(0.56, 0.8, text, transform=ax4.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=14)
    min_d, max_d = min(mass_dist), max(mass_dist)
    bw = 0.4
    ax4.hist(mass_dist, bins=np.arange(min_d, max_d + bw, bw))

    ax5 = plt.subplot(gs[4:6, 0:2])
    ax5.set_title('IMF sampling')
    min_x, max_x = min(isoch_mass[N_fc[0]]), max(isoch_mass[N_fc[0]])
    plt.xlim(min_x - 0.75, max_x + 0.75)
    ax5.invert_yaxis()
    ax5.minorticks_on()
    ax5.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    plt.text(0.05, 0.92, 'e', transform=ax5.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
    text1 = 'N=%d\n' % len(isoch_mass[0])
    text2 = '$M = {:.0f}\,M_{{\odot}}$'.format(sum(isoch_mass[m_ini]))
    text = text1 + text2
    plt.text(0.6, 0.87, text, transform=ax5.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=14)
    ax5.scatter(isoch_mass[N_fc[0]], isoch_mass[0], s=30, c='steelblue',
                lw=0.5)

    ax6 = plt.subplot(gs[4:6, 2:4])
    ax6.set_title('Binarity')
    plt.xlim(min_x - 0.75, max_x + 0.75)
    ax6.invert_yaxis()
    ax6.minorticks_on()
    ax6.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    plt.text(0.05, 0.92, 'f', transform=ax6.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
    text1 = 'N=%d\n' % len(isoch_binar[0])
    text2 = '$b_{{frac}} = {}$\n'.format(bin_frac)
    text3 = '$M = {:.0f}\,M_{{\odot}}$'.format(sum(isoch_binar[m_ini]))
    text = text1 + text2 + text3
    plt.text(0.6, 0.81, text, transform=ax6.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=14)
    binar_x, binar_y = isoch_binar[N_fc[0]][binar_idx],\
        isoch_binar[0][binar_idx]
    sing_x, sing_y = np.delete(isoch_binar[N_fc[0]], binar_idx),\
        np.delete(isoch_binar[0], binar_idx)
    ax6.scatter(sing_x, sing_y, s=30, c='steelblue', lw=0.5)
    ax6.scatter(binar_x, binar_y, s=30, c='red', lw=0.5)

    if isoch_compl.any():
        ax7 = plt.subplot(gs[6:8, 0:2])
        ax7.set_title('Completeness')
        plt.xlim(min_x - 0.75, max_x + 0.75)
        ax7.invert_yaxis()
        ax7.minorticks_on()
        ax7.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
        plt.text(0.05, 0.92, 'g', transform=ax7.transAxes,
                 bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
        text1 = 'N=%d\n' % len(isoch_compl[0])
        text2 = '$M = {:.0f}\,M_{{\odot}}$'.format(sum(isoch_compl[m_ini]))
        text = text1 + text2
        plt.text(0.6, 0.87, text, transform=ax7.transAxes,
                 bbox=dict(facecolor='white', alpha=0.5), fontsize=14)
        binar_x, binar_y = isoch_compl[N_fc[0]][compl_idx],\
            isoch_compl[0][compl_idx]
        sing_x, sing_y = np.delete(isoch_compl[N_fc[0]], compl_idx),\
            np.delete(isoch_compl[0], compl_idx)
        ax7.scatter(binar_x, binar_y, s=30, c='red', lw=0.5)
        ax7.scatter(sing_x, sing_y, s=30, c='steelblue', lw=0.5)

        ax8 = plt.subplot(gs[6:8, 2:4])
        ax8.set_title('Errors')
        plt.xlim(min_x - 0.75, max_x + 0.75)
        ax8.invert_yaxis()
        ax8.minorticks_on()
        ax8.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
        plt.text(0.05, 0.92, 'h', transform=ax8.transAxes,
                 bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
        plt.text(
            0.6, 0.92, 'N=%d' % len(synth_clust[0][0]),
            transform=ax8.transAxes, bbox=dict(facecolor='white', alpha=0.5),
            fontsize=14)
        # ax8.scatter(synth_clust[0][N_fc[0]], synth_clust[0][0], marker='o',
        #             s=30, c='#4682b4', lw=0.5)
        binar_x, binar_y = synth_clust[0][N_fc[0]][compl_idx],\
            synth_clust[0][0][compl_idx]
        sing_x, sing_y = np.delete(synth_clust[0][N_fc[0]], compl_idx),\
            np.delete(synth_clust[0][0], compl_idx)
        ax8.scatter(sing_x, sing_y, s=30, c='steelblue', lw=0.5)
        ax8.scatter(binar_x, binar_y, s=30, c='red', lw=0.5)

    for ax in [ax2, ax3, ax5, ax6, ax7, ax8]:
        ax.set_xlabel('$color$', fontsize=15)
        ax.set_ylabel('$magnitude$', fontsize=15)

    # plt.show()
    fig.tight_layout()
    plt.savefig(path, dpi=300)
    # Close to release memory.
    plt.clf()
    plt.close()
    print 'Synthetic cluster plotted'
    import pdb; pdb.set_trace()  # breakpoint 0290b401 //

