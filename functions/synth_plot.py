# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 21:26:28 2014

@author: gabriel
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def synth_clust_plot(mass_dist, isochrone, model, isoch_moved, isoch_cut,
    isoch_mass, isoch_binar, isoch_compl, isoch_error, path):
    '''
    Plot several diagrams related with the synthetic clusters.
    '''

    m, a, e, d, mass, bin_f = model

    # figsize(x1, y1), GridSpec(y2, x2) --> To have square plots: x1/x2 =
    # y1/y2 = 2.5
    fig = plt.figure(figsize=(10, 20))  # create the top-level container
    gs = gridspec.GridSpec(8, 4)  # create a GridSpec object

    ax1 = plt.subplot(gs[0:2, 0:2])
    ax1.set_title('Isochrone (interpolated)')
    ax1.invert_yaxis()
    ax1.set_xlabel('$(B-V)_o$', fontsize=15)
    ax1.set_ylabel('$M_{V}$', fontsize=15)
    ax1.minorticks_on()
    ax1.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    text1 = '$z = %0.4f$' '\n' % m
    text2 = '$log(age) = %0.2f$' % a
    text = text1 + text2
    plt.text(0.1, 0.1, text, transform=ax1.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=15)
    plt.text(0.05, 0.92, 'a', transform=ax1.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
    plt.text(0.75, 0.92, 'N=%d' % len(isochrone[0]), transform=ax1.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=14)
    ax1.scatter(isochrone[0], isochrone[1], s=30, c='steelblue', lw=0.1)

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
    ax2.scatter(isoch_moved[0], isoch_moved[1], s=30, c='steelblue', lw=0.1)

    ax3 = plt.subplot(gs[2:4, 0:2])
    ax3.set_title('Max magnitude cut')
    ax3.invert_yaxis()
    ax3.minorticks_on()
    ax3.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    plt.text(0.05, 0.92, 'c', transform=ax3.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
    plt.text(0.75, 0.92, 'N=%d' % len(isoch_cut[0]), transform=ax3.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=14)
    ax3.scatter(isoch_cut[0], isoch_cut[1], s=30, c='steelblue', lw=0.3)

    ax4 = plt.subplot(gs[2:4, 2:4])
    ax4.set_title('Mass distribution')
    ax4.set_xlabel('$M_{\odot}$', fontsize=15)
    ax4.set_ylabel('$N$', fontsize=15)
    ax4.minorticks_on()
    plt.xlim(0, 10)
    ax4.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    plt.arrow(0.83, 0.085, 0.1, 0., transform=ax4.transAxes, fc="k", ec="k",
        lw=1.5, head_width=0.01)
    plt.text(0.83, 0.033, '(100)', transform=ax4.transAxes, fontsize=12)
    plt.text(0.05, 0.92, 'd', transform=ax4.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
    text1 = 'N=%d\n' % len(mass_dist)
    text2 = '$M_T = {}\,M_{{\odot}}$\n'.format(mass)
    text3 = '$M = {:.1f}\,M_{{\odot}}$'.format(sum(mass_dist))
    text = text1 + text2 + text3
    plt.text(0.56, 0.8, text, transform=ax4.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=14)
    min_d, max_d = min(mass_dist), max(mass_dist)
    bw = 0.4
    ax4.hist(mass_dist, bins=np.arange(min_d, max_d + bw, bw))

    ax5 = plt.subplot(gs[4:6, 0:2])
    ax5.set_title('IMF masses')
    min_x, max_x = min(isoch_mass[0]), max(isoch_mass[0])
    plt.xlim(min_x - 0.75, max_x + 0.75)
    ax5.invert_yaxis()
    ax5.minorticks_on()
    ax5.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    plt.text(0.05, 0.92, 'e', transform=ax5.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
    text1 = 'N=%d\n' % len(isoch_mass[0])
    text2 = '$M = {:.1f}\,M_{{\odot}}$'.format(sum(isoch_mass[2]))
    text = text1 + text2
    plt.text(0.6, 0.87, text, transform=ax5.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=14)
    ax5.scatter(isoch_mass[0], isoch_mass[1], s=30, c='steelblue', lw=0.5)

    ax6 = plt.subplot(gs[4:6, 2:4])
    ax6.set_title('Binarity')
    plt.xlim(min_x - 0.75, max_x + 0.75)
    ax6.invert_yaxis()
    ax6.minorticks_on()
    ax6.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    plt.text(0.05, 0.92, 'f', transform=ax6.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
    text1 = 'N=%d\n' % len(isoch_binar[0])
    text2 = '$b_{{frac}} = {}$\n'.format(bin_f)
    text3 = '$M = {:.1f}\,M_{{\odot}}$'.format(sum(isoch_binar[2]))
    text = text1 + text2 + text3
    plt.text(0.6, 0.81, text, transform=ax6.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=14)
    ax6.scatter(isoch_binar[0], isoch_binar[1], s=30, c='steelblue', lw=0.5)

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
        text2 = '$M = {:.1f}\,M_{{\odot}}$'.format(sum(isoch_compl[2]))
        text = text1 + text2
        plt.text(0.6, 0.87, text, transform=ax7.transAxes,
            bbox=dict(facecolor='white', alpha=0.5), fontsize=14)
        ax7.scatter(isoch_compl[0], isoch_compl[1], s=30, c='steelblue', lw=0.5)

        ax8 = plt.subplot(gs[6:8, 2:4])
        ax8.set_title('Errors')
        plt.xlim(min_x - 0.75, max_x + 0.75)
        ax8.invert_yaxis()
        ax8.minorticks_on()
        ax8.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
        plt.text(0.05, 0.92, 'h', transform=ax8.transAxes,
                 bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
        plt.text(0.75, 0.92, 'N=%d' % len(isoch_error[0]),
            transform=ax8.transAxes, bbox=dict(facecolor='white', alpha=0.5),
            fontsize=14)
        ax8.scatter(isoch_error[0], isoch_error[2], marker='o', s=30,
            c='#4682b4', lw=0.5)

    for ax in [ax2, ax3, ax5, ax6, ax7, ax8]:
        ax.set_xlabel('$(B-V)$', fontsize=15)
        ax.set_ylabel('$V$', fontsize=15)

    plt.show()
    #fig.tight_layout()
    #plt.savefig(path, dpi=300)
    print 'Synthetic cluster plotted'
    raw_input()
