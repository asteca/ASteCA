# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 21:26:28 2014

@author: gabriel
"""

import matplotlib.pyplot as plt


def synth_clust_plot(mass_dist, isoch_inter, params, isoch_moved, isoch_cut,
                     isoch_m_d0, isoch_m_d, clust_compl, clust_error):
    '''
    Plot several diagrams related with the synthetic clusters.
    '''

    m, a, e, d = params    
    
    f, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8)) = plt.subplots(2, 4)
    for ax in [ax2, ax3, ax5, ax6, ax7, ax8]:
        ax.set_xlabel('$(C-T_1)$', fontsize=15)
        ax.set_ylabel('$T_1$', fontsize=15)
        
        
    ax1.set_title('Isochrone (interp)')
    ax1.invert_yaxis()
    ax1.set_xlabel('$(C-T_1)_o$', fontsize=15)
    ax1.set_ylabel('$M_{T_1}$', fontsize=15)
    text1 = '$z = %0.4f$' '\n' % m
    text2 = '$log(age) = %0.2f$' % a
    text=text1+text2
    plt.text(0.1, 0.1, text, transform=ax1.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=15)
    plt.text(0.1, 0.9, 'a)', transform=ax1.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
    ax1.scatter(isoch_inter[0], isoch_inter[1], s=10, c='aqua')


    ax2.set_title('Shifted isochrone')
    ax2.invert_yaxis()
    text1 = '$E_{(B-V)} = %0.2f$' '\n' % e
    text2 = '$(m-M)_o = %0.2f$' % d
    text=text1+text2
    plt.text(0.1, 0.1, text, transform=ax2.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=15)
    plt.text(0.1, 0.9, 'b)', transform=ax2.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
    ax2.scatter(isoch_moved[0], isoch_moved[1], s=10, c='azure')  



    ax3.set_title('Max magnitude cut')
    ax3.invert_yaxis()
    plt.text(0.1, 0.9, 'c)', transform=ax3.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
    ax3.scatter(isoch_cut[0], isoch_cut[1], s=10, c='bisque')    


    ax4.set_title('Mass distribution')
    ax4.set_xlabel('$M_{\odot}$', fontsize=15)
    ax4.set_ylabel('$N$', fontsize=15)
    plt.text(0.1, 0.9, 'd)', transform=ax4.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
    plt.text(0.7, 0.9, 'N=%d' % len(mass_dist), transform=ax4.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=14)
    ax4.hist(mass_dist, bins=50, normed=True)
    
    
    ax5.set_title('IMF masses')
    ax5.invert_yaxis()
    plt.text(0.1, 0.9, 'e)', transform=ax5.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
    ax5.scatter(isoch_m_d0[0], isoch_m_d0[1], s=10, c='teal')    


    ax6.set_title('Binarity')
    ax6.invert_yaxis()
    plt.text(0.1, 0.9, 'f)', transform=ax6.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
    ax6.scatter(isoch_m_d[0], isoch_m_d[1], s=10, c='steelblue')


    if isoch_m_d.any():
        ax7.set_title('Completeness')
        ax7.invert_yaxis()
        plt.text(0.1, 0.9, 'g)', transform=ax7.transAxes,
                 bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
        ax7.scatter(clust_compl[0], clust_compl[1], s=10, c='blue')


        ax8.set_title('Add errors')
        ax8.invert_yaxis()
        plt.text(0.1, 0.9, 'h)', transform=ax8.transAxes,
                 bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
        ax8.scatter(clust_error[0], clust_error[2], s=10, c='black')


    plt.show()
    print 'Synthetic cluster plotted'
    raw_input()
