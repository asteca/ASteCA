# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 21:02:33 2014

@author: gabriel
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.ndimage.filters import gaussian_filter


def GA_plot(i, mm_m, mm_a, mm_e, mm_d, params_ga, lkl, lkl_old, ext_imm_indx, isoch_done, generation):
    '''
    Plots several diagrams related to the Genetic Algorithm.
    '''
    
    n_pop, n_gen, fdif, p_cross, cr_sel, p_mut, n_el, n_ei, n_es = params_ga    
    
    fig = plt.figure(figsize=(12, 6))
    gs = gridspec.GridSpec(2, 4)
                              
    ax0 = plt.subplot(gs[0:1, 0:1])
    plt.xlim(-5, n_gen+int(0.05*n_gen))
    plt.ylim(0, 4000)
    ax0.tick_params(axis='both', which='major', labelsize=7)
    plt.xlabel('Generation', fontsize=10)
    plt.ylabel('Likelihood', fontsize=10)
    text1 = '$N = %d\,;\,L_{min}=%0.2f$' '\n' % (i, lkl[0])
    text2 = '$n_{gen}=%d\,;\,n_{pop}=%d$' '\n' % (n_gen, n_pop)
    text3 = '$f_{dif}=%0.2f\,;\,cr_{sel}=%s$' '\n' % (fdif, cr_sel)
    text4 = '$p_{cross}=%0.2f\,;\,p_{mut}=%0.2f$' '\n' % (p_cross, p_mut)
    text5 = '$n_{el}=%d\,;\,n_{ei}=%d\,;\,n_{es}=%d$' % (n_el, n_ei, n_es)
    text = text1+text2+text3+text4+text5
    plt.text(0.41, 0.7, text, transform = ax0.transAxes, \
    bbox=dict(facecolor='white', alpha=0.5), fontsize=8)
    ax0.plot(range(i+1), lkl_old[0], lw=1., c='black')
    ax0.plot(range(i+1), lkl_old[1], lw=1., c='blue')
    for lin in ext_imm_indx:
        plt.axvline(x=lin, linestyle='--', lw=0.6, color='black')


    ax1 = plt.subplot(gs[0:1, 1:2])
    plt.xlim(0, len(isoch_done[1]))
    plt.ylim(lkl[0]-(lkl[0]*0.2), 4000)
    ax1.tick_params(axis='both', which='major', labelsize=7)
    plt.rc('font', size=7)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0), labelsize=7)
    plt.xlabel('Solutions', fontsize=10)
    plt.ylabel('Likelihoods', fontsize=10)
    # Plot points in generation.
    ax1.scatter(range(len(isoch_done[1])), isoch_done[1], s=5, lw=0.2, c='blue')
    plt.axhline(y=lkl[0], linestyle='--', lw=0.8, color='red')
    

    ax2 = plt.subplot(gs[0:1, 2:3])
    plt.xlim(mm_m[0], mm_m[1])
    plt.ylim(mm_a[0], mm_a[1])
    ax2.tick_params(axis='both', which='major', labelsize=7)
    plt.xlabel('$z$', fontsize=12)
    plt.ylabel('$log(age)$', fontsize=12)
    # Plot points in generation.
    ax2.scatter(zip(*generation)[0], zip(*generation)[1], s=7, lw=0.2, c='blue')
    ax2.scatter(generation[0][0], generation[0][1], s=20, c='red')
    
    
    ax3 = plt.subplot(gs[0:1, 3:4])
    plt.xlim(mm_e[0], mm_e[1])
    plt.ylim(mm_d[0], mm_d[1])
    ax3.tick_params(axis='both', which='major', labelsize=7)
    plt.xlabel('$E_{(B-V)}$', fontsize=12)
    plt.ylabel('$dist\;mod$', fontsize=12)
    # Plot points in generation.
    ax3.scatter(zip(*generation)[2], zip(*generation)[3], s=8, lw=0.2, c='blue')
    ax3.scatter(generation[0][2], generation[0][3], s=20, c='red')
    
    
    ax4 = plt.subplot(gs[1:2, 0:1])
    plt.ylim(lkl[0], 3000)
    plt.xlim(mm_m[0], mm_m[1])
    ax4.tick_params(axis='both', which='major', labelsize=7)
    plt.ylabel('Likelihood', fontsize=10)
    plt.xlabel('$z$', fontsize=12)
    text = '$z = %0.4f$' % generation[0][0]
    plt.text(0.67, 0.9, text, transform = ax4.transAxes, \
    bbox=dict(facecolor='white', alpha=0.5), fontsize=8)
    hist, xedges, yedges = np.histogram2d(zip(*isoch_done[0])[0], isoch_done[1], bins=100)
    # H_g is the 2D histogram with a gaussian filter applied
    h_g = gaussian_filter(hist, 2, mode='constant')        
    plt.imshow(h_g.transpose(), origin='lower',
               extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],\
               cmap=plt.get_cmap('gist_yarg'), aspect='auto')
    plt.axvline(x=generation[0][0], linestyle='--', color='red')

    
    ax5 = plt.subplot(gs[1:2, 1:2])
    plt.ylim(lkl[0], 3000)
    plt.xlim(mm_a[0], mm_a[1])
    ax5.tick_params(axis='both', which='major', labelsize=7)
    plt.ylabel('Likelihood', fontsize=10)
    plt.xlabel('$log(age)$', fontsize=12)
    text = '$log(age) = %0.2f$' % generation[0][1]
    plt.text(0.6, 0.9, text, transform = ax5.transAxes, \
    bbox=dict(facecolor='white', alpha=0.5), fontsize=8)
    hist, xedges, yedges = np.histogram2d(zip(*isoch_done[0])[1], isoch_done[1], bins=100)
    # H_g is the 2D histogram with a gaussian filter applied
    h_g = gaussian_filter(hist, 2, mode='constant')        
    plt.imshow(h_g.transpose(), origin='lower',
               extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],\
               cmap=plt.get_cmap('gist_yarg'), aspect='auto')
    plt.axvline(x=generation[0][1], linestyle='--', color='red')


    ax6 = plt.subplot(gs[1:2, 2:3])
    plt.ylim(lkl[0], 3000)
    plt.xlim(mm_e[0], mm_e[1])
    ax6.tick_params(axis='both', which='major', labelsize=7)
    plt.ylabel('Likelihood', fontsize=10)
    plt.xlabel('$E_{(B-V)}$', fontsize=12)
    text = '$E_{(B-V)} = %0.2f$' % generation[0][2]
    plt.text(0.65, 0.9, text, transform = ax6.transAxes, \
    bbox=dict(facecolor='white', alpha=0.5), fontsize=8)
    hist, xedges, yedges = np.histogram2d(zip(*isoch_done[0])[2], isoch_done[1], bins=100)
    # H_g is the 2D histogram with a gaussian filter applied
    h_g = gaussian_filter(hist, 2, mode='constant')        
    plt.imshow(h_g.transpose(), origin='lower',
               extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],\
               cmap=plt.get_cmap('gist_yarg'), aspect='auto')
    plt.axvline(x=generation[0][2], linestyle='--', color='red')
    
    
    ax7 = plt.subplot(gs[1:2, 3:4])
    plt.ylim(lkl[0], 3000)
    plt.xlim(mm_d[0], mm_d[1])
    ax7.tick_params(axis='both', which='major', labelsize=7)
    plt.ylabel('Likelihood', fontsize=10)
    plt.xlabel('$dist\;mod$', fontsize=12)
    text = '$dist\;mod = %0.2f$' % generation[0][3]
    plt.text(0.5, 0.9, text, transform = ax7.transAxes, \
    bbox=dict(facecolor='white', alpha=0.5), fontsize=8)
    hist, xedges, yedges = np.histogram2d(zip(*isoch_done[0])[3], isoch_done[1], bins=100)
    # H_g is the 2D histogram with a gaussian filter applied
    h_g = gaussian_filter(hist, 2, mode='constant')        
    plt.imshow(h_g.transpose(), origin='lower',
               extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],\
               cmap=plt.get_cmap('gist_yarg'), aspect='auto')
    plt.axvline(x=generation[0][3], linestyle='--', color='red')
    
    
    fig.tight_layout()
    # Generate output file for each data file.
    plt.savefig('/home/gabriel/Descargas/GA/GA_'+str(i).zfill(3)+'.png', dpi=150)
    # Close to release memory.
    plt.clf()
    plt.close()
    