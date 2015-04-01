# -*- coding: utf-8 -*-
"""
Created on Wed Dic 10 12:00:00 2014

@author: gabriel
"""

#import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.offsetbox as offsetbox


#def pl_likel_dens():
    #'''
    #Produce output image for top tier models.
    #'''
    #ax = plt.subplot(gs[10:12, 0:2])
    ## Axis limits.
    #plt.xlim(m_min, m_max)
    #plt.ylim(a_min, a_max)
    #plt.xlabel('$z$', fontsize=16)
    #plt.ylabel('$log(age)$', fontsize=16)
    #plt.minorticks_on()
    ## Plot best fit point.
    #plt.scatter(m, a, marker='o', c='r', s=30)
    ## Check if errors in both dimensions are defined.
    #if all([i > 0. for i in [e_m, e_a]]):
        ## Plot ellipse error.
        #plt.gca()
        #ellipse = Ellipse(xy=(m, a), width=2 * e_m, height=2 * e_a,
                                #edgecolor='r', fc='None', lw=1.)
        #ax.add_patch(ellipse)
    #elif e_m < 0.:
        #plt.errorbar(m, a, yerr=e_a, color='r')
    #elif e_a < 0.:
        #plt.errorbar(m, a, xerr=e_m, color='r')
    ## Plot density map.
    #hist, xedges, yedges = np.histogram2d(zip(*model_done[0])[0],
        #zip(*model_done[0])[1], bins=100)
    ## H_g is the 2D histogram with a gaussian filter applied
    #h_g = gaussian_filter(hist, 2, mode='constant')
    #plt.imshow(h_g.transpose(), origin='lower',
               #extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
               #cmap=plt.get_cmap('Blues'), aspect='auto')
    ## Plot top tiers.
    #top_tiers = [x for y, x in sorted(zip(model_done[1], model_done[0]))]
    #print top_tiers[:10]
    #top_t = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    #for i, txt in enumerate(top_t):
        #ax.scatter(zip(*top_tiers)[0][i], zip(*top_tiers)[1][i], s=1)
        #ax.annotate(txt, (zip(*top_tiers)[0][i], zip(*top_tiers)[1][i]),
            #size=9)

def pl_bf_synth_cl(N, gs, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax,
    y_ax, synth_clst, cp_r, shift_isoch):
    '''
    Top tiers synthetic clusters obtained.
    '''
    gs_map = {1: gs[0:2, 0:2], 2: gs[0:2, 2:4], 3: gs[0:2, 4:6],
        4: gs[0:2, 6:8], 5: gs[0:2, 8:10], 6: gs[2:4, 0:2], 7: gs[2:4, 2:4],
        8: gs[2:4, 4:6], 9: gs[2:4, 6:8], 10: gs[2:4, 8:10]}
    ax = plt.subplot(gs_map.get(N))
    #Set plot limits
    plt.xlim(x_min_cmd, x_max_cmd)
    plt.ylim(y_min_cmd, y_max_cmd)
    #Set axis labels
    plt.xlabel('$' + x_ax + '$', fontsize=18)
    plt.ylabel('$' + y_ax + '$', fontsize=18)
    # Set minor ticks
    ax.minorticks_on()
    ax.xaxis.set_major_locator(MultipleLocator(1.0))
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    # Add text box
    text1 = '\n$N_{{synth}} = {}$\n'.format(len(synth_clst[0]))
    text2 = '$z = {}$\n'.format(cp_r[0])
    text3 = '$log(age) = {}$\n'.format(cp_r[1])
    text4 = '$E_{{(B-V)}} = {}$\n'.format(cp_r[2])
    text5 = '$(m-M)_o = {}$\n'.format(cp_r[3])
    text6 = '$M_{{\odot}} = {}$\n'.format(cp_r[4])
    text7 = '$b_{{frac}} = {}$'.format(cp_r[5])
    text = text1 + text2 + text3 + text4 + text5 + text6 + text7
    ob = offsetbox.AnchoredText(text, pad=.2, loc=1, prop=dict(size=12))
    ob.patch.set(alpha=0.6)
    ax.add_artist(ob)
    # Number of top tier.
    ob = offsetbox.AnchoredText('{}'.format(N), pad=0.2, loc=2,
        prop=dict(size=12))
    ob.patch.set(alpha=0.85)
    ax.add_artist(ob)
    # Plot isochrone.
    plt.plot(shift_isoch[0], shift_isoch[1], 'r', lw=1.2)
    # Plot synth clust.
    plt.scatter(synth_clst[0], synth_clst[2], marker='o', s=40,
                c='#4682b4', lw=0.5)


def plot(N, *args):
    '''
    Handle each plot separately.
    '''
    plt_map = dict.fromkeys([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [pl_bf_synth_cl,
        'synthetic cluster top tier number'])

    fxn = plt_map.get(N, None)[0]
    if fxn is None:
        raise ValueError("  ERROR: there is no plot {}.".format(N))

    try:
        fxn(N, *args)
    except:
        print("  WARNING: error when plotting {} {}.".format(plt_map.get(N)[1],
            N))
        if not args[7].any():
            print ("           (synthetic cluster is empty).")
