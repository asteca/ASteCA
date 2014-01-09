# -*- coding: utf-8 -*-
"""
Created on Sat Sep 21 17:10:04 2013

@author: gabriel
"""

import scipy.spatial.distance as sp

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator


def track_distance(t_ind):
    '''
    Function that takes as input an index pointing to a track/isochrone
    and returns the 'score' for that track, that is: the number obtained
    after summing over all the minimum distances from each cluster star
    to the point in the track closest to it. The smaller this value is,
    the better the track fits the cluster.
    '''
    
    # Define track according to index passed.
    track = funcs[t_ind]

    # Interpolate the track so all tracks will have equidistant
    # points. This is VERY important, otherwise the block that
    # finds the point located at the minimum distance with the
    # cluster stars will not work correctly.
    t = np.linspace(0, 1, len(track[0]))
    # Generate twice as many interpolating points as cluster's stars present.
    # Use a value of 500 if the number of stars is smaller than that and
    # 2000 if it's larger.
    num_inter = min(2000, max(500, 2*len(data[0][0])))
    t2 = np.linspace(0, 1, num_inter)
    
    # One-dimensional linear interpolation.
    x2 = np.interp(t2, t, track[0])
    y2 = np.interp(t2, t, track[1])
    # Store track interpolated values.
    track_inter = [x2,y2]

    # Get distances of *every* star to *every* interpolated point in this
    # isochrone/function. The list is made of N sub-lists where N is the
    # number of cluster stars and M items in each sub-list where M is
    # the number of points in the function/track/isochrone.
    # dist_m = [[star_1], [star_2], ..., [star_N]]
    # [star_i] = [dist_1, dist_2, ..., dist_M]
    dist_m = np.array(sp.cdist(zip(*data[0]), zip(*track_inter), 'euclidean'))
                
    # Identify closest point in track for each star, add weight to this
    # minimum distance and save the weighted distance.
    # The parameter 'mag_weight' is another weight added so brighter stars 
    # will have a bigger impact in the fitting process than dimmer ones.
    mag_weight = 1./np.array(data[0][1])
    func_dist = dist_m.min(axis=1)*np.array(weights[0])*mag_weight
    
    # Sum all the weighted minimum distances to every star for this
    # function/isochrone and append value to 'score' list. This final value
    # is a measure of how good the isochrone fits the data (group of stars),
    # the smaller the value, the better the fit.
    func_score = sum(func_dist)
    print 'orig', func_score
    
    return func_score
        
        
        
def track_distance_err_weight(t_ind):
    '''
    Function that takes as input an index pointing to a track/isochrone
    and returns the 'score' for that track, that is: the number obtained
    after summing over all the *error weighted* distances from each
    cluster star to the point in the track closest to it. The smaller this
    value is, the better the track fits the cluster.
    '''
    
    # Define track according to index passed.
    track = funcs[t_ind]

    # Interpolate the track so all tracks will have equidistant
    # points. This is VERY important, otherwise the block that
    # finds the point located at the minimum dist with the
    # cluster stars will not work correctly.
    t = np.linspace(0, 1, len(track[0]))
    # Generate twice as many interpolating points as cluster's stars present.
    # Use a value of 500 if the number of stars is smaller than that and
    # 2000 if it's larger.
    num_inter = min(2000, max(500, 2*len(data[0][0])))
    t2 = np.linspace(0, 1, num_inter)
    
    # One-dimensional linear interpolation.
    x2 = np.interp(t2, t, track[0])
    y2 = np.interp(t2, t, track[1])
    # Store interpolated values.
    track_inter = [x2,y2]

    # Get distances of *every* star to *every* point in this
    # isochrone/function. The list is made of N sub-lists where N is the
    # number of cluster stars and M items in each sub-list where M is
    # the number of points in the function/track/isochrone.
    # dist_m = [[star_1], [star_2], ..., [star_N]]
    # [star_i] = [dist_1, dist_2, ..., dist_M]

    dist_m = []
    # Iterate through every star in cluster.
    for indx, x_coo in enumerate(data[0][0]):
        y_coo = data[0][1][indx]
        dist_star = []
        # Get weighted distances for this point to every point in track.
        # Iterate through every point in track.
        for indx2, x_coo2 in enumerate(track_inter[0]):
            y_coo2 = track_inter[1][indx2]
            # Weighted distance in x.
            err_x = errors[0][0][indx] if errors[0][0][indx] != 0. else 0.0001
            err_x = 1.
            x_dist_weight = (x_coo-x_coo2)/err_x
            # Weighted distance in y.
            err_y = errors[0][1][indx] if errors[0][1][indx] != 0. else 0.0001
            err_y = 1.
            y_dist_weight = (y_coo-y_coo2)/err_y 
            # Weighted distance between star in cluster passed and this point
            # from this track.
            dist = np.sqrt(x_dist_weight**2 + y_dist_weight**2)
            # Append weighted distance value to list.
            dist_star.append(round(dist, 8))
        dist_m.append(dist_star)

    # Identify closest point in track for each star, add weight to this
    # minimum distance and save the weighted distance.
    func_dist = np.array(dist_m).min(axis=1)*np.array(weights[0])
    
    # Sum all the weighted minimum distances to every star for this
    # function/isochrone and append value to 'score' list. This final value
    # is a measure of how good the isochrone fits the data (group of stars),
    # the smaller the value, the better the fit.
    func_score = sum(func_dist)
    print 'slow', func_score
    
    return func_score        
        
        
        
def fast_wdist(t_ind):
    """
    Compute the weighted euclidean distance between two arrays of points:

    D{i,j} = 
    sqrt( (A{0,i} - B{0,j} / W{0,i})^2 + ... + (A{k,i} - B{k,j} / W{k,i})^2 )

    inputs:
        A is an (k, m) array of coordinates
        B is an (k, n) array of coordinates
        W is an (k, m) array of weights

    returns:
        D is an (m, n) array of weighted euclidean distances
    """
    
    A = data[0]
    
    # Define track according to index passed.
    track = funcs[t_ind]

    # Interpolate the track so all tracks will have equidistant
    # points. This is VERY important, otherwise the block that
    # finds the point located at the minimum dist with the
    # cluster stars will not work correctly.
    t = np.linspace(0, 1, len(track[0]))
    # Generate twice as many interpolating points as cluster's stars present.
    # Use a value of 500 if the number of stars is smaller than that and
    # 2000 if it's larger.
    num_inter = min(2000, max(500, 2*len(data[0][0])))
    t2 = np.linspace(0, 1, num_inter)
    
    # One-dimensional linear interpolation.
    x2 = np.interp(t2, t, track[0])
    y2 = np.interp(t2, t, track[1])
    # Store interpolated values.
    track_inter = np.array([x2,y2])
    
    B = track_inter
    
    W = errors[0]

    # compute the differences and apply the weights in one go using
    # broadcasting jujitsu. the result is (n, k, m)
    wdiff = (A[np.newaxis,...] - B[np.newaxis,...].T) / W[np.newaxis,...]

    # square and sum over the second axis, take the sqrt and transpose. the
    # result is an (m, n) array of weighted euclidean distances
    D = np.sqrt((wdiff*wdiff).sum(1)).T

    # Identify closest point in track for each star, add weight to this
    # minimum distance and save the weighted distance.
    func_dist = np.array(D).min(axis=1)*np.array(weights[0])

    func_score = sum(func_dist)
    print 'fast', func_score
    
    return func_score
        


    for t_ind in range(len(funcs)):
        # Call function that returns the score for a given track.
        track_score2 = track_distance_err_weight(t_ind)
        # Store the scores for each function/track into list.
        score2.append(track_score2)

        # Call function that returns the score for a given track.
        track_score3 = fast_wdist(t_ind)
        # Store the scores for each function/track into list.
        score3.append(track_score3)



    
    
    
    # Make plot.
    # figsize(x1, y1), GridSpec(y2, x2) --> To have square plots: x1/x2 = 
    # y1/y2 = 2.5 
    fig = plt.figure(figsize=(20, 30)) # create the top-level container
    gs = gridspec.GridSpec(12, 8)  # create a GridSpec object
    
    # Set global x limits for the parameters.
    span = max(score) - min(score)
    x_min, x_max = min(score)-span/100., min(score)+span*0.1    
    

    # Plot most probable members stars with best fir isochrone.
    ax0 = plt.subplot(gs[8:10, 0:2])
    plt.xlim(max(-0.9,min(data[0][0])-0.2), min(3.9, max(data[0][0])+0.2))
    plt.ylim(max(data[0][1])+1., min(data[0][1])-0.5)
    #Set axis labels
    plt.xlabel(r'$(C-T_1)$', fontsize=18)
    plt.ylabel(r'$T_1$', fontsize=18)
    # Add text box
    text1 = r'$E_{(B-V)} = %0.2f$' '\n' % params[best_func][2]
    text2 = r'$Age (Gyr) = %0.3f$' '\n' % params[best_func][1]
    text3 = r'$z = %0.3f$' '\n' % z_met
    text4 = r'$(m-M)_o = %0.2f$''\n' % params[best_func][3]
    text5 = r'$d (Kpc) = %0.2f$' '\n' % dist_kpc
    text6 = r'$score = %0.2f$' % min_score
    text = text1+text2+text3+text4+text5+text6
    plt.text(0.05, 0.05, text, transform = ax0.transAxes, 
             bbox=dict(facecolor='white', alpha=0.6), fontsize=12)
    # Set minor ticks and grid.
    ax0.minorticks_on()
    ax0.xaxis.set_major_locator(MultipleLocator(1.0))
    ax0.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    # Create new list with inverted values so higher prob stars are on top.
    m_p_m_temp = [data[0][0], data[0][1], weights[0]]
    cm = plt.cm.get_cmap('RdYlBu_r')
    # Invert values.
    m_p_m_temp_inv = [i[::-1] for i in m_p_m_temp]
    # Plot stars.
    plt.scatter(m_p_m_temp_inv[0], m_p_m_temp_inv[1], marker='o', 
                c=m_p_m_temp_inv[2], s=30, cmap=cm, lw=0.5, vmin=0, vmax=1)
    # Plot best fit isochrone.
    plt.plot(funcs[best_func][0], funcs[best_func][1], c='g',lw=1.8)
    plt.scatter(funcs[best_func][0], funcs[best_func][1], marker='o', color='k',
                s=2)
    

    # Plot most probable members stars with best fit isochrone.
#    ax5 = plt.subplot(gs[8:10, 2:6])
#    plt.ylabel(r'$score$', fontsize=18)
#    plt.xlabel(r'$tracks$', fontsize=18)
#    plt.xlim(0, len(funcs))
#    ax5.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
#    ebv_iso_temp = []
#    for isoch in params:
#        ebv_iso_temp.append(isoch[2])
#    cm = plt.cm.get_cmap('jet')
#    plt.scatter(range(len(funcs)), score, marker='o', s=1, c=ebv_iso_temp, cmap=cm) 
#    ax5.plot(best_func, min(score), 'or')
    
    

    # Plot for E(B-V) values.
    ax1 = plt.subplot(gs[10:12, 0:2])
    ebv_iso_temp = []
    for isoch in params:
        ebv_iso_temp.append(isoch[2])
    #Set axis labels
    plt.xlabel(r'$score$', fontsize=18)
    plt.ylabel(r'$E_{(B-V)}$', fontsize=18)
    plt.xlim(x_min, x_max)
    plt.ylim(e_bv_min-0.05, e_bv_max+0.05)
    delta_ebv = ebv_var[1]-ebv_var[0]
    erp = 100.*(delta_ebv)/params[best_func][2] if params[best_func][2] != 0.\
    else 100.*(delta_ebv)
    text1 = r'$re_{%d\%%} = %d \%%$' '\n' % (j_ebv, erp)
    text2 = r'$\Delta_{E_{(B-V)}}=(%0.2f;\, %0.2f)$' % (ebv_var[0], ebv_var[1])
    text = text1+text2
    plt.text(0.05, 0.9, text, transform = ax1.transAxes, 
             bbox=dict(facecolor='white', alpha=0.6), fontsize=12)
    # Set minor ticks and grid.
    ax1.minorticks_on()
#    ax1.xaxis.set_major_locator(MultipleLocator(2.0))
    ax1.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    plt.vlines(x=min_score+j_ebv*perc_err, ymin=ebv_var[0], ymax=ebv_var[1], 
               color='r', linestyles='dashed', linewidth=3)
    plt.scatter(score, ebv_iso_temp, s=10, edgecolor='None')
    ax1.plot(min_score, params[best_func][2], 'or')
    
    

    # Plot for age values.
    ax2 = plt.subplot(gs[10:12, 2:4])
    age_iso_temp = []
    for isoch in params:
        age_iso_temp.append(isoch[1])
    #Set axis labels
    plt.xlabel(r'$score$', fontsize=18)
    plt.ylabel(r'$Age (Gyr)$', fontsize=18)
    plt.xlim(x_min, x_max)
    plt.ylim(age_min-0.001, age_max+0.1)
    delta_age = age_var[1]-age_var[0]
    erp = 100.*(delta_age)/params[best_func][1]
    text1 = r'$re_{%d\%%} = %d \%%$' '\n' % (j_age, erp)
    text2 = r'$\Delta_{age}=(%0.3f;\, %0.3f)$' % (age_var[0], age_var[1])
    text = text1+text2
    plt.text(0.05, 0.9, text, transform = ax2.transAxes, 
             bbox=dict(facecolor='white', alpha=0.6), fontsize=12)
    # Set minor ticks and grid.
    ax2.minorticks_on()
#    ax2.xaxis.set_major_locator(MultipleLocator(2.0))
    ax2.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    plt.vlines(x=min_score+j_age*perc_err, ymin=age_var[0], ymax=age_var[1], 
               color='r', linestyles='dashed', linewidth=3)
    plt.scatter(score, age_iso_temp, s=10, edgecolor='None')
    ax2.plot(min_score, params[best_func][1], 'or')
    


    # Plot for z values.
    ax3 = plt.subplot(gs[10:12, 4:6])
    # Reformat list.
    dis_iso_temp = []
    for isoch in params:
        dis_iso_temp.append(isoch[0])
    #Set axis labels
    plt.xlabel(r'$score$', fontsize=18)
    plt.ylabel(r'$z$', fontsize=18)
    plt.xlim(x_min, x_max)
    plt.ylim(z_min-0.0005, z_max+0.005)
    delta_met = met_var[1]-met_var[0]
    erp = 100.*(delta_met)/params[best_func][0]
    text1 = r'$re_{%d\%%} = %d \%%$' '\n' % (j_met, erp)
    text2 = r'$\Delta_{z}=(%0.4f;\, %0.4f)$' % (met_var[0], met_var[1])
    text = text1+text2
    plt.text(0.05, 0.9, text, transform = ax3.transAxes, 
             bbox=dict(facecolor='white', alpha=0.6), fontsize=12)
    # Set minor ticks and grid.
    ax3.minorticks_on()
#    ax3.xaxis.set_major_locator(MultipleLocator(2.0))
    ax3.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    plt.vlines(x=min_score+j_met*perc_err, ymin=met_var[0], ymax=met_var[1], 
               color='r', linestyles='dashed', linewidth=3)
    plt.scatter(score, dis_iso_temp, s=10, edgecolor='None')
    ax3.plot(min_score, params[best_func][0], 'or')
    
    

    # Plot for distance modulus values.
    ax4 = plt.subplot(gs[10:12, 6:8])
    # Reformat list.
    dis_iso_temp = []
    for isoch in params:
        dis_iso_temp.append(isoch[3])
    #Set axis labels
    plt.xlabel(r'$score$', fontsize=18)
    plt.ylabel(r'$(m-M)_o$', fontsize=18)
    plt.xlim(x_min, x_max)
    plt.ylim(dis_mod_min-0.1, dis_mod_max+0.1)
    delta_dis = dis_var[1]-dis_var[0]
    erp = 100.*(delta_dis)/params[best_func][3] if params[best_func][3] != 0.\
    else 100.*(delta_dis)
    text1 = r'$re_{%d\%%} = %d \%%$' '\n' % (j_dis, erp)
    text2 = r'$\Delta_{d}=(%0.4f;\, %0.4f)$' % (dis_var[0], dis_var[1])
    text = text1+text2
    plt.text(0.05, 0.9, text, transform = ax4.transAxes, 
             bbox=dict(facecolor='white', alpha=0.6), fontsize=12)
    # Set minor ticks and grid.
    ax4.minorticks_on()
#    ax4.xaxis.set_major_locator(MultipleLocator(2.0))
    ax4.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    plt.vlines(x=min_score+j_dis*perc_err, ymin=dis_var[0], ymax=dis_var[1], 
               color='r', linestyles='dashed', linewidth=3)
    plt.scatter(score, dis_iso_temp, s=10, edgecolor='None')
    ax4.plot(min_score, params[best_func][3], 'or')


    fig.tight_layout()

   
    
print '\nEnd.'