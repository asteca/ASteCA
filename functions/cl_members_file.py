# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 17:04:26 2013

@author: gabriel
"""

from os.path import join

def cluster_members_file(output_dir, sub_dir, clust_name, membership_prob_list,
                         membership_prob_avrg_sort):
    '''
    Create output data file with most probable members for the cluster.
    '''
    
    # Create new file.
    memb_file = join(output_dir, sub_dir, clust_name+'_memb.dat')
    
    with open(memb_file, 'w') as out_data_file:
        out_data_file.write("#Region ID x y T1 e_T1 CT1 e_CT1 cont_index\n")
    
    # "a" opens the file for appending
    # Save all regions used by the decontamination algorithm.
    for region in membership_prob_list:
        with open(memb_file, "a") as f_out:
            for line in region:
                f_out.write('{:<3} {:>10} {:>10} {:>8} {:>8} {:>8} {:>8} {:>8} \
{:>8}'.format(*line))
                f_out.write('\n')
                
    # Save average region obtained with regions used by the decont algor.
    with open(memb_file, "a") as f_out:
        for line in membership_prob_avrg_sort:
            line_2 = [99.] + line
            f_out.write('{:<3} {:>10} {:>10} {:>8} {:>8} {:>8} {:>8} {:>8} \
{:>8}'.format(*line_2))
            f_out.write('\n')