# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 17:04:26 2013

@author: gabriel
"""

from os.path import join


def cluster_members_file(output_dir, sub_dir, clust_name, memb_prob_avrg_sort):
    '''
    Create output data file with stars inside the cluster radius along with
    their membership probabilities.
    '''

    # Create new file.
    memb_file = join(output_dir, sub_dir, clust_name + '_memb.dat')

    with open(memb_file, 'w') as out_data_file:
        out_data_file.write("#ID x y T1 e_T1 CT1 e_CT1 memb_prob\n")

    # Save average region obtained with regions used by the decont algor.
    with open(memb_file, "a") as f_out:
        for line in memb_prob_avrg_sort:
            f_out.write('{:<8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} \
{:>8}'.format(*line))
            f_out.write('\n')