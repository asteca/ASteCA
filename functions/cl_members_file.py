# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 17:04:26 2013

@author: gabriel
"""


def cluster_members_file(memb_file_out, memb_prob_avrg_sort):
    '''
    Create output data file with stars inside the cluster radius along with
    their membership probabilities.
    '''

    with open(memb_file_out, 'w') as out_data_file:
        out_data_file.write("#ID x y mag e_mag col1 e_col1 memb_prob\n")

    # Save average region obtained with regions used by the decont algor.
    with open(memb_file_out, "a") as f_out:
        for line in memb_prob_avrg_sort:
            f_out.write('{:<8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} \
{:>8}'.format(*line))
            f_out.write('\n')