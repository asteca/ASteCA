# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 14:11:22 2013

@author: gabriel
"""

'''
Calculate the QQ-plot for the distribution of p-values obtained comparing
the cluster's KDE with the field region's KDEs.
'''


def ppoints(vector):
    '''
    Mimics R's function 'ppoints'.
    '''

    m_range = int(vector[0]) if len(vector)==1 else len(vector)
        
    n = vector[0] if len(vector)==1 else len(vector)
    a = 3./8. if n <= 10 else 1./2
         
    m_value =  n if len(vector)==1 else m_range
    pp_list = [((m+1)-a)/(m_value+(1-a)-a) for m in range(m_range)]
    
    return pp_list
    
    
def quantile():
    '''
    Mimics R's function 'quantile'.
    '''
    
    
    
def qqplot(p_vals_cl, p_vals_f):
    
    # Call ppoints function.
    B = ppoints(p_vals_cl)
    print B
    
    C = ppoints(p_vals_f)
    print C

    raw_input()