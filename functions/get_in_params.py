# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 14:03:44 2014

@author: gabriel
"""

from os.path import join


def get_in_params(mypath):
    '''
    This function reads the input data parameters stored in the 'ocaat_input.dat'
    file and returns them packaged for each function to use.
    '''
    
# Allows to work with columns data files.

    data_file = join(mypath, 'ocaat_input.dat')
    
    with open(data_file, mode="r") as f_dat:
        
        # Iterate through each line in the file.
        for line in f_dat:
            
            if not line.startswith("#") and line.strip() != '':
                reader = line.split()
                
                # Read folder paths where clusters are stored.
                if reader[0] == 'MO':
                    mode = str(reader[1])
                elif reader[0] == 'CP0':
                    mypath2 = str(reader[1])
                elif reader[0] == 'CP1':
                    mypath3 = str(reader[1])
                elif reader[0] == 'CP2':
                    output_dir = str(reader[1])
                elif reader[0] == 'PD':
                    gd_params = map(int, reader[1:])
                elif reader[0] == 'CC':
                    gc_params = map(float, reader[1:])
                elif reader[0] == 'BR':
                    br_params = map(float, reader[1:])
                elif reader[0] == 'CR':
                    cr_params = map(float, reader[1:])
                elif reader[0] == 'ER':
                    er_params = map(float, reader[1:])
                elif reader[0] == 'GR':
                    gr_params = map(int, reader[1:])
                elif reader[0] == 'PV':
                    pv0_params = True if reader[1] == 'True' else False
                    pv1_params = str(reader[2])
                    pv2_params = int(reader[3])
                elif reader[0] == 'DA':
                    da0_params = str(reader[1])
                    da1_params = int(reader[2])
                    da2_params = int(reader[3])
                elif reader[0] == 'PS_p':
                    iso_path = str(reader[1])
                elif reader[0] == 'PS_l':
                    line_start = str(reader[1])
                elif reader[0] == 'PS_i':
                    indx_cols = map(int, reader[1:])
                elif reader[0] == 'PS_m':
                    m_rs = map(float, reader[1:])
                elif reader[0] == 'PS_a':
                    a_rs = map(float, reader[1:])
                elif reader[0] == 'PS_e':
                    e_rs = map(float, reader[1:])
                elif reader[0] == 'PS_d':
                    d_rs = map(float, reader[1:])
                elif reader[0] == 'PS_xy':
                    axes_names = map(str, reader[1:])

                elif reader[0] == 'SC':
                    sys_sel = str(reader[1])
                    IMF_name = str(reader[2])
                    tot_mass = float(reader[3])
                    f_bin = float(reader[4])
                    q_bin = float(reader[5])
                    
                    
    in_dirs = [mypath2, mypath3, output_dir]
    pv_params = [pv0_params, pv1_params, pv2_params]
    da_params = [da0_params, da1_params, da2_params, mypath2]
    ps_params = [iso_path, line_start, indx_cols, m_rs, a_rs, e_rs, d_rs,
                 axes_names]
    sc_params = [sys_sel, IMF_name, tot_mass, f_bin, q_bin]
    
    return mode, in_dirs, gd_params, gc_params, br_params, cr_params, er_params,\
    gr_params, pv_params, da_params, ps_params, sc_params