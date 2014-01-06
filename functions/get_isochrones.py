# -*- coding: utf-8 -*-
"""
Created on Sun Aug  4 10:31:00 2013

@author: gabriel
"""

#from os import getcwd
from os.path import join
import numpy as np

def get_isochrones(mypath, clust_name):
    '''
    Get isochrone for a given cluster, moved by the corresponding E(B-V) and
    (m-M)o values.
    '''
    
    # Path to isochrones dir.
    iso_path = '/media/rest/github/isochrones/isochrones'
    
    # Get extinction, age, metallicity and dist module for this cluster from
    # the file that stores this data.
    myfile = 'clusters_data_isos.dat'
    with open(join(mypath, myfile), mode="r") as f_cl_dt:
       
        clust_found = False
        for line in f_cl_dt:
            li=line.split()
            # li[0]=name, [1]=E(B-V), [2]=Age, [3]=[Fe/H], [4]=(m-M)o
            if li[0] == clust_name:
                # Obtain E(B-V), Age, [Fe/H] and (m-M)o values for this cluster.
                cl_e_bv, cl_age, cl_feh, cl_dmod = float(li[1]), float(li[2]),\
                float(li[3]), float(li[4])
                # Change flag if cluster found in list.
                clust_found = True


    if clust_found:
        # If cluster was found in list.
        
        # Transform [Fe/H] to Z to match the isochrones files names.
        # A, B constants obtained through zunzun.com fitting of points:
        # z=(0.001,0.004,0.008,0.019), [Fe/H]=(-1.3,-0.7,-0.4,0.0) using a
        # y = a + b.log(x) fitting function.
        A, B = 1.7354259305164, 1.013629121876
        cl_z = round(10**((cl_feh-A)/B), 6)
        str_cl_z = str(cl_z)
        if cl_z == 0.001012:
            str_cl_z = '0.001'
        elif cl_z == 0.003957:
            str_cl_z = '0.004'
        elif cl_z == 0.007821:
            str_cl_z = '0.008'
        elif cl_z == 0.019405:
            str_cl_z = '0.019'
        
        # Obtain ages from the isochrone file with the corresponding metallicity.
        cl_file = str_cl_z+'.dat'
        # Open the isochrone for this cluster's metallicity.
        with open(join(iso_path, cl_file), mode="r") as f_iso:
    
            # List that will hold all the ages in the file.        
            ages = []
            # List that will hold all the isochrones.
            iso_list = [[[], []] for _ in range(150)]
    
            age_indx = -1
            for line in f_iso:
    
                li=line.strip()
                
                if li.startswith("#\tIsochrone\tZ = "):
                    age_indx += 1
                    for i, part in enumerate(line.split("Age =")):
                        if i==1:
                            # Store all ages in list.
                            ages.append(float(part[:-3])/1.e09)            
                
                # Save T1 and (C-T1) values for each isochrone.
                if not li.startswith("#"):
                    reader = li.split()
                    iso_list[age_indx][0].append(float(reader[9]))
                    iso_list[age_indx][1].append(float(reader[7]) -
                    float(reader[9]))
    
        # Get age in isochrone file closest to cluster's age.
        idx = (np.abs(np.array(ages)-cl_age)).argmin()
        # Update age variable with the closest value found.
        cl_age = min(ages, key=lambda x:np.abs(x-cl_age))
        
        # Isolate isochrone for this cluster.
        isochrone = [[], []]
        isochrone[0], isochrone[1] = iso_list[idx][0], iso_list[idx][1]
        
        # Move isochrone for this cluster given the E(B-V) and (m-M)o values.
        iso_moved = [[], []]
        for item in isochrone[0]:
            V_Mv = cl_dmod + 3.2*cl_e_bv
            iso_moved[0].append(item - 0.58*cl_e_bv + V_Mv)
        for item in isochrone[1]:
            iso_moved[1].append(item + 1.97*cl_e_bv)
            
    else:
        # If cluster was NOT found in list.
        iso_moved = [[0.], [0.]]
        # Set ZAMS values.
        cl_e_bv, cl_age, cl_feh, cl_dmod = 0.0, 0.0, 0.019, 18.90


    # Obtain ZAMS.
    cl_file = 'zams.dat'
    
    zams = [[], []]
    # Open the ZAMS isochrone.
    with open(join(iso_path, cl_file), mode="r") as f_zams:
        # Iterate through all lines.
        for line in f_zams:

            li=line.strip()
            
            # Save T1 and (C-T1) values for each isochrone.
            if not li.startswith("#"):
                reader = li.split()
                zams[0].append(float(reader[9]))
                zams[1].append(float(reader[7]) -
                float(reader[9]))
                
    # Move ZAMS. Use E(B-V)=0.0 and (m-M)o=18.90 if cluster was not found in
    # the list with assigned parameters.
    zams_iso = [[], []]
    for item in zams[0]:
        V_Mv = cl_dmod + 3.2*cl_e_bv
        zams_iso[0].append(item - 0.58*cl_e_bv + V_Mv)
    for item in zams[1]:
        zams_iso[1].append(item + 1.97*cl_e_bv)
            

    return cl_e_bv, cl_age, cl_feh, cl_dmod, iso_moved, zams_iso   