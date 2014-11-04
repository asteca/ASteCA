# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 10:14:00 2014

@author: gabriel
"""

from collections import OrderedDict as odict


def phot_mags():
    '''
    Dictionary (ordered) that stores the names and column indexes for each
    magnitude defined in each photometric system, as presented in the CMD
    Girardi et al. files.
    '''

    # - 0: UBVRIJKH	   (cf. Maiz-Apellaniz 2006 + Bessell 1990)
    # - 1: JHKs        (2MASS)
    # - 2: CMT1T2BVRI  (Washington)

    d_sys = [
        ('ubvrijhk', ['U', 'B', 'V', 'R', 'I', 'J', 'H', 'K']),
        ('2mass', ['J', 'H', 'Ks']),
        ('washington', ['C', 'M', 'T1', 'T2', 'B', 'V', 'R', 'I'])
        ]

    # Store as ordered dictionary.
    all_systs = odict(d_sys)

    return all_systs


def phot_wavelengths(sys, mag):
    '''
    Return filter effective wavelength for the magnitude passed.
    Values are given in Armstrongs.
    '''

    all_systs = phot_mags()

    wave_dict = {
        # Effective wavelengths in Angstroms from Girardi CMD table.
        '2mass': (12329.79, 16395.59, 21522.05,),
        'ubvrijhk': (3641.89, 4460.62, 5501.70, 6557.09, 8036.57, 12314.46,
            16369.53, 21937.19),
        'washington': (3982.34, 5120.46, 6420.73, 8077.89, 4487.10, 5523.22,
            6557.09, 8036.57)
        }

    # Get index of magnitude as stored in the dictionary.
    m_idx = all_systs[sys].index(mag)
    # Get effective wavelength for this magnitude in inverse microns.
    eff_wave = 10000. / wave_dict[sys][m_idx]

    return eff_wave