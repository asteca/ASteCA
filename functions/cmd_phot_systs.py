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
        ('ubvrijhk', [8, ['U', 'B', 'V', 'R', 'I', 'J', 'H', 'K']]),
        ('2mass', [8, ['J', 'H', 'Ks']]),
        ('washington', [8, ['C', 'M', 'T1', 'T2', 'B', 'V', 'R', 'I']])
        ]

    # Store as ordered dictionary.
    all_systs = odict(d_sys)

    return all_systs


def phot_wavelengths(sys, mag):
    '''
    Return filter effective wavelength for the magnitude passed.
    Values are given in inverse microns.
    '''

    all_systs = phot_mags()

    wave_dict = {
        # CCM 89, Table 3
        'ubvrijhk': (2.78, 2.27, 1.82, 1.43, 1.11, 0.8, 0.63, 0.46),
        '2mass': (),
        'washington': ()
        }

    # Get index of magnitude as stored in the dictionary.
    m_idx = all_systs[sys][1].index(mag)
    print mag, m_idx
    # Get effective wavelength for this magnitude.
    eff_wave = wave_dict[sys][m_idx]

    return eff_wave