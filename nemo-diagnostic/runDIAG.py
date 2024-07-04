#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 12:37:07 2024

@author: ftucciar
"""

from classDIAG import DIAG_NEMO

# %%
gridfile = "/home/ftucciar/def_JAMES/nemo-domaincfg/domain_cfg_R3.nc"
basefile = "/home/ftucciar/Stockage12T/These_Data/data_R3/R3_dmd/FULL/GYRE_dmd_5d_00010101_00151230"

outfile = basefile + "_diagnostic.nc"

# Set param
param = {
    'dl': 72*5, # start delay
    'dt': 5.*86400, # deta step [s]
    'nb': 4*3+1, # pts. near boundary (including boundary)
    'infiles_grid' : gridfile, 
    'infiles_grid_U': [basefile + "_grid_U.nc"],
    'infiles_grid_V': [basefile + "_grid_V.nc"],
    'infiles_grid_W': [basefile + "_grid_W.nc"],
    'infiles_grid_T': [basefile + "_grid_T.nc"],
    'outfile': outfile,
    'dtype': 'float32', # 'float32' or 'float64'
}   

dia = DIAG_NEMO(param)
