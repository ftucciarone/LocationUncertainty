#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 11:02:13 2021

@author: ftucciar
"""

# Load modules
# ------------
import numpy as np
import scipy.io.netcdf as nc
import sys
from netCDF4 import Dataset
import matplotlib.pyplot as plt

# Set parammeters
# ---------------
# Inputs
base_dir = '/Users/ftucciar/LU_NEMO/'
subs_dir = ['TO_EXP_TUCCIA/']
basefile = 'spm'
ext = '.nc'
ingrid = '/Users/ftucciar/LU_prepocessing/CoarseGrain_def/domain_cfg_out.nc'

fU_in = Dataset(base_dir + subs_dir[0] + "/" + basefile + 'u' + ext, "r",
                format="NETCDF4")
fV_in = Dataset(base_dir + subs_dir[0] + "/" + basefile + 'v' + ext, "r",
                format="NETCDF4")
fW_in = Dataset(base_dir + subs_dir[0] + "/" + basefile + 'w' + ext, "r",
                format="NETCDF4")
grid = nc.netcdf_file(ingrid, 'r')

# Get domain dimension from domain file
# -------------------------------------
nx = grid.variables['jpiglo'].data
ny = grid.variables['jpjglo'].data
nz = grid.variables['jpkglo'].data

# U-points (for integration)
e1u = grid.variables['e1u'].data
e2u = grid.variables['e2u'].data
# V-points (for integration)
e1v = grid.variables['e1v'].data
e2v = grid.variables['e2v'].data
# W-points (for integration)
e1w = grid.variables['e1t'].data
e2w = grid.variables['e2t'].data
# U-points (for integration)
e1t = grid.variables['e1t'].data
e2t = grid.variables['e2t'].data
e3t = grid.variables['e3t_0'].data
# Mid-point rule area at u and v points
dAu = np.multiply(e1u, e2u)
dAv = np.multiply(e1v, e2v)
dAw = np.multiply(e1w, e2w)
# Free memory
del e1u, e2u, e1v, e2v, e1w, e2w

# Allocation of empty variables
# -----------------------------
wpu = np.diag(dAu.reshape((nx * ny)))
wpv = np.diag(dAv.reshape((nx * ny)))
wpw = np.diag(dAw.reshape((nx * ny)))
# Free memory
del dAu, dAv, dAw

# %%
# Read the fields
# ---------------

u = fU_in.variables['spat_basis_u_020'][:, :, :]
v = fV_in.variables['spat_basis_v_020'][:, :, :]
w = fW_in.variables['spat_basis_w_020'][:, :, :]



dxU = np.divide(u[:-1, :-1, 1:] - u[:-1, :-1, :-1],
                np.repeat(e1t[:, 1:, 1:], nz-1, axis=0))

dyV = np.divide(v[:-1, 1:, :-1] - v[:-1, :-1, :-1],
                np.repeat(e2t[:, 1:, 1:], nz-1, axis=0))

dzW = np.divide(w[:-1, 1:, 1:] - w[1:, 1:, 1:], e3t[0, :-1, 1:, 1:])

divU = dxU + dyV + dzW

print('Maximum divergence', np.abs(divU).max())