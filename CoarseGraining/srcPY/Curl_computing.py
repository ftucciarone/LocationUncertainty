#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 11:02:13 2021

         Prints the output in 3 different files

@author: ftucciar
"""

# Load modules
# ------------
import numpy as np
import scipy.io.netcdf as nc
import sys


# Set parammeters
# ---------------
# Inputs
base_dir = ''  # '/Volumes/LaCie/Nemo/Data_Tuccia/R27/'
subs_dir = ['100-102y/', '102-104y/', '104-106y/', '106-108y/', '108-110y/']
infile = 'ocref_r3.nc'
ingrid = 'domain_cfg_out.nc'
# Outputs
sfx = ['u', 'v', 'w']
wnd = ['uwnd', 'vwnd']
drc = ['zonal', 'meridional', 'vertical']


outfile = base_dir + 'oceof_r3_f2_twind.nc'

single_file = 0
check_opt = 1

fin = nc.netcdf_file(base_dir + subs_dir[0] + infile, 'r')
grid = nc.netcdf_file(ingrid, 'r')

# Get domain dimension from domain file
# -------------------------------------
nx = grid.variables['jpiglo'].data
ny = grid.variables['jpjglo'].data
nz = grid.variables['jpkglo'].data

# Evaluation of the number of modes (i.e. of the times instants)
# --------------------------------------------------------------
nt = 0
for s in subs_dir:
    file1 = base_dir + s + infile
    fin = nc.netcdf_file(file1, 'r')
    tmp = fin.dimensions['time_counter']
    nt += tmp
    fin.close()
nm = nt

# Check for dimensions
# --------------------
if nx != fin.dimensions['x']:
    sys.exit('x-grid dimension not matching between domain and input files')
if ny != fin.dimensions['y']:
    sys.exit('y-grid dimension not matching between domain and input files')
if nz != fin.dimensions['nav_lev']:
    sys.exit('z-grid dimension not matching between domain and input files')
fin.close()

# U-points (for integration)
e1u = grid.variables['e1u'].data
e2u = grid.variables['e2u'].data
# V-points (for integration)
e1v = grid.variables['e1v'].data
e2v = grid.variables['e2v'].data
# W-points (for integration)
e1w = grid.variables['e1t'].data
e2w = grid.variables['e2t'].data
# f-points (for integration)
e1f = grid.variables['e1f'].data
e2f = grid.variables['e2f'].data
v3 = np.multiply(e1f, e2f)
# Mid-point rule area at u and v points
dAu = np.multiply(e1u, e2u)
dAv = np.multiply(e1v, e2v)
dAw = np.multiply(e1w, e2w)
# Free memory
del e2u, e1v, e1w, e2w

# Allocation of empty variables
# ----------------------------
wpu = np.diag(dAu.reshape((nx * ny)))
wpv = np.diag(dAv.reshape((nx * ny)))
wpw = np.diag(dAw.reshape((nx * ny)))
# Free memory
del dAu, dAv, dAw

curl = np.empty((nt, nz, ny, nx))
    
# POD procedure layer by layer
# ----------------------------
for k in range(nz - 1):
    print()
    print('Layer = ', k+1)
    print('--------')
    print('Collecting data')
    # Allocate temporary matrices
    u = np.empty((0, ny, nx))
    v = np.empty((0, ny, nx))
    w = np.empty((0, ny, nx))
    # Append velocity data from every .nc file
    for s in subs_dir:
        file1 = base_dir + s + infile
        print('Opening file: ', file1)
        fin = nc.netcdf_file(file1, 'r')
        tmp = fin.variables['ur'][:, k, :, :].copy()
        u = np.append(u, tmp, axis=0)
        tmp = fin.variables['vr'][:, k, :, :].copy()
        v = np.append(v, tmp, axis=0)
        tmp = fin.variables['wr'][:, k, :, :].copy()
        w = np.append(w, tmp, axis=0)
        del tmp
        fin.close()
    print('Compute curl of 2D field')
    v1 = np.multiply(e2v, v)
    v2 = np.multiply(e1u, u)
    curl[:, k, :-1, :-1] = (v1[:, :-1, 1:] - v1[:, :-1, :-1]) - \
        (v2[:, 1:, :-1] - v2[:, :-1, :-1])
    curl[:, k, :, :] = np.divide(curl[:, k, :, :], v3)
    del v1, v2

# Save output curl
# ----------------
outfile = base_dir + 'curl.nc'
print('Writing outputs in file: ', outfile)
fout = nc.netcdf_file(outfile, 'w')
# Create dimensions
fout.createDimension('time', nt)
fout.createDimension('y', ny)
fout.createDimension('x', nx)
fout.createDimension('z', nz)
# Create variables
curlout = fout.createVariable('curl', 'f4', ('time', 'z', 'y', 'x',))
# Add attributes
curlout.long_name = 'Curl on f grid'
curlout.units = '1/s'
# Write data
curlout[:, :, :, :] = curl
# Close file
fout.close()

