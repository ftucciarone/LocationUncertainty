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
import scipy.io as sio
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os
import glob
import sys
import snapshot_POD as pod


# Set parammeters
# ---------------
# Inputs
base_dir = '/Volumes/LaCie/Nemo/Data_Tuccia/R27/'
subs_dir = ['100-102y/', '102-104y/', '104-106y/', '106-108y/', '108-110y/']
infile = 'data/ocref_r3.nc'
ingrid = 'domain_cfg_out.nc'
# Outputs
outfile = base_dir + 'oceof_r3_f2_t.nc'

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
    tmp = fin.dimensions['time']
    nt += tmp
    fin.close()
nm = nt

# Check for dimensions
# --------------------
if nx != fin.dimensions['xp']:
    sys.exit('x-grid dimension not matching between domain and input files')
if ny != fin.dimensions['yp']:
    sys.exit('y-grid dimension not matching between domain and input files')
if nz != fin.dimensions['z']:
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
# Mid-point rule area at u and v points
dAu = np.multiply(e1u, e2u)
dAv = np.multiply(e1v, e2v)
dAw = np.multiply(e1w, e2w)
# Free memory
del e1u, e2u, e1v, e2v, e1w, e2w

# Allocation of empty variables
# ----------------------------
wpu = np.diag(dAu.reshape((nx * ny)))
wpv = np.diag(dAv.reshape((nx * ny)))
wpw = np.diag(dAw.reshape((nx * ny)))
# Free memory
del dAu, dAv, dAw

pc = np.empty((nm, nz))

umean = np.empty((nz, ny, nx))
vmean = np.empty((nz, ny, nx))
wmean = np.empty((nz, ny, nx))

umode = np.empty((nm, nz, ny, nx))
vmode = np.empty((nm, nz, ny, nx))
wmode = np.empty((nm, nz, ny, nx))

tmode = np.empty((nm, nt, nz))  # eigenvectors
tmeco = np.empty((nm, nz))      #

# POD procedure layer by layer
# ----------------------------
for k in range(nz-1):
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
    print('Performing POD procedure')
    u = u.reshape((nt, nx * ny))
    v = v.reshape((nt, nx * ny))
    w = w.reshape((nt, nx * ny))
    pc[:, k], ume, vme, wme, umo, vmo, wmo, tmode[:, :, k] = \
        pod.spot_POD(u, v, w, wpu, wpv, wpw, check_opt)
    del u, v, w
    print('Project mean on spatial modes')
    tmeco[:, k] = pod.inner_prod(ume, umo, wpu) + \
        pod.inner_prod(vme, vmo, wpv) + \
        pod.inner_prod(wme, wmo, wpw)
    umean[k, :, :] = ume.reshape((ny, nx))
    vmean[k, :, :] = vme.reshape((ny, nx))
    wmean[k, :, :] = wme.reshape((ny, nx))
    del ume, vme, wme

    # Multlipication of modes by eigenvalue
    umo = np.matmul(np.diag(np.sqrt(pc[:, k])), umo)
    vmo = np.matmul(np.diag(np.sqrt(pc[:, k])), vmo)
    wmo = np.matmul(np.diag(np.sqrt(pc[:, k])), wmo)

    umode[:, k, :, :] = umo.reshape((nm, ny, nx))
    vmode[:, k, :, :] = vmo.reshape((nm, ny, nx))
    wmode[:, k, :, :] = wmo.reshape((nm, ny, nx))
    del umo, vmo, wmo

'''
# Compute spectrum of modes
# -------------------------
ric = np.cumsum(pc, axis=0) / np.sum(pc, axis=0)
for k in range(nz):
    plt.figure()
    plt.plot(range(nm), ric[:,k])
    plt.xlabel(r'Number of modes')
    plt.ylabel(r'Relative information content')
plt.show()
#nm = input("Please enter the number of modes to be saved:\n")
'''

# Save output data
# ----------------
print('Writing outputs in file: ', outfile)
fout = nc.netcdf_file(outfile, 'w')
# Create dimensions
fout.createDimension('time', nt)
fout.createDimension('mode', nt)
fout.createDimension('y', ny)
fout.createDimension('x', nx)
fout.createDimension('z', nz)
# Create variables
# tid = fout.createVariable('time', 'f4', ('time',))
mid = fout.createVariable('mode', 'i4', ('mode',))
# yid = fout.createVariable('y', 'f4', ('y',))
# xid = fout.createVariable('x', 'f4', ('x',))
# zid = fout.createVariable('z', 'f4', ('z',))
pcid = fout.createVariable('pc', 'f4', ('mode', 'z',))
tcid = fout.createVariable('tmeco', 'f4', ('mode', 'z',))

umid = fout.createVariable('umean', 'f4', ('z', 'y', 'x',))
vmid = fout.createVariable('vmean', 'f4', ('z', 'y', 'x',))
wmid = fout.createVariable('wmean', 'f4', ('z', 'y', 'x',))

uoid = fout.createVariable('umode', 'f4', ('mode', 'z', 'y', 'x',))
void = fout.createVariable('vmode', 'f4', ('mode', 'z', 'y', 'x',))
woid = fout.createVariable('wmode', 'f4', ('mode', 'z', 'y', 'x',))
toid = fout.createVariable('tmode', 'f4', ('mode', 'time', 'z',))
# Add attributes
mid.long_name = 'Mode index'
pcid.long_name = 'Principle Components'
pcid.units = 'm^2/s^2'

tcid.long_name = 'Projection coef. of mean on modes'
umid.long_name = 'Mean zonal velocity'
vmid.long_name = 'Mean meridional velocity'
wmid.long_name = 'Mean vertical velocity'
umid.units = 'm/s'
vmid.units = 'm/s'
wmid.units = 'm/s'
uoid.long_name = 'Zonal velocity modes'
void.long_name = 'Meridional velocity modes'
woid.long_name = 'Vertical velocity modes'
toid.long_name = 'Temporal modes'
uoid.units = 'm/s'
void.units = 'm/s'
woid.units = 'm/s'
# Write data
mid[:] = range(nt)
pcid[:, :] = pc
tcid[:, :] = tmeco
umid[:, :, :] = umean
vmid[:, :, :] = vmean
wmid[:, :, :] = wmean
uoid[:, :, :, :] = umode
void[:, :, :, :] = vmode
woid[:, :, :, :] = wmode
toid[:, :, :] = tmode
# Close file
fout.close()

