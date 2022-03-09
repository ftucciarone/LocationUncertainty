#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 11:02:13 2021

         Computes the EOF through an Proper Orthogonal Decomposition 
         procedure starting of a 2D velocity field.
         Prints the output in 3 different files

@author: ftucciar
"""

# Load modules
# ------------
import numpy as np
from netCDF4 import Dataset
import sys
import snapshot_POD as pod
import matplotlib.pyplot as plt


# Set parammeters
# ---------------
# Inputs
base_dir = '/Users/ftucciar/LocationUncertainty/CoarseGraining/data/'
subs_dir = ['']
infile = '2Dvel.nc'
ingrid = 'domain_cfg_out.nc'
# Outputs
sfx = ['u', 'v', 'w']
drc = ['zonal', 'meridional', 'vertical']


single_file = 0
check_opt = 1

fin = Dataset(base_dir + subs_dir[0] + infile, 'r', format="NETCDF4")
grid = Dataset(base_dir + ingrid, 'r', format="NETCDF4")

# Get domain dimension from domain file
# -----------------------------------------------------------------------
nx = grid.variables['jpiglo'][:]
ny = grid.variables['jpjglo'][:]
nz = grid.variables['jpkglo'][:]

# Evaluation of the number of modes (i.e. of the times instants)
# -----------------------------------------------------------------------
nt = fin.dimensions['time'].size
nm = nt

fin.close()

# U-points (for integration) --------------------------------------------
e1u = grid.variables['e1u'][0, :, :].data
e2u = grid.variables['e2u'][0, :, :].data
# V-points (for integration) --------------------------------------------
e1v = grid.variables['e1v'][0, :, :].data
e2v = grid.variables['e2v'][0, :, :].data
# W-points (for integration) --------------------------------------------
e1w = grid.variables['e1t'][0, :, :].data
e2w = grid.variables['e2t'][0, :, :].data
# f-points (for integration) --------------------------------------------
e1f = grid.variables['e1f'][0, :, :].data
e2f = grid.variables['e2f'][0, :, :].data
v3 = np.multiply(e1f, e2f)
# Mid-point rule area at u and v points ---------------------------------
dAu = np.multiply(e1u, e2u)
dAv = np.multiply(e1v, e2v)
dAw = np.multiply(e1w, e2w)
# Free memory
del e1u, e2u, e1v, e2v, e1w, e2w

# Allocation of empty variables
# ----------------------------
wpu = np.diag(dAu.reshape((nx * ny)))
wpv = np.diag(dAv.reshape((nx * ny)))
# Free memory
del dAu, dAv, dAw

pc = np.empty((nm, nz))

umean = np.zeros((nz, ny, nx))
vmean = np.zeros((nz, ny, nx))
wmean = np.zeros((nz, ny, nx))

umode = np.zeros((nm, nz, ny, nx))
vmode = np.zeros((nm, nz, ny, nx))
wmode = np.zeros((nm, nz, ny, nx))

tmode = np.zeros((nm, nt, nz))  # eigenvectors
tmeco = np.zeros((nm, nz))      #

# POD procedure layer by layer
# ----------------------------
for k in range(nz-1):
    print()
    print('Layer = ', k+1)
    print('--------')
    print('Collecting data')
    # Allocate temporary matrices
    u = np.empty((nt, ny, nx))
    v = np.empty((nt, ny, nx))
    # Append velocity data from every .nc file
    file1 = base_dir + infile
    print('Opening file: ', file1)
    fin = Dataset(file1, 'r')
    u = fin.variables['udf'][:, k, :, :].data.astype('float64')
    v = fin.variables['vdf'][:, k, :, :].data.astype('float64')
    fin.close()
    print('Performing POD procedure')
    u = u.reshape((nt, nx * ny))
    v = v.reshape((nt, nx * ny))
    pc[:, k], ume, vme, umo, vmo, tmode[:, :, k] = \
        pod.spot_POD2F(u, v, wpu, wpv, check_opt)
    del u, v
    print('Project mean on spatial modes')
    tmeco[:, k] = pod.inner_prod(ume, umo, wpu) + \
        pod.inner_prod(vme, vmo, wpv)
    # Construction of the bias by projection
    ume = ume.reshape((ny, nx))
    vme = vme.reshape((ny, nx))
    umean[k, :, :] = np.matmul(umo.transpose(), tmeco[:, k]).reshape((ny, nx))
    vmean[k, :, :] = np.matmul(vmo.transpose(), tmeco[:, k]).reshape((ny, nx))
    # Multlipication of modes by eigenvalue
    umo = np.matmul(np.diag(np.sqrt(pc[:, k])), umo)
    vmo = np.matmul(np.diag(np.sqrt(pc[:, k])), vmo)

    umode[:, k, :, :] = umo.reshape((nm, ny, nx))
    vmode[:, k, :, :] = vmo.reshape((nm, ny, nx))
    del umo, vmo

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


# %% Saving .n files
# Save output data
# ----------------
print('Preparing output file: ', base_dir + 'spmXXXDF.nc')
for dim in range(3):
    outfile = base_dir + 'spm' + sfx[dim] + 'DF.nc'
    print('Writing outputs in file: ', outfile)
    fout = Dataset(outfile, 'w', format='NETCDF4')
    # Create dimensions
    fout.createDimension('y', ny)
    fout.createDimension('x', nx)
    fout.createDimension('z', nz)

    umid = fout.createVariable(sfx[dim] + '_mean', 'f4', ('z', 'y', 'x',))
    umid.units = 'm/s'
    if dim == 0:
        umid[:, :, :] = umean
    elif dim == 1:
        umid[:, :, :] = vmean
    else:
        umid[:, :, :] = wmean

    for i in range(nm):
        idx = 'spat_basis_' + sfx[dim] + '_' + str(i+1).zfill(3)
        uoid = fout.createVariable(idx, 'f4', ('z', 'y', 'x'))
        uoid.units = 'm/s'

        # Write data
        if dim == 0:
            uoid[:, :, :] = umode[i, :, :, :]
        elif dim == 1:
            uoid[:, :, :] = vmode[i, :, :, :]
        else:
            uoid[:, :, :] = wmode[i, :, :, :]

    # Close file
    fout.close()

# %% Saving .dat file
np.savetxt(base_dir + 'eigen_ModByLevsDF.dat', pc[:, :-1], delimiter=' ')
