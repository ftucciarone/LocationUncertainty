#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 11:02:13 2021

         Computes the EOF through an Proper Orthogonal Decomposition
         procedure starting of a 3D velocity field.
         Prints the output in 3 different files

@author: ftucciar
"""

# Load modules
# ------------
import numpy as np
import scipy.io.netcdf as nc
import sys
import snapshot_POD as pod
import matplotlib.pyplot as plt


# Set parammeters
# ---------------
# Inputs
base_dir = '/Users/ftucciar/LocationUncertainty/CoarseGraining/data/'
subs_dir = ['100-102y/', '102-104y/', '104-106y/', '106-108y/', '108-110y/']
infile = 'ocref_r3.nc'
ingrid = 'domain_cfg_out.nc'
# Outputs
sfx = ['u', 'v', 'w']
drc = ['zonal', 'meridional', 'vertical']


single_file = 0
check_opt = 1

fin = nc.netcdf_file(base_dir + subs_dir[0] + infile, 'r')
grid = nc.netcdf_file(base_dir + ingrid, 'r')

# Get domain dimension from domain file
# -----------------------------------------------------------------------
nx = grid.variables['jpiglo'].data
ny = grid.variables['jpjglo'].data
nz = grid.variables['jpkglo'].data

# Evaluation of the number of modes (i.e. of the times instants)
# -----------------------------------------------------------------------
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
        pod.spot_POD3F(u, v, w, wpu, wpv, wpw, check_opt)
    del u, v, w
    print('Project mean on spatial modes')
    tmeco[:, k] = pod.inner_prod(ume, umo, wpu) + \
        pod.inner_prod(vme, vmo, wpv) + \
        pod.inner_prod(wme, wmo, wpw)
    ume = ume.reshape((ny, nx))
    vme = vme.reshape((ny, nx))
    wme = wme.reshape((ny, nx))
    umean[k, :, :] = np.matmul(umo.transpose(), tmeco[:, k]).reshape((ny, nx))
    vmean[k, :, :] = np.matmul(vmo.transpose(), tmeco[:, k]).reshape((ny, nx))
    wmean[k, :, :] = np.matmul(wmo.transpose(), tmeco[:, k]).reshape((ny, nx))
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


# %% Saving .n files
# Save output data
# ----------------
for dim in range(3):
    outfile = base_dir + 'spm' + sfx[dim] + '.nc'
    print('Writing outputs in file: ', outfile)
    fout = nc.netcdf_file(outfile, 'w')
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
np.savetxt(base_dir + 'eigen_ModByLevs.dat', pc[:, :-1], delimiter=' ')
