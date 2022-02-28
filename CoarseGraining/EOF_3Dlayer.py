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
import snapshot_POD as pod


# Set parammeters
# ---------------
# Inputs
base_dir = '/Volumes/LaCie/Nemo/Data_Tuccia/R27/'
subs_dir = ['100-102y/', '102-104y/', '104-106y/', '106-108y/', '108-110y/']
infile = 'data/ocref_r3.nc'
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

# %% Ocean velocity POD
# Ocean velocity structures
# ----------------------------
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
for k in range(1):  # range(nz-1):
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

    if k == 0:
        val = pod.inner_prod(umo, umo, wpu) + \
            pod.inner_prod(vmo, vmo, wpv) + \
            pod.inner_prod(wmo, wmo, wpw)
        ans = np.diag(np.diag(val)) - np.diag(pc[:, k])


"""
    del umo, vmo, wmo

# %% Wind velocity POD
# Wind velocity structures
# ----------------------------
nm_wnd = nm
pcwnd = np.empty((nm))

uwnd_mean = np.empty((ny, nx))
vwnd_mean = np.empty((ny, nx))

uwnd_mode = np.empty((nm, ny, nx))
vwnd_mode = np.empty((nm, ny, nx))

twnd_mode = np.empty((nm, nt))  # eigenvectors
twnd_meco = np.empty((nm))      #

# POD procedure layer by layer
# ----------------------------
print()
print('Wind stress POD')
print('--------')
print('Collecting data')
# Allocate temporary matrices
uwnd = np.empty((0, ny, nx))
vwnd = np.empty((0, ny, nx))
# Append velocity data from every .nc file
for s in subs_dir:
    file1 = base_dir + s + infile
    print('Opening file: ', file1)
    fin = nc.netcdf_file(file1, 'r')
    tmp = fin.variables['uwndr'][:, :, :].copy()
    uwnd = np.append(uwnd, tmp, axis=0)
    tmp = fin.variables['vwndr'][:, :, :].copy()
    vwnd = np.append(vwnd, tmp, axis=0)
    del tmp
    fin.close()
print('Performing POD procedure')
uwnd = uwnd.reshape((nt, nx * ny))
vwnd = vwnd.reshape((nt, nx * ny))
pcwnd, uwnd_me, vwnd_me, uwnd_mo, vwnd_mo, twnd_mode = \
    pod.spot_POD2F(uwnd, vwnd, wpu, wpv, check_opt)
del uwnd, vwnd

# Reshaping array to neglect unphysical modes
pos_ind = np.where(pcwnd >= 0)
nm_wnd = pos_ind[-1][-1]
uwnd_mode = uwnd_mode[:nm_wnd, :, :]
vwnd_mode = vwnd_mode[:nm_wnd, :, :]

twnd_mode = twnd_mode[:nm_wnd, :]  # eigenvectors
twnd_meco = twnd_meco[:nm_wnd]     #

uwnd_mo = uwnd_mo[:nm_wnd, :]
vwnd_mo = vwnd_mo[:nm_wnd, :]

pcwnd = pcwnd[:nm_wnd]

print('Project mean on spatial modes')
twnd_meco = pod.inner_prod(uwnd_me, uwnd_mo, wpu) + \
    pod.inner_prod(vwnd_me, vwnd_mo, wpv)
uwnd_mean[:, :] = uwnd_me.reshape((ny, nx))
vwnd_mean[:, :] = vwnd_me.reshape((ny, nx))
del uwnd_me, vwnd_me

# Multlipication of modes by eigenvalue
uwnd_mo1 = np.matmul(np.diag(np.sqrt(pcwnd)), uwnd_mo)
uwnd_mo1 = np.matmul(np.diag(np.sqrt(pcwnd)), uwnd_mo)

uwnd_mode[:, :, :] = uwnd_mo.reshape((nm_wnd, ny, nx))
vwnd_mode[:, :, :] = vwnd_mo.reshape((nm_wnd, ny, nx))

del uwnd_mo, vwnd_mo
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
    # Create variables
    # yid = fout.createVariable('y', 'f4', ('y',))
    # xid = fout.createVariable('x', 'f4', ('x',))
    # zid = fout.createVariable('z', 'f4', ('z',))

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

    if dim < 2:
        umid = fout.createVariable(wnd[dim] + '_mean', 'f4', ('y', 'x',))
        umid.units = 'm/s'
        if dim == 0:
            umid[:, :] = uwnd_mean
        elif dim == 1:
            umid[:, :] = vwnd_mean

        for i in range(nm_wnd):
            idx = 'spat_basis_' + wnd[dim] + '_' + str(i+1).zfill(3)
            uoid = fout.createVariable(idx, 'f4', ('y', 'x'))
            uoid.units = 'm/s'

            # Write data
            if dim == 0:
                uoid[:, :] = uwnd_mode[i, :, :]
            elif dim == 1:
                uoid[:, :] = vwnd_mode[i, :, :]

    # Close file
    fout.close()

# %% Saving .dat file
np.savetxt(base_dir + 'eigen_ModByLevs.dat', pc[:, :-1], delimiter=' ')
np.savetxt(base_dir + 'eigen_ModWind.dat', pcwnd[:], delimiter=' ')
"""