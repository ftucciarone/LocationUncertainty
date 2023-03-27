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
import os
import sys
import glob
import pathlib
import numpy as np
from numpy import linalg as LA
import scipy.io as sio
import scipy.io.netcdf as nc
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from netCDF4 import Dataset
# Set directories
# ---------------
# home_dir: depending on the system
# work_dir: folder containing the LU procedures
# data_dir: containing the i/o data
#  dom_dir: containing the domain files from NEMO
#  out_dir: output parent directory
home_dir = str(pathlib.Path.home())
work_dir = home_dir + "/__future_LocationUncertainty"
data_dir = work_dir + "/data/CsGrain_Output"
dom_dir = work_dir + "/data"
out_dir = work_dir + "/data/DMD_Output"
print(" isdir=", os.path.isdir(home_dir), "home_dir = ", home_dir)
print(" isdir=", os.path.isdir(work_dir), "work_dir = ", work_dir)
print(" isdir=", os.path.isdir(data_dir), "data_dir = ", data_dir)
print(" isdir=", os.path.isdir(dom_dir), " dom_dir = ", dom_dir)
print(" isdir=", os.path.isdir(out_dir), " out_dir = ", out_dir)
print("")

if len(sys.argv[1:])<4:
        print("                 ERROR: not enough arguments to run the script")
        print("Correct syntax is:",sys.argv[0], "-prexix PR -weight TreatMean_XXw1_XXw2")
        sys.exit(0)
# Set name of the input and output files
# --------------------------------------
# prefix: choice of the operation on the coarse-grained fluctuation (ke for rescaled)
prefix = sys.argv[2] #"keLR" #"keLR"   # Argument of first input
suffix = sys.argv[4] #"36w1_144w2"     # Argument of second input
infile = data_dir + "/" + prefix + "_CsGrained_r9_" + suffix + ".nc"
ingrid = dom_dir + "/domain_cfg_R9.nc"

print("isfile=", os.path.isfile(ingrid), "          Opening domain file:", ingrid)
print("isfile=", os.path.isfile(infile), "Opening coarse-grained fields:", infile)
print("")

# %%
outfile = prefix + "_DMD_r9_" + suffix + ".nc"

# Outputs
sfx = ['u', 'v', 'w']
drc = ['zonal', 'meridional', 'vertical']

dt = 5*86400
single_file = 0
check_opt = 1
enpp = 0.99


# Opening files
# -------------
fin = Dataset(infile, "r")
grid = Dataset(ingrid, "r")

# Get domain dimension from domain file
# -----------------------------------------------------------------------
nx = grid.variables['jpiglo'][:]
ny = grid.variables['jpjglo'][:]
nz = grid.variables['jpkglo'][:]

# Evaluation of the number of modes (i.e. of the times instants)
# -----------------------------------------------------------------------
nt = fin.variables["time_counter"][:].size
nm = nt - 1

pc = np.empty((nm, nz))

umean = np.empty((nz, ny, nx))
vmean = np.empty((nz, ny, nx))
wmean = np.empty((nz, ny, nx))

Uavg = np.empty((nz, ny, nx))
Vavg = np.empty((nz, ny, nx))
Wavg = np.empty((nz, ny, nx))

Uavg = fin.variables["umean"][:, :, :].data
Vavg = fin.variables["vmean"][:, :, :].data
Wavg = fin.variables["wmean"][:, :, :].data
Uavg[-1, :, :] = 0
Vavg[-1, :, :] = 0
Wavg[-1, :, :] = 0



umode = np.empty((nm, nz, ny, nx))
vmode = np.empty((nm, nz, ny, nx))
wmode = np.empty((nm, nz, ny, nx))

tmode = np.empty((nm, nt, nz))  # eigenvectors
tmeco = np.empty((nm, nz))      #

# %%
U = np.empty((0, nz, ny, nx))
V = np.empty((0, nz, ny, nx))
W = np.empty((0, nz, ny, nx))
 
# Append velocity data from every .nc file
U = fin.variables["uband"][:, :, :, :].data
V = fin.variables["vband"][:, :, :, :].data
W = fin.variables["wband"][:, :, :, :].data
U[:, -1, :, :] = 0
V[:, -1, :, :] = 0
W[:, -1, :, :] = 0

Um = np.mean(U, axis=0, dtype=np.float64)
Vm = np.mean(V, axis=0, dtype=np.float64)
Wm = np.mean(W, axis=0, dtype=np.float64)
U -= Um[None] 
V -= Vm[None] 
W -= Wm[None] 

X1 = np.concatenate((U.reshape((nt,nz*ny*nx)).T,
                     V.reshape((nt,nz*ny*nx)).T,
                     W.reshape((nt,nz*ny*nx)).T), axis=0)

# POD
U1, S1, V1 = LA.svd(X1[:,:-1], full_matrices=False)

"""
ric = np.cumsum(S1, dtype=np.float64)/np.sum(S1, dtype=np.float64)
plt.figure()
plt.plot(ric)
plt.show()
idm = np.abs(ric - enpp).argmin()
print(r'%d EOFs is used to capture %4.3f of total energy'%(idm,enpp))
# Low-rank
U1 = U1[:,:idm]
S1 = S1[:idm]
V1 = V1[:,:idm]
nm = len(S1)
"""

# DMD
A1 = np.matmul(np.matmul(np.matmul(U1.transpose(), X1[:,1:]), V1), np.diag(1./S1))
lamb, W1 = LA.eig(A1) # 'lamb' is discrete-time eigenvalues
Phi = np.matmul(np.matmul(np.matmul(X1[:,1:], V1), np.diag(1./S1)), W1) # DMD modes
ct_lamb = np.log(lamb)/dt # continuous-time eigenvalues
b = np.matmul(LA.pinv(Phi), X1[:,0]) # DMD amptitude




# Create ncfile
#fout = nc.netcdf_file(base_dir + outfile, 'w')
fout = Dataset(out_dir + "/" + outfile, 'w', format='NETCDF4')

# Create dimensions
fout.createDimension('time', nt)
fout.createDimension('mode', nm)
fout.createDimension('y', ny)
fout.createDimension('x', nx)
fout.createDimension('z', nz)
# Create variables
tid = fout.createVariable('time', 'f8', ('time',))
mid = fout.createVariable('mode', 'i4', ('mode',))
yid = fout.createVariable('y', 'f8', ('y',))
xid = fout.createVariable('x', 'f8', ('x',))
zid = fout.createVariable('z', 'f8', ('z',))

lreid = fout.createVariable('lambda_real', 'f8', ('mode',))
limid = fout.createVariable('lambda_imag', 'f8', ('mode',))

ctlreid = fout.createVariable('ct_lamb_real', 'f8', ('mode',))
ctlimid = fout.createVariable('ct_lamb_imag', 'f8', ('mode',))

breid = fout.createVariable('b_real', 'f8', ('mode',))
bimid = fout.createVariable('b_imag', 'f8', ('mode',))

umid = fout.createVariable('umean', 'f8', ('z','y','x',))
vmid = fout.createVariable('vmean', 'f8', ('z','y','x',))
wmid = fout.createVariable('wmean', 'f8', ('z','y','x',))

ureid = fout.createVariable('umode_real', 'f8', ('mode','z','y','x',))
uimid = fout.createVariable('umode_imag', 'f8', ('mode','z','y','x',))

vreid = fout.createVariable('vmode_real', 'f8', ('mode','z','y','x',))
vimid = fout.createVariable('vmode_imag', 'f8', ('mode','z','y','x',))

wreid = fout.createVariable('wmode_real', 'f8', ('mode','z','y','x',))
wimid = fout.createVariable('wmode_imag', 'f8', ('mode','z','y','x',))
# Add attributes
tid.long_name = 'Time axis'
tid.units = 'years'
mid.long_name = 'Mode index'
yid.long_name = 'Ocean Y axis (T-grid)'
yid.units = 'km'
xid.long_name = 'Ocean X axis (T-grid)'
xid.units = 'km'
zid.long_name = 'Ocean mid-layer axis'
zid.units = 'km'
umid.long_name = 'Mean of zonal velocity'
umid.units = 'm/s'
vmid.long_name = 'Mean of meridional velocity'
vmid.units = 'm/s'
wmid.long_name = 'Mean of vertical velocity'
wmid.units = 'm/s'
lreid.long_name = 'DMD eigenvalues (real part)'
limid.long_name = 'DMD eigenvalues (imag part)'
ctlreid.long_name = 'Continuous time eigenvalues (real part)'
ctlimid.long_name = 'Continuous time eigenvalues (imag part)'
breid.long_name = 'DMD amplitudes (real part)'
bimid.long_name = 'DMD amplitudes (imag part)'
ureid.long_name = 'Zonal DMD modes (real part)'
vreid.long_name = 'Meridional DMD modes (real part)'
wreid.long_name = 'Vertical DMD modes (real part)'
uimid.long_name = 'Zonal DMD modes (imag part)'
vimid.long_name = 'Meridional DMD modes (imag part)'
wimid.long_name = 'Vertical DMD modes (imag part)'
# Write data
tid[:] = range(nt)
mid[:] = range(nm)
yid[:] = range(ny)
xid[:] = range(nx)
zid[:] = range(nz)
lreid[:] = lamb.real
limid[:] = lamb.imag
ctlreid[:] = ct_lamb.real
ctlimid[:] = ct_lamb.imag
breid[:] = b.real
bimid[:] = b.imag
umid[:,:,:] = Uavg
vmid[:,:,:] = Vavg
wmid[:,:,:] = Wavg
ureid[:,:,:,:] = (Phi[:nx*ny*nz,:].real.transpose()).reshape((nm,nz,ny,nx))
vreid[:,:,:,:] = (Phi[nx*ny*nz:2*nx*ny*nz,:].real.transpose()).reshape((nm,nz,ny,nx))
wreid[:,:,:,:] = (Phi[2*nx*ny*nz:,:].real.transpose()).reshape((nm,nz,ny,nx))
uimid[:,:,:,:] = (Phi[:nx*ny*nz,:].imag.transpose()).reshape((nm,nz,ny,nx))
vimid[:,:,:,:] = (Phi[nx*ny*nz:2*nx*ny*nz,:].imag.transpose()).reshape((nm,nz,ny,nx))
wimid[:,:,:,:] = (Phi[2*nx*ny*nz:,:].imag.transpose()).reshape((nm,nz,ny,nx))
# Close file
fout.close()

