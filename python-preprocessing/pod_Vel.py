#!/usr/bin/env python3

import os 
import copy
import numpy as np


import classDMD as dmd
import classPHY as phy
import classFILT as filt
import snapshot_POD as pod

from netCDF4 import Dataset


from os import environ
N_THREADS = '48'
environ['OMP_NUM_THREADS'] = N_THREADS
environ['OPENBLAS_NUM_THREADS'] = N_THREADS
environ['MKL_NUM_THREADS'] = N_THREADS
environ['VECLIB_MAXIMUM_THREADS'] = N_THREADS
environ['NUMEXPR_NUM_THREADS'] = N_THREADS

outfile = "trial.nc"
# %%
res = 3

# Path of R27 det. data
high_res = 27
dir_nemo = '/home/ftucciar/dataNEMO/Deter/'
dir_mode = 'R27d/'
dir_subs = ['100-102y/']#,'102-104y/','104-106y/','106-108y/','108-110y/']
filename = 'GYRE_5d_00010101_00021230'

dtype = 'float32'

# 
outfile = "spfilt_R" + str(high_res) + "toR" + str(3)

# %%



ratio= int(high_res/res)
flt_sigma = 0.25*ratio
dt = 5*86400
fc = 1/(60*86400)
levs = [i for i in range(31)]
nt = 0

# %%
if not os.path.isfile(outfile + ".nc"): 
    # 
    # If the given netCDF file does not exist, Filter the high resolution data
    infiles_uvel = [dir_nemo + dir_mode + subdir + filename + "_grid_U.nc" for subdir in dir_subs]
    infiles_vvel = [dir_nemo + dir_mode + subdir + filename + "_grid_V.nc" for subdir in dir_subs]
    infiles_wvel = [dir_nemo + dir_mode + subdir + filename + "_grid_W.nc" for subdir in dir_subs]

    uvel = np.zeros([0, len(levs), 20*high_res, 30*high_res])
    vvel = np.zeros([0, len(levs), 20*high_res, 30*high_res])
    wvel = np.zeros([0, len(levs), 20*high_res, 30*high_res])

    # Collect data from sequence of U files
    for i in range(len(infiles_uvel)):
        print(f'Reading data from {infiles_uvel[i]}')
        # U files
        nc = Dataset(infiles_uvel[i],'r')
        nt += nc.dimensions['time_counter'].size
        uvel = np.append(uvel, nc.variables['vozocrtx'][:,levs,1:-1, 1:-1].data.astype(dtype), axis=0)
        nc.close()

    # Collect data from sequence of V files
    for i in range(len(infiles_vvel)):
        print(f'Reading data from {infiles_vvel[i]}')
        # V files
        nc = Dataset(infiles_vvel[i],'r')
        vvel = np.append(vvel, nc.variables['vomecrty'][:,levs,1:-1, 1:-1].data.astype(dtype), axis=0) 
        nc.close()

    # Collect data from sequence of W files
    for i in range(len(infiles_vvel)):
        print(f'Reading data from {infiles_wvel[i]}')
        # W files
        nc = Dataset(infiles_wvel[i],'r')
        wvel = np.append(wvel, nc.variables['vovecrtz'][:,levs,1:-1, 1:-1].data.astype(dtype), axis=0) 
        nc.close()

    # Compute filtered fields and save them
    filt_vels = filt.sp_filt(nt, len(levs), high_res, res, 0.25, [{"u": uvel}, {"v": vvel}, {"w": wvel}])
    filt_vels.SpatialNetCDFSaveOut(outfile + ".nc")

    print(filt_vels)
    utopod = filt_vels.uflt["data"]
    vtopod = filt_vels.vflt["data"]
    wtopod = filt_vels.wflt["data"]

    uavg_ds = filt_vels.uavg_ds["data"]
    uavg_ds = filt_vels.vavg_ds["data"]
    wavg_ds = filt_vels.wavg_ds["data"]

# %%
if os.path.isfile(outfile + ".nc"): 
    #
    # Read file
    print("Reading data from " + outfile + ".nc")
    # W files
    nc = Dataset(outfile + ".nc",'r')
    utopod = nc.variables["uflt"][:, levs, 1:-1, 1:-1].data.astype(dtype)
    vtopod = nc.variables["vflt"][:, levs, 1:-1, 1:-1].data.astype(dtype)
    wtopod = nc.variables["wflt"][:, levs, 1:-1, 1:-1].data.astype(dtype)

    uavg_ds = nc.variables["uavg_ds"][levs, 1:-1, 1:-1].data.astype(dtype)
    vavg_ds = nc.variables["vavg_ds"][levs, 1:-1, 1:-1].data.astype(dtype)
    wavg_ds = nc.variables["wavg_ds"][levs, 1:-1, 1:-1].data.astype(dtype)


# %% Read R3 det for mean
infiles_uvel_R3 = ["/home/ftucciar/Stockage12T/September23/R9_det/EXP00/GYRE_dt_5d_01010101_01110101_grid_U.nc"]
infiles_vvel_R3 = ["/home/ftucciar/Stockage12T/September23/R9_det/EXP00/GYRE_dt_5d_01010101_01110101_grid_V.nc"]

infiles_uvel_R3 = ["/home/ftucciar/Stockage12T/These_Data/data_R3/R3_det/FULL/GYRE_det_5d_00010101_00151230_grid_U.nc"]
infiles_vvel_R3 = ["/home/ftucciar/Stockage12T/These_Data/data_R3/R3_det/FULL/GYRE_det_5d_00010101_00151230_grid_V.nc"]


nc = Dataset(infiles_uvel_R3[0],'r')
uavg_R3 = np.mean(nc.variables['vozocrtx'][:720,levs,1:-1, 1:-1].data.astype(dtype), axis=0) 
nc.close()

nc = Dataset(infiles_vvel_R3[0],'r')
vavg_R3 = np.mean(nc.variables['vomecrty'][:720,levs,1:-1, 1:-1].data.astype(dtype), axis=0) 
nc.close()

umean = uavg_ds - uavg_R3
vmean = vavg_ds - vavg_R3


nx = np.shape(utopod)[-1]
ny = np.shape(utopod)[-2] 
nz = np.shape(utopod)[-3] 
nm = np.shape(utopod)[-4]
nt = np.shape(utopod)[-4]
# %%



# Allocation of empty variables
# ----------------------------
print("-------------------------------------------------------------")
print("                                  Building 3D weight matrices")
#dVu = np.ones([nz, ny, nx])
#wpu = np.diag(dVu.reshape((nx * ny * nz))) / np.sum(dVu, axis=(-3,-2,-1))
#wpv = wpu # np.diag(dVv.reshape((nx * ny * nz))) / np.sum(dVv, axis=(-3,-2,-1))
#wpw = wpu # np.diag(dVw.reshape((nx * ny * nz))) / np.sum(dVw, axis=(-3,-2,-1))


# %%
pc = np.empty(nm)

# Initialise modes
umode = np.empty((nm, nz, ny, nx))
vmode = np.empty((nm, nz, ny, nx))
wmode = np.empty((nm, nz, ny, nx))

# Initialise temporal coefficient
tmode = np.empty((nm, nt))  # eigenvectors

# Initialise bias projection coefficient 
tmeco = np.empty((nm, nz))      #

# POD procedure
# ----------------------------
print("-------------------------------------------------------------")


print("                                       Removing mean velocity")
uprime = utopod - uavg_ds
vprime = vtopod - vavg_ds
wprime = wtopod - wavg_ds

# %%
print("                                     Performing POD procedure")
u = uprime.reshape((nt, nx * ny * nz))
v = vprime.reshape((nt, nx * ny * nz))
w = wprime.reshape((nt, nx * ny * nz))

pc, ume, vme, wme, umo, vmo, wmo, tmode = pod.spot_POD3F(u, v, w, nm, False)
del u, v, w


print("                         Multiplying modes by the eigenvalues")
# Multlipication of modes by eigenvalue
umo = np.matmul(np.diag(np.sqrt(pc)), umo)
vmo = np.matmul(np.diag(np.sqrt(pc)), vmo)
wmo = np.matmul(np.diag(np.sqrt(pc)), wmo)

umode = umo.reshape((nm, nz, ny, nx))
vmode = vmo.reshape((nm, nz, ny, nx))
wmode = wmo.reshape((nm, nz, ny, nx))
del umo, vmo, wmo


# %% Saving .n files
# Save output data
# ----------------
print("Preparing output file: ", "/" + outfile + "_spmXXX.nc")
sfx = ["u", "v", "w"]
for dim in range(3):
    ncout = outfile + "_spm" + sfx[dim] +  ".nc"
    print("Writing outputs in file: ", ncout)
    fout = Dataset(ncout, "w", format="NETCDF4")
    # Create dimensions
    fout.createDimension("y", ny+2)
    fout.createDimension("x", nx+2)
    fout.createDimension("z", nz)

    # Baroclinic Mean
    umid = fout.createVariable(sfx[dim] + "_mean", "f4", ("z", "y", "x",))
    umid.units = "m/s"
    if dim == 0:
        tmp = np.zeros([nz, ny+2, nx+2])
        tmp[...,1:-1,1:-1] = umean
        umid[:, :, :] = tmp
    elif dim == 1:
        tmp = np.zeros([nz, ny+2, nx+2])
        tmp[...,1:-1,1:-1] = vmean
        umid[:, :, :] = tmp
    else:
        umid[:, :, :] = np.zeros([nz, ny+2, nx+2])
        
    for i in range(nm):
        # Baroclinic
        idx = "spat_basis_" + sfx[dim] + "_" + str(i+1).zfill(3)
        uoid = fout.createVariable(idx, "f4", ("z", "y", "x"))
        uoid.units = "m/s"

        # Write data
        if dim == 0:
            tmp = np.zeros([nz, ny+2, nx+2])
            tmp[...,1:-1,1:-1] = umode[i, :, :, :]
            uoid[:, :, :] = tmp
        elif dim == 1:
            tmp = np.zeros([nz, ny+2, nx+2])
            tmp[...,1:-1,1:-1] = vmode[i, :, :, :]
            uoid[:, :, :] = tmp
        else:
            uoid[:, :, :] = np.zeros([nz, ny+2, nx+2])
            
    # Close file
    fout.close()

