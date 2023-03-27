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
import sys
import pathlib
import numpy as np
import snapshot_POD as pod
from netCDF4 import Dataset
import matplotlib.pyplot as plt


# Set parameters
# ---------------
# Output directory
home_dir = str(pathlib.Path.home())
OUTdir =  home_dir + "/LocationUncertainty/data/POD_Output/"

base_dir = home_dir + "/LocationUncertainty/data/"
subs_dir = ["100-102y/"] #, "102-104y/", "104-106y/", "106-108y/", "108-110y/"]
insuffix = "_0w1_36w2"
infile = "ocref_r3" + insuffix + ".nc"
dom_dir = home_dir + "/LocationUncertainty/data/"
ingrid = "domain_cfg_R3.nc"
# Outputs
sfx = ["u", "v", "w"]
drc = ["zonal", "meridional", "vertical"]
# Inner Product type
IPT = "weighted"
IPT = "euclidean"

single_file = 0
check_opt = 1

fin = Dataset(base_dir + subs_dir[0] + infile, "r")
grid = Dataset(dom_dir + ingrid, "r")

# Get domain dimension from domain file
# -----------------------------------------------------------------------
nx = grid.variables["jpiglo"][:]
ny = grid.variables["jpjglo"][:]
nz = grid.variables["jpkglo"][:]

# Evaluation of the number of modes (i.e. of the times instants)
# -----------------------------------------------------------------------
nt = 0
for s in subs_dir:
    file1 = base_dir + s + infile
    fin = Dataset(file1, "r")
    tmp = fin.variables["time_counter"][:].size
    nt += tmp
    fin.close()
nm = nt-1

# %%
print("-------------------------------------------------------------")

print("                                Collecting T-grid information")
e1t = grid.variables["e1t"][0, :, :].data
e2t = grid.variables["e2t"][0, :, :].data
dAt = np.multiply(e1t, e2t)
del e1t, e2t
e3t = grid.variables["e3t_0"][0, :, :].data
e3t_1d = grid.variables["e3t_1d"][0, :].data
depth = np.cumsum(np.append([5], e3t_1d[:-1]))

print("                                Collecting U-grid information")
e1u = grid.variables["e1u"][0, :, :].data
e2u = grid.variables["e2u"][0, :, :].data
dAu = np.multiply(e1u, e2u)
del e1u, e2u
e3u = grid.variables["e3u_0"][0, :, :, :].data

print("                                Collecting V-grid information")
e1v = grid.variables["e1v"][0, :, :].data
e2v = grid.variables["e2v"][0, :, :].data
dAv = np.multiply(e1v, e2v)
del e1v, e2v
e3v = grid.variables["e3v_0"][0, :, :].data

print("                                Collecting W-grid information")
# e3w = grid.variables["e3w_0"][0, :, :].data

# # f-points (for integration) --------------------------------------------
# e1f = grid.variables["e1f"][0, :, :].data
# e2f = grid.variables["e2f"][0, :, :].data
# v3 = np.multiply(e1f, e2f)


print("                                        Building volumes grid")
# # Mid-point rule area at u and v points ---------------------------------
if IPT == "weighted":
    dVu = np.multiply(dAu[None, :, :], e3u)
    # dVv = np.multiply(dAv[None, :, :], e3v)
    # dVw = np.multiply(dAw[None, :, :], e3w)
elif IPT == "euclidean":
    dVu = np.ones([nz, ny, nx])
    # dVv = np.tile(dAv, (nz,1,1))
    # dVw = np.tile(dAw, (nz,1,1))


# Allocation of empty variables
# ----------------------------
print("-------------------------------------------------------------")
print("                                  Building 3D weight matrices")
wpu = np.diag(dVu.reshape((nx * ny * nz))) / np.sum(dVu, axis=(-3,-2,-1))
#wpv = wpu # np.diag(dVv.reshape((nx * ny * nz))) / np.sum(dVv, axis=(-3,-2,-1))
#wpw = wpu # np.diag(dVw.reshape((nx * ny * nz))) / np.sum(dVw, axis=(-3,-2,-1))


# %%
pc = np.empty(nm)

# Initialize mean flow (Reynolds dec.) 
Uavg = np.empty((nz, ny, nx))
Vavg = np.empty((nz, ny, nx))
Wavg = np.empty((nz, ny, nx))

# Initialise fluctuations (Reynolds dec.)
Ufltz = np.empty((nt, ny * nx))
Vfltz = np.empty((nt, ny * nx))
Wfltz = np.empty((nt, ny * nx))

# Initialise modes
umode = np.empty((nm, nz, ny, nx))
vmode = np.empty((nm, nz, ny, nx))
wmode = np.empty((nm, nz, ny, nx))

# Initialise temporal coefficient
tmode = np.empty((nm, nt))  # eigenvectors

# Initialise bias projection coefficient 
tmeco = np.empty((nm, nz))      #

# Initialize bias term
ubias = np.empty((nz, ny, nx))
vbias = np.empty((nz, ny, nx))
wbias = np.empty((nz, ny, nx))

vTKEt = 0.0
vTKEz = np.zeros((nz))

mTKEt = 0.0
mTKEz = np.zeros((nz))

lTKEt = 0.0



# POD procedure
# ----------------------------
print("-------------------------------------------------------------")
print("                                              Collecting data")

# Allocate temporary matrices
u = np.empty((0, nz, ny, nx))
v = np.empty((0, nz, ny, nx))
w = np.empty((0, nz, ny, nx))

# Append velocity data from every .nc file
for s in subs_dir:
    file1 = base_dir + s + infile
    print("Opening file: ", file1)
    fin = Dataset(file1, "r")
    tmp = fin.variables["uband"][:, :, :, :].data
    u = np.append(u, tmp, axis=0)
    tmp = fin.variables["vband"][:, :, :, :].data
    v = np.append(v, tmp, axis=0)
    tmp = fin.variables["wband"][:, :, :, :].data
    tmp[:, -1, :, :] = 0
    w = np.append(w, tmp, axis=0)
    del tmp
    fin.close()
print("                                      Computing mean velocity")
Uavg = np.mean(u, axis=0)
Vavg = np.mean(v, axis=0)
Wavg = np.mean(w, axis=0)


print("                            Computing Baroclinic fluctuations")
Uflt = u - Uavg
Vflt = v - Vavg
Wflt = w - Wavg


print("                         Computin energy of the mean velocity")
Uavg = Uavg.reshape((nx * ny * nz)) 
Vavg = Vavg.reshape((nx * ny * nz))
Wavg = Wavg.reshape((nx * ny * nz))
avgKE = pod.inner_prod(Uavg, Uavg, wpu) + \
	pod.inner_prod(Vavg, Vavg , wpu) + \
	pod.inner_prod(Wavg, Wavg , wpu)




# # %%
# print("Computing energy of fluctuations (by layers)")
# for k in range(nz-1):
#     Ufltz = Uflt[:, k, :, :].reshape((nt, nx * ny))
#     Vfltz = Vflt[:, k, :, :].reshape((nt, nx * ny))
#     Wfltz = Wflt[:, k, :, :].reshape((nt, nx * ny))
#     vTKEz[k] = np.mean(np.diag(pod.inner_prod(Ufltz, Ufltz, wpuA)) + \
#         np.diag(pod.inner_prod(Vfltz, Vfltz, wpuA)) + \
#         np.diag(pod.inner_prod(Wfltz, Wfltz, wpuA)), axis=0)
# del Ufltz, Vfltz, Wfltz

# Uflt = Uflt.reshape((nt, nx * ny * nz))
# Vflt = Vflt.reshape((nt, nx * ny * nz))
# Wflt = Wflt.reshape((nt, nx * ny * nz))

# print("Computing Total energy of fluctuations")
# vTKEt = np.mean(np.diag(pod.inner_prod(Uflt, Uflt, wpu)) + \
#     np.diag(pod.inner_prod(Vflt, Vflt , wpv)) + \
#     np.diag(pod.inner_prod(Wflt, Wflt , wpw)), axis=0)
    
# %%
print("                          Performing Baroclinic POD procedure")
u = u.reshape((nt, nx * ny * nz))
v = v.reshape((nt, nx * ny * nz))
w = w.reshape((nt, nx * ny * nz))

pc, ume, vme, wme, umo, vmo, wmo, tmode = \
    pod.spot_POD3F(u, v, w, wpu, wpu , wpu , nm, check_opt)
del u, v, w

# %%
print("                     Project Baroclinic mean on spatial modes")
tmeco = pod.inner_prod(ume, umo, wpu) + \
    pod.inner_prod(vme, vmo, wpu) + \
    pod.inner_prod(wme, wmo, wpu)
ume = ume.reshape((nz, ny, nx))
vme = vme.reshape((nz, ny, nx))
wme = wme.reshape((nz, ny, nx))
ubias = np.matmul(umo.transpose(), tmeco)
vbias = np.matmul(vmo.transpose(), tmeco)
wbias = np.matmul(wmo.transpose(), tmeco)
del ume, vme, wme
biasKE = pod.inner_prod(ubias, ubias, wpu) + \
	pod.inner_prod(vbias, vbias , wpu) + \
	pod.inner_prod(wbias, wbias , wpu)
ubias = ubias.reshape((nz, ny, nx))
vbias = vbias.reshape((nz, ny, nx))
wbias = wbias.reshape((nz, ny, nx))


print("Scaling the bias with original kinetic energy of time average")
ubias = ubias * avgKE * biasKE
vbias = vbias * avgKE * biasKE
wbias = wbias * avgKE * biasKE

print("                         Multiplying modes by the eigenvalues")
# Multlipication of modes by eigenvalue
umo = np.matmul(np.diag(np.sqrt(pc)), umo)
vmo = np.matmul(np.diag(np.sqrt(pc)), vmo)
wmo = np.matmul(np.diag(np.sqrt(pc)), wmo)

umode = umo.reshape((nm, nz, ny, nx))
vmode = vmo.reshape((nm, nz, ny, nx))
wmode = wmo.reshape((nm, nz, ny, nx))
del umo, vmo, wmo

# %%
Ufltz = np.empty((nm, ny * nx))
Vfltz = np.empty((nm, ny * nx))
Wfltz = np.empty((nm, ny * nx))


print("-------------------------------------------------------------")
print("                        Computing energy of modes (by layers)")
for k in range(nz-1):
    Ufltz = umode[:, k, :, :].reshape((nm, nx * ny))
    Vfltz = vmode[:, k, :, :].reshape((nm, nx * ny))
    Wfltz = wmode[:, k, :, :].reshape((nm, nx * ny))
    mTKEz[k] = np.sum(np.diag(pod.inner_prod(Ufltz, Ufltz, wpuA) + \
        pod.inner_prod(Vfltz, Vfltz, wpuA) + \
        pod.inner_prod(Wfltz, Wfltz, wpuA)), axis=0)
del Ufltz, Vfltz, Wfltz

Uflt = umode.reshape((nm, nx * ny * nz))
Vflt = vmode.reshape((nm, nx * ny * nz))
Wflt = wmode.reshape((nm, nx * ny * nz))

print("                              Computing Total energy of modes")
mTKEt = np.sum(np.diag(pod.inner_prod(Uflt, Uflt, wpu)) + \
    np.diag(pod.inner_prod(Vflt, Vflt , wpv)) + \
    np.diag(pod.inner_prod(Wflt, Wflt , wpw)), axis=0)
lTKEt = np.sum(pc)


"""
# Compute spectrum of modes
# -------------------------
ric = np.cumsum(pc, axis=0) / np.sum(pc, axis=0)
for k in range(nz):
    plt.figure()
    plt.plot(range(nm), ric[:,k])
    plt.xlabel(r"Number of modes")
    plt.ylabel(r"Relative information content")
plt.show()
#nm = input("Please enter the number of modes to be saved:\n")
"""

# %% Saving .n files
# Save output data
# ----------------
print("Preparing output file: ", base_dir + "3Dpod/spmXXX" + insuffix + ".nc")
for dim in range(3):
    outfile = base_dir + "3Dpod/spm" + sfx[dim] + insuffix + "_" + IPT +".nc"
    print("Writing outputs in file: ", outfile)
    fout = Dataset(outfile, "w", format="NETCDF4")
    # Create dimensions
    fout.createDimension("y", ny)
    fout.createDimension("x", nx)
    fout.createDimension("z", nz)

    # Baroclinic Mean
    umid = fout.createVariable(sfx[dim] + "_mean", "f4", ("z", "y", "x",))
    umid.units = "m/s"
    if dim == 0:
        umid[:, :, :] = ubias
    elif dim == 1:
        umid[:, :, :] = vbias
    else:
        umid[:, :, :] = wbias
        
    for i in range(nm):
        # Baroclinic
        idx = "spat_basis_" + sfx[dim] + "_" + str(i+1).zfill(3)
        uoid = fout.createVariable(idx, "f4", ("z", "y", "x"))
        uoid.units = "m/s"

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
np.savetxt(base_dir + "3Dpod/eigen_ModByLevs3D" + insuffix + \
           "_" + IPT + ".dat", pc[:-1], delimiter=" ")
