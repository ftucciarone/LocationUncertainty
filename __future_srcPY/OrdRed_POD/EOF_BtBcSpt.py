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
import pathlib
import numpy as np
import snapshot_POD as pod
from netCDF4 import Dataset
import matplotlib.pyplot as plt


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
out_dir = work_dir + "/data/POD_Output/new_weight"
print(" isdir=", os.path.isdir(home_dir), "home_dir = ", home_dir)
print(" isdir=", os.path.isdir(work_dir), "work_dir = ", work_dir)
print(" isdir=", os.path.isdir(data_dir), "data_dir = ", data_dir)
print(" isdir=", os.path.isdir(dom_dir), " dom_dir = ", dom_dir)
print(" isdir=", os.path.isdir(out_dir), " out_dir = ", out_dir)
print("")

if len(sys.argv[1:])<4:
        print("                 ERROR: not enough arguments to run the script")
        print("Correct syntax is:",sys.argv[0], "-prexix PR -weight XXw1_XXw2")
        sys.exit(0)
# Set name of the input and output files
# --------------------------------------
# prefix: choice of the operation on the coarse-grained fluctuation (ke for rescaled)
prefix = sys.argv[2] #"keLR" #"keLR"   # Argument of first input
suffix = sys.argv[4] #"36w1_144w2"     # Argument of second input
infile = data_dir + "/" + prefix + "_CsGrained_r3_BtBcSpt_" + suffix + ".nc"
ingrid = dom_dir + "/domain_cfg_R3.nc"

print("isfile=", os.path.isfile(ingrid), "          Opening domain file:", ingrid)
print("isfile=", os.path.isfile(infile), "Opening coarse-grained fields:", infile)
print("")

# Define the type of inner product weighting
# ------------------------------------------
#   weighted: usual L2 product in R^3
# xyL2_zeucl: L2 in xy plane, euclidean on z
#  euclidean: not weighted in x,y,z
IPT = "weighted"
#IPT = "xyL2_zeucl"
IPT = "euclidean"

# Outputs
sfx = ["u", "v", "w"]
drc = ["zonal", "meridional", "vertical"]

# Options for POD procedures
single_file = 0
check_opt = 0

# Opening files
# -------------
fin = Dataset(infile, "r")
grid = Dataset(ingrid, "r")

# Get domain dimension from domain file
# -------------------------------------
nx = grid.variables["jpiglo"][:]
ny = grid.variables["jpjglo"][:]
nz = grid.variables["jpkglo"][:]

# Evaluation of the number of modes (i.e. of the times instants)
# --------------------------------------------------------------
nt = fin.variables["time_counter"][:].size
nm = nt - 1



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
elif IPT == "xyL2_zeucl":
    dVu = np.multiply(dAu[None, :, :], np.ones([nz, ny, nx]))
    # dVv = np.tile(dAv, (nz,1,1))
    # dVw = np.tile(dAw, (nz,1,1))
elif IPT == "euclidean":
    dVu = np.ones([nz, ny, nx]) * nx * ny * nz
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


# POD procedure
# ----------------------------
print("-------------------------------------------------------------")
print("                                              Collecting data")
u = fin.variables["uband"][:, :, :, :].data
v = fin.variables["vband"][:, :, :, :].data
w = fin.variables["wband"][:, :, :, :].data
u[:, -1, :, :] = 0
v[:, -1, :, :] = 0
w[:, -1, :, :] = 0

idt = 618
umean = fin.variables["uf2"][idt, :, :, :].data
vmean = fin.variables["vf2"][idt, :, :, :].data
wmean = fin.variables["wf2"][idt, :, :, :].data
umean [-1, :, :] = 0
vmean [-1, :, :] = 0
wmean [-1, :, :] = 0

print("                                      Computing mean velocity")
Uavg = np.mean(u, axis=0)
Vavg = np.mean(v, axis=0)
Wavg = np.mean(w, axis=0)

# %%
print("                                     Performing POD procedure")
u = u.reshape((nt, nx * ny * nz))
v = v.reshape((nt, nx * ny * nz))
w = w.reshape((nt, nx * ny * nz))
umean = umean.reshape((nx * ny * nz))
vmean = vmean.reshape((nx * ny * nz))
wmean = wmean.reshape((nx * ny * nz))

pc, ume, vme, wme, umo, vmo, wmo, tmode = \
    pod.spot_POD3F(u, v, w, wpu, wpu , wpu , nm, check_opt)
del u, v, w

# %%
ubias = umean.reshape((nz, ny, nx))
vbias = vmean.reshape((nz, ny, nx))
wbias = wmean.reshape((nz, ny, nx))

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
print("Preparing output file: ", out_dir + "/" + IPT + "/spmXXX" + "_BtBcSpt_" + suffix + "_" + IPT + ".nc")
for dim in range(3):
    outfile = out_dir + "/" + IPT + "/" + prefix + "spm" + sfx[dim] + "_BtBcSpt_" + suffix + "_" + IPT +".nc"
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

