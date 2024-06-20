#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 31 10:57:18 2024

@author: ftucciar
"""
import os
import sys
import copy
import glob
import pathlib
import numpy as np
from netCDF4 import Dataset
from numpy import linalg as LA
import classUtils as uts
import classDMD as dmd
import classPHY as phy

user = "ftucciar"
rootdir = "/home/" + user + "/def_JAMES/"

filtrat = "48d1_54d2"
res = "9"


# Define input files
# ------------------
infile = rootdir + "srcF90/Output/keLR_CsGrained_r" + res + "_ZeroMean_" + filtrat + ".nc"
ingrid = rootdir + "/nemo-domaincfg/domain_cfg_R" + res + ".nc"
dmdoutf = rootdir + "output-DMD/keLR_aDMD_r" + res + "_ZeroMean_" + filtrat + ".nc"
diaoutf = rootdir + "output-DMD/keLR_diag_r" + res + "_ZeroMean_" + filtrat + ".nc"

# Define sampling time
# --------------------
dt = 5*86400

# Opening files
# -------------
fin = Dataset(infile, "r")
grid = Dataset(ingrid, "r")

# Get domain dimension from domain file
# --------------------------------------
nx = grid.variables['jpiglo'][:].data
ny = grid.variables['jpjglo'][:].data
nz = grid.variables['jpkglo'][:].data

# Evaluation of the number of modes (i.e. of the times instants)
# --------------------------------------------------------------
nt = fin.variables["time_counter"][:].size
nt = 360

# Define data-type collecting physiscal data
# ------------------------------------------
ids = "f2"
uvw3D = phy.Phys_data(nx, ny, nz, nt)
for vv in ["U", "V", "W"]:
    uvw3D.append(vv+"data", fin.variables[vv.lower() + ids][:nt, :, :, :].data )
    uvw3D.time_detrend(vv+"data",vv)

# Construct snapshot matrix
# -------------------------
X1 = np.concatenate(( uvw3D.unroll("Udata").T,
                      uvw3D.unroll("Vdata").T,
                      uvw3D.unroll("Wdata").T), axis=0) 
# Define the frequencies for plotting 
# -----------------------------------
freq_time = np.array([5, 10, 15, 20, 25, 30, 45, 60, 90, 120], dtype=float)

# Perform DMD and clean the result from exponential modes
# -------------------------------------------------------
DMD_vecStruct = dmd.DMD_vecmodes( uvw3D.attr, uvw3D.sizes, uvw3D.shapes, X1, 2000, dt )
dmd.PlotCircle(DMD_vecStruct.ct_eig, DMD_vecStruct.dt_eig, "(Base state)", freq_time=freq_time)
dmd.PlotAmplitudes(DMD_vecStruct.ct_eig, DMD_vecStruct.amp, "(Base state)" )

# %%
# Cleaning and Plotting 
# ---------------------
DMD_vecStruct.CleanUnitCirc(0.025)
DMD_vecStruct.CleanExpModes(0.005) 
dmd.PlotCircle(DMD_vecStruct.ct_eig, DMD_vecStruct.dt_eig, "(Cleaned)", freq_time=freq_time)
dmd.PlotAmplitudes(DMD_vecStruct.ct_eig, DMD_vecStruct.amp, "(Cleaned)" )

# %%
# Check for tresholding random/correlated modes 
# ---------------------------------------------
tt = 20
dmd.PlotCircle(DMD_vecStruct.ct_eig, DMD_vecStruct.dt_eig, "(Cleaned, $t=20$d)", freq_time=freq_time, split=tt )
dmd.PlotAmplitudes(DMD_vecStruct.ct_eig, DMD_vecStruct.amp, "(Cleaned, $t=20$d)", freq_time=freq_time[::2], split=tt )

# %%
# Gramiam rescaling procedure
# ---------------------------
prova = copy.deepcopy(DMD_vecStruct)
Psi = np.matmul( DMD_vecStruct.modes, LA.pinv(np.matmul(DMD_vecStruct.modes.conj().transpose(), DMD_vecStruct.modes)) )
prova.amp = np.matmul(Psi.conj().transpose(), X1[:,0].T)
dmd.PlotAmplitudes(prova.ct_eig, prova.amp, "(Rescaled, $t=20$d)", freq_time=freq_time[::2], split=tt )

# %%
# Splitting with threshold
# ------------------------
DMD_corr, DMD_uncorr = DMD_vecStruct.SplitByTreshold( tt )
diaoutf = diaoutf[:-3] + "_" + str(tt) + "d.nc"

#%%


Psi = np.matmul( DMD_corr.modes, LA.pinv(np.matmul(DMD_corr.modes.conj().transpose(), DMD_corr.modes)) )
DMD_corr.amp = np.matmul(Psi.conj().transpose(), X1[:,0].T)
dmd.PlotAmplitudes(DMD_corr.ct_eig, DMD_corr.amp, "(Rescaled, $t=20$d)", freq_time=freq_time[::2], split=tt )


Psi = np.matmul( DMD_uncorr.modes, LA.pinv(np.matmul(DMD_uncorr.modes.conj().transpose(), DMD_uncorr.modes)) )
DMD_uncorr.amp = np.matmul(Psi.conj().transpose(), X1[:,0].T)
dmd.PlotAmplitudes(DMD_uncorr.ct_eig, DMD_uncorr.amp, "(Rescaled, $t=20$d)", freq_time=freq_time[::2], split=tt )

# %%
DMD_corrphy = DMD_corr.ToSpatial()
DMD_uncorrphy = DMD_uncorr.ToSpatial()
dmd.SaveNetCDF(DMD_corrphy, DMD_uncorrphy, uvw3D.Utrnd, uvw3D.Vtrnd, diaoutf )

