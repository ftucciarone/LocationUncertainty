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
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt

# Operations
from OpenYAML import yaml2dict
from Operations import ke_compute
from mask import findOverlap
import DiagnosticMetrics as dgms

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Helvetica"],
    "figure.dpi": 300
    })

# %% Read dictionary of paths
# ---------------------------
HRfiles, LRfiles, IO_dict = yaml2dict(r'paths.yaml')

# %% Read Low resolution and High resolution domain files
# -------------------------------------------------------
LRgrid = Dataset(IO_dict["dom_dir"] + "/domain_cfg_R3.nc", "r", format="NETCDF4")
HRgrid = Dataset(IO_dict["dom_dir"] + "/domain_cfg_R27.nc", "r", format="NETCDF4")

# High resolution
HRlon = HRgrid.variables['nav_lon'][:, :].data
HRlat = HRgrid.variables['nav_lat'][:, :].data
HRe1 = HRgrid.variables['e1t'][0, :, :].data
HRe2 = HRgrid.variables['e2t'][0, :, :].data
HRe3 = HRgrid.variables['e3t_0'][0, :, :].data
HRdA = np.multiply(HRe1, HRe2); del HRe1, HRe2
HRdV = np.multiply(HRdA, HRe3); del HRe3
HRny = HRgrid.dimensions["y"].size
HRnx = HRgrid.dimensions["x"].size

# Low resolution
LRlon = LRgrid.variables['nav_lon'][:, :].data
LRlat = LRgrid.variables['nav_lat'][:, :].data
LRe1 = LRgrid.variables['e1t'][0, :, :].data
LRe2 = LRgrid.variables['e2t'][0, :, :].data
LRe3 = LRgrid.variables['e3t_0'][0, :, :].data
LRdA = np.multiply(LRe1, LRe2); del LRe1, LRe2
LRdV = np.multiply(LRdA, LRe3); del LRe3
LRny = LRgrid.dimensions["y"].size
LRnx = LRgrid.dimensions["x"].size

# Compute integration weights
# HRdA = np.ones_like(HRdA)
# HRdV = np.ones_like(HRdV)
# LRdA = np.ones_like(LRdA)
# LRdV = np.ones_like(LRdV)
r1_HRdA = 1 / np.sum(HRdA, axis=(-2, -1))
r1_HRdV = 1 / np.sum(HRdV, axis=(-3, -2, -1))
r1_LRdA = 1 / np.sum(LRdA, axis=(-2, -1))
r1_LRdV = 1 / np.sum(LRdV, axis=(-3, -2, -1))

# %% Find overlap of the two grids for downsampling purposes
cols, rows = findOverlap(HRlat, HRlon, LRlat, LRlon)

# %% Read corresponding to the High resolution models
# ------------------------------------------------
HRdata = Dataset(HRfiles[0]["uvel_files"][0], "r", format="NETCDF4")
HRvardims = HRdata.variables["vozocrtx"].dimensions
if np.size(HRvardims) == 3:
    HRvar = np.empty([0, 542, 812])
    HRcmvar_1m = np.zeros([542, 812])  # first moment cumulated
    HRcmvar_2m = np.zeros([542, 812])  # second moment cumulated
elif np.size(HRvardims) == 4:
    HRvar = np.empty([0, 31, 542, 812])
    HRcmvar_1m = np.zeros([31, 542, 812])  # first moment cumulated
    HRcmvar_2m = np.zeros([31, 542, 812])  # second moment cumulated

nt = 0
HRfct = np.empty((0,31,542,812))
for i in range(len(HRfiles[0]["curl_files"])):
    print(HRfiles[0]["curl_files"][i])
    if np.size(HRvardims) == 3:
        fin = Dataset(HRfiles[0]["curl_files"][i], "r", format="NETCDF4")
        temp = fin.variables["socurlt"][:, :, :]
        HRvar = np.append(HRvar, temp, axis=0)
        nt += fin.dimensions[HRvardims[0]].size
        HRcmvar_1m += np.sum(HRvar, axis=0)
        HRcmvar_2m += np.sum(HRvar**2, axis=0)
        # del HRvar
    elif np.size(HRvardims) == 4:
        fin = Dataset(HRfiles[0]["curl_files"][i], "r", format="NETCDF4")
        temp = fin.variables["socurlt"][:, :, :]
        HRvar = np.append(HRvar, temp, axis=0)
        HRcmvar_1m += np.sum(HRvar, axis=0)
        HRcmvar_2m += np.sum(HRvar**2, axis=0)
        # del HRvar
        nt += fin.dimensions[HRvardims[0]].size
    fin.close    
    
    
# %%
HRvar_1m = HRcmvar_1m / nt # First moment
HRvar_2m = HRcmvar_2m / nt # Second moment

HRavg = HRvar_1m; del HRvar_1m             # Mean flow
HRfct = HRvar - HRavg                      # Fluctuations
HRvnc = HRvar_2m - HRavg**2; del HRvar_2m  # Variance
HRstd = np.sqrt(HRvnc)                     # Standard deviation

HRavg_ds = np.delete(np.delete(HRavg, cols.astype(int), axis=-1), rows.astype(int), axis=-2)
HRfct_ds = np.delete(np.delete(HRfct, cols.astype(int), axis=-1), rows.astype(int), axis=-2)
HRvnc_ds = np.delete(np.delete(HRvnc, cols.astype(int), axis=-1), rows.astype(int), axis=-2)
HRstd_ds = np.delete(np.delete(HRstd, cols.astype(int), axis=-1), rows.astype(int), axis=-2)

# %%
HRavg_norm = np.sum(np.multiply(HRavg_ds, LRdV), axis=(-3,-2,-1))
HRvnc_norm = np.sum(np.multiply(HRvnc_ds, LRdV), axis=(-3,-2,-1))
HRfct_norm = np.mean(np.sum(np.multiply(HRfct_ds, LRdV), axis=(-3,-2,-1)))


# %%Read corresponding to the Low resolution models
# -----------------------------------------------
LRdata = Dataset(LRfiles[0]["curl_files"][0], "r", format="NETCDF4")
LUdata = Dataset(LRfiles[1]["curl_files"][0], "r", format="NETCDF4")

minTime = 0
maxTime = 720
idk = 1
# Get domain dimension from domain file
# -------------------------------------
LRvardims = LRdata.variables["vorticity"].dimensions
if np.size(LRvardims) == 3:
    LRvar = LRdata.variables["vorticity"][minTime:maxTime, :, :]
    LUvar = LUdata.variables["vorticity"][minTime:maxTime, :, :]
elif np.size(LRvardims) == 4:
    LRvar = LRdata.variables["vorticity"][minTime:maxTime, :, :, :]
    LUvar = LUdata.variables["vorticity"][minTime:maxTime, :, :, :]
            

LRvar_1m = np.mean(LRvar, axis=0)
LUvar_1m = np.mean(LUvar, axis=0)
LRvar_2m = np.mean(LRvar**2, axis=0) 
LUvar_2m = np.mean(LUvar**2, axis=0) 


LRavg = LRvar_1m                  # Mean flow
LRfct = LRvar - LRavg             # Fluctuations
LRvnc = LRvar_2m - LRavg**2       # Variance
LRstd = np.sqrt(LRvnc)            # Standard deviation
del LRvar_1m, LRvar_2m

LUavg = LUvar_1m                  # Mean flow
LUfct = LUvar - LUavg             # Fluctuations
LUvnc = LUvar_2m - LUavg**2       # Variance
LUstd = np.sqrt(LUvnc)            # Standard deviation
del LUvar_1m, LUvar_2m

LRavg_norm = np.sum(np.multiply(LRavg, LRdV), axis=(-3,-2,-1))
LUavg_norm = np.sum(np.multiply(LUavg, LRdV), axis=(-3,-2,-1))

LRvnc_norm = np.sum(np.multiply(LRvnc, LRdV), axis=(-3,-2,-1))
LUvnc_norm = np.sum(np.multiply(LUvnc, LRdV), axis=(-3,-2,-1))

LRfct_norm = np.mean(np.sum(np.multiply(LRfct, LRdV), axis=(-3,-2,-1)))
LUfct_norm = np.mean(np.sum(np.multiply(LUfct, LRdV), axis=(-3,-2,-1)))


# %% Compute and print RMS(avg)
HRavg_rms = dgms.RootMeanSquare(HRavg, HRdV, r1_HRdV)
LRavg_rms = dgms.RootMeanSquare(LRavg, LRdV, r1_LRdV)
LUavg_rms = dgms.RootMeanSquare(LUavg, LRdV, r1_LRdV)
print("RMS(avg),                  HR,                  LR,                  LU")
print("         ", HRavg_rms, LRavg_rms, LUavg_rms)


# %%
HRfct_rms = dgms.RootMeanSquare(HRfct, HRdV, r1_HRdV)
LRfct_rms = dgms.RootMeanSquare(LRfct, LRdV, r1_LRdV)
LUfct_rms = dgms.RootMeanSquare(LUfct, LRdV, r1_LRdV)

# %%
title_string = r"$\mathrm{RMS}(\zeta^{\prime}) = " + \
    "\int_{V}(\zeta^{\prime}(x,y,z,t))^{2}\,\mathrm{d}x\mathrm{d}y\mathrm{d}z$"

fig = plt.figure()
plt.plot(HRfct_rms, label = r"R27d")
plt.plot(LRfct_rms, label = r"R3d")
plt.plot(LUfct_rms, label = r"R3LU")
plt.title(title_string, fontsize=15.5)
plt.legend()

out_name = IO_dict["out_dir"] + "/" + os.path.basename(__file__)[:-3] +"_fctRMS"
plt.savefig(out_name + ".pdf", dpi=fig.dpi, bbox_inches='tight')
plt.savefig(out_name + ".eps", dpi=fig.dpi, bbox_inches='tight')
plt.show()

# %% Compute and print RMS(avg)
HRfct_rms_z = dgms.RootMeanSquare(HRfct, HRdA, r1_HRdA)
LRfct_rms_z = dgms.RootMeanSquare(LRfct, LRdA, r1_LRdA)
LUfct_rms_z = dgms.RootMeanSquare(LUfct, LRdA, r1_LRdA)
# %%
var2plt = np.stack((HRfct_rms_z[:600, :-1], 
                     LRfct_rms_z[:600, :-1], 
                     LUfct_rms_z[:600, :-1]), axis=0)
# %%
fig, axes_list = plt.subplots(np.size(var2plt, axis=0), 1, figsize=(5,2.5))
fig.tight_layout()
titles = [r"$\mathrm{R27d}$",
          r"$\mathrm{R3d}$",
          r"$\mathrm{R3LU}$"]
title_string = r"$\mathrm{RMS}(\zeta^{\prime})(z) = " + \
    "\int_{A}(\zeta^{\prime}(x,y,z,t))^{2}\,\mathrm{d}x\mathrm{d}y$"

plt.suptitle(title_string)
for idx, ax in enumerate(axes_list):
    # Plot command
    ax.imshow(np.transpose(var2plt[idx, :, :]), cmap='RdBu_r', interpolation='none') 

    ax.tick_params(which='major', top=False, bottom=False, left=False, right=False)
    ax.tick_params(which='minor', top=False, bottom=False, left=False, right=False)
    ax.xaxis.grid(False, which='major')
    ax.xaxis.grid(False, which='minor')
    
    # Y-axis parameters
    ax.set_ylabel(titles[idx] + "\n" + "depth [m]", fontsize = 5.5)
    ax.set_xticks([])
    ax.set_yticks([])
    
out_name = IO_dict["out_dir"] + "/" + os.path.basename(__file__)[:-3] +"_fctRMS_z"
plt.savefig(out_name + ".pdf", dpi=fig.dpi, bbox_inches='tight')
plt.savefig(out_name + ".eps", dpi=fig.dpi, bbox_inches='tight')
plt.show()

# %%
LRavg_rmse = dgms.RootMeanSquareERROR(LRavg, HRavg_ds, LRdV, r1_LRdV)
LUavg_rmse = dgms.RootMeanSquareERROR(LUavg, HRavg_ds, LRdV, r1_LRdV)
LRstd_rmse = dgms.RootMeanSquareERROR(LRstd, HRstd_ds, LRdV, r1_LRdV)
LUstd_rmse = dgms.RootMeanSquareERROR(LUstd, HRstd_ds, LRdV, r1_LRdV)

print("RMSE(avg),                  LR,                  LU")
print("         ", LRavg_rmse, LUavg_rmse)
print("RMSE(std),                  LR,                  LU")
print("         ", LRstd_rmse, LUstd_rmse)

# %
def MeanflowPCC(estimator, observable, weights):
    
    r1_weight_sum = 1. / np.sum(weights, axis=(-3,-2,-1))
    titles[0]
    mu_est = np.sum(np.multiply(estimator, weights), axis=(-3, -2, -1))
    mu_est = mu_est * r1_weight_sum
    
    mu_obs = np.sum(np.multiply(observable, weights), axis=(-3, -2, -1))   
    mu_obs = mu_obs * r1_weight_sum
    
    est_fct = estimator - mu_est
    obs_fct = observable - mu_obs
    
    temp = np.multiply(est_fct,obs_fct)
    
    temp3 = np.sum(np.multiply(est_fct**2, weights), axis=(-3, -2, -1)) * r1_weight_sum
    temp4 = np.sum(np.multiply(obs_fct**2, weights), axis=(-3, -2, -1)) * r1_weight_sum
    
    num = np.sum(np.multiply(temp, weights), axis=(-3, -2, -1)) * r1_weight_sum
    den = np.sqrt(temp3 * temp4)
    

    mPC = num/den

    return mPC


LRpc = MeanflowPCC(LRavg, HRavg_ds, LRdV)
LUpc = MeanflowPCC(LUavg, HRavg_ds, LRdV)
print("mPC,                  LR,                  LU")
print("     ", LRpc, LUpc)


# %% Gaussian Relative Entropy
# -------------------------

HRavg_norm = np.sum(np.multiply(HRavg_ds, LRdV), axis=(-3,-2,-1))
LRavg_norm = np.sum(np.multiply(LRavg, LRdV), axis=(-3,-2,-1))
LUavg_norm = np.sum(np.multiply(LUavg, LRdV), axis=(-3,-2,-1))

HRvnc_norm = np.sum(np.multiply(HRvnc_ds, LRdV), axis=(-3,-2,-1))
LRvnc_norm = np.sum(np.multiply(LRvnc, LRdV), axis=(-3,-2,-1))
LUvnc_norm = np.sum(np.multiply(LUvnc, LRdV), axis=(-3,-2,-1))

LRgre_int = dgms.AveragedGaussianRelativeEntropy(LRavg, HRavg_ds, LRvnc, HRvnc_ds, LRdV)
LUgre_int = dgms.AveragedGaussianRelativeEntropy(LUavg, HRavg_ds, LUvnc, HRvnc_ds, LRdV)

print("GRE,                  LR,                  LU")
print("     ", LRgre_int, LUgre_int)
# Plotting and plotting settings
# ------------------------------
titles = [r"$\mathrm{RMS}(\zeta_{_{\mathrm{R27d}}}) = " + "{:.2e}".format(HRavg_rms) + "$",
          r"$\mathrm{RMS}(\zeta_{_{\mathrm{R3d}}}) = " + "{:.2e}".format(LRavg_rms) + "$",
          r"$\mathrm{RMS}(\zeta_{_{\mathrm{R3LU}}}) = " +"{:.2e}".format(LUavg_rms) + "$"]
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Helvetica"],
    "figure.figsize" : (12/2.54, 8/2.54),
    "figure.dpi": 300
    })

nx = np.zeros(3)
ny = np.zeros(3)
rx = np.zeros(3)
ry = np.zeros(3)

nx[1:] = LRnx
ny[1:] = LRny
nx[0] = HRnx
ny[0] = HRny

rx[1:] = 1
ry[1:] = 1
rx[0] = (nx[0] - 2)/(nx[1] - 2)
ry[0] = (ny[0] - 2)/(nx[1] - 2)


maxval = 0.000004
minval = -maxval

fig, (ax1, ax2, ax3) = plt.subplots(1, 3)


im0 = ax1.imshow(np.rot90(HRavg[idk, ::-1, :]), cmap='RdBu_r', vmin = minval, vmax = maxval)
im1 = ax2.imshow(np.rot90(LRavg[idk, ::-1, :]), cmap='RdBu_r', vmin = minval, vmax = maxval)
im2 = ax3.imshow(np.rot90(LUavg[idk, ::-1, :]), cmap='RdBu_r', vmin = minval, vmax = maxval)

ax1.set_title(titles[0], fontsize=4.5)
ax2.set_title(titles[1], fontsize=4.5)
ax3.set_title(titles[2], fontsize=4.5)


box1 = ax1.get_position()
box2 = ax2.get_position()
box3 = ax3.get_position()

axes_list = [ax1, ax2, ax3]
for idx, ax in enumerate(axes_list):
    ax.tick_params(which='major', label1On=False)
    ax.tick_params(which='major', axis='both', width=0.25,  direction="in")
    ax.tick_params(which='minor', axis='both', width=0.125, direction='in')

    ax.tick_params(which='major', top=True, bottom=True, left=True, right=True)
    ax.tick_params(which='minor', top=True, bottom=True, left=True, right=True)

    # Major ticks
    ax.set_xticks(np.arange(ry[idx]*10, ny[idx]-1, ry[idx]*10))
    ax.set_yticks(np.arange(rx[idx]*10, nx[idx]-1, rx[idx]*10))

    # Minor ticks
    ax.set_xticks(np.arange(5, ny[idx]-1, 5), minor=True)
    ax.set_yticks(np.arange(5, nx[idx]-1, 5), minor=True)

    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(0.5)

#fig.text(1.75, 0.5,
#         'Preliminary version', transform=ax1.transAxes,
#         fontsize=30, color='gray', alpha=0.75,
#         ha='center', va='center', rotation='20')

# cbar_ax = fig.add_axes([box1.p0[0]-0.05, box1.p0[1], 0.02, box1.height])
# ticks = np.linspace(-maxval, maxval, 5, endpoint=True)
# cb = fig.colorbar(im2, cax=cbar_ax, ticks=ticks)
# cb.ax.tick_params(labelsize=5, direction='in', width=0.25,
#                   labelleft=False, labelright=True)

# cb.ax.set_yticklabels([r'-0.2', r'-0.1', r'0.0', r'0.1', r'0.2'],
#                       rotation=90, va="center")

# fig.text(0.5, 0.5, "Preliminary", transform=ax.transAxes,
#         fontsize=20, color='gray', alpha=0.5,
#         ha='center', va='center', rotation=45)

out_name = IO_dict["out_dir"] + "/" + os.path.basename(__file__)[:-3] +"_RMS"
plt.savefig(out_name + ".pdf", dpi=fig.dpi, bbox_inches='tight')
plt.savefig(out_name + ".eps", dpi=fig.dpi, bbox_inches='tight')
plt.show()

# %%

file = IO_dict["out_dir"] + "/" + os.path.basename(__file__)[:-3] +"_RMS_avg.txt"
np.savetxt(file, np.array([HRavg_rms, LRavg_rms, LUavg_rms]), header="HR, LR, LU")


file = IO_dict["out_dir"] + "/" + os.path.basename(__file__)[:-3] +"_RMS_fct.txt"
np.savetxt(file, np.array([HRfct_rms, LRfct_rms, LUfct_rms]), header="HR, LR, LU")


file = IO_dict["out_dir"] + "/" + os.path.basename(__file__)[:-3] +"_GRE.txt"
np.savetxt(file, np.array([LRgre_int, LUgre_int]), header="LR, LU")
