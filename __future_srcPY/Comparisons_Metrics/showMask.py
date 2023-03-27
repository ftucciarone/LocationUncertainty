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
import DiagnosticMetrics as dgms

grid2plt = "T"
ext = ".nc" 
HRvarname = "socurlt"
LRvarname = "vorticity"

HRfiles, LRfiles, IO_dict = yaml2dict(r'paths.yaml')

LRgrid = Dataset(IO_dict["dom_dir"] + "/domain_cfg_R3.nc", "r", format="NETCDF4")
HRgrid = Dataset(IO_dict["dom_dir"] + "/domain_cfg_R27.nc", "r", format="NETCDF4")

HRlon = HRgrid.variables['nav_lon'][:, :].data
HRlat = HRgrid.variables['nav_lat'][:, :].data
HRe1 = HRgrid.variables['e1t'][0, :, :].data
HRe2 = HRgrid.variables['e2t'][0, :, :].data
HRe3 = HRgrid.variables['e3t_0'][0, :, :].data
HRdA = np.multiply(HRe1, HRe2); del HRe1, HRe2
HRdV = np.multiply(HRdA, HRe3); del HRe3
r1_HRdA = 1 / np.sum(HRdA, axis=(-2, -1))
r1_HRdV = 1 / np.sum(HRdV, axis=(-3, -2, -1))
HRny = HRgrid.dimensions["y"].size
HRnx = HRgrid.dimensions["x"].size

LRlon = LRgrid.variables['nav_lon'][:, :].data
LRlat = LRgrid.variables['nav_lat'][:, :].data
LRe1 = LRgrid.variables['e1t'][0, :, :].data
LRe2 = LRgrid.variables['e2t'][0, :, :].data
LRe3 = LRgrid.variables['e3t_0'][0, :, :].data
LRdA = np.multiply(LRe1, LRe2); del LRe1, LRe2
LRdV = np.multiply(LRdA, LRe3); del LRe3
r1_LRdA = 1 / np.sum(LRdA, axis=(-2, -1))
r1_LRdV = 1 / np.sum(LRdV, axis=(-3, -2, -1))
LRny = LRgrid.dimensions["y"].size
LRnx = LRgrid.dimensions["x"].size

# Find mask for overlap (not elegant but working)
tr1 = np.repeat(np.repeat(LRlon[1:-1, 1:-1],
                          repeats=9*np.ones(60).astype(int), axis=0),
                repeats=9*np.ones(90).astype(int), axis=1) - HRlon[1:-1, 1:-1]
tr2 = np.repeat(np.repeat(LRlat[1:-1, 1:-1],
                          repeats=9*np.ones(60).astype(int), axis=0),
                repeats=9*np.ones(90).astype(int), axis=1) - HRlat[1:-1, 1:-1]
mask = tr1 == tr2
cols = np.linspace(0, HRnx-1, HRnx)
rows = np.linspace(0, HRny-1, HRny)
cols = cols[1:-1]
rows = rows[1:-1]

cols = cols[np.sum(mask, axis=0) < 60]
rows = rows[np.sum(mask, axis=1) < 1]

# %%
minTime = 0
maxTime = 400
idk = 1
# Read corresponding to the High resolution models
# ------------------------------------------------
HRdata = Dataset(HRfiles[0]["curl_files"][0], "r", format="NETCDF4")
HRvardims = HRdata.variables[HRvarname].dimensions
if np.size(HRvardims) == 3:
    HRcmvar_1m = np.zeros([542, 812])  # first moment cumulated
    HRcmvar_2m = np.zeros([542, 812])  # second moment cumulated
elif np.size(HRvardims) == 4:
    HRcmvar_1m = np.zeros([31, 542, 812])  # first moment cumulated
    HRcmvar_2m = np.zeros([31, 542, 812])  # second moment cumulated


fin = Dataset(HRfiles[0]["curl_files"][0], "r", format="NETCDF4")
HRvar = fin.variables[HRvarname][0, :, :]






# %% Plotting and plotting settings
# ------------------------------



def createJetStreamMask(lat_in, lon_in):
    latct = 29
    lonct = -85
    latwd = 5
    lonwd = 20
    
    
    lat_out = lat_in
    lon_out = lon_in
    
    lon_out = np.where(lon_in > lonct - lonwd, lon_in, 361)
    lon_out = np.where(lon_out < lonct + lonwd, lon_out, 361)
    
    lat_out = np.where(lat_in > latct - latwd, lat_in, 91)
    lat_out = np.where(lat_out < latct + latwd, lat_out, 91)
    mask_out = np.where(lat_out != 91, np.ones_like(lat_out), 0)
    mask_out = np.multiply(np.where(lon_out != 361, np.ones_like(lat_out), 0), mask_out)
    return lat_out, lon_out, mask_out













HRlat_mask, HRlon_mask, mask = createJetStreamMask(HRlat, HRlon)
HR2show = np.multiply(HRvar, mask)
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Helvetica"],
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

maxval = 0.00002 #np.max([np.nanmax(HRstd), np.nanmax(LRstd), np.nanmax(LUstd)])
minval = -maxval

fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
fig.set_size_inches(4.8, 3)

im0 = ax1.imshow(np.rot90(HR2show[idk, ::-1, :]), cmap='RdBu', vmin = minval, vmax = maxval)
im1 = ax2.imshow(np.rot90(HRlat_mask[ ::-1, :]), cmap='RdBu')
im2 = ax3.imshow(np.rot90(HRlon_mask[ ::-1, :]), cmap='RdBu')


"""
ax1.set_title(titles[0], fontsize=5.5)
ax2.set_title(titles[1], fontsize=5.5)
ax3.set_title(titles[2], fontsize=5.5)


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

cbar_ax = fig.add_axes([box1.p0[0]-0.05, box1.p0[1], 0.02, box1.height])
ticks = np.linspace(-maxval, maxval, 5, endpoint=True)
cb = fig.colorbar(im2, cax=cbar_ax, ticks=ticks)
cb.ax.tick_params(labelsize=5, direction='in', width=0.25,
                  labelleft=False, labelright=True)

cb.ax.set_yticklabels([r'-0.2', r'-0.1', r'0.0', r'0.1', r'0.2'],
                      rotation=90, va="center")

fig.text(0.5, 0.5, "Preliminary", transform=ax.transAxes,
        fontsize=20, color='gray', alpha=0.5,
        ha='center', va='center', rotation=30)

#out_name = IO_dict["out_dir"] + "/" + os.path.basename(__file__)[:-3] +"_RMS"
#plt.savefig(out_name + ".pdf", dpi=fig.dpi, bbox_inches='tight')
#plt.savefig(out_name + ".eps", dpi=fig.dpi, bbox_inches='tight')
plt.show()
"""