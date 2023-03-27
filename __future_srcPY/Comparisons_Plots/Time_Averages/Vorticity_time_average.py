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
data_dir = home_dir + "/dataNEMO"
dom_dir = home_dir + "/dataNEMO/Deter"
out_dir = work_dir + "/netCDFplot/Outputs"

sub_dir = ["/Deter/R27d",
           "/../NEMO_downloads/proms23/DET00",
           "/../NEMO_downloads/Schmidt_1000/proms23/EXP00/"]

HRsubs_dir = ['/100-102y', '/102-104y', '/104-106y', '/106-108y', '/108-110y']

filename = ['curl_',
            'GYRE_5d_00010101_00101230_grid_',
            'GYRE_5d_00010101_00101230_grid_']
# Print labels
LRprnt_lab = [r"R27d", 
              r"R3d",
              r"R3LU"]

grid2plt = 'T'
ext = '.nc'

HRvarname = 'socurlt'
LRvarname = 'vorticity'

print("")
print("Sanity check on existence of folders and files")
print("")
print(" isdir=", os.path.isdir(home_dir), "home_dir = ", home_dir)
print(" isdir=", os.path.isdir(work_dir), "work_dir = ", work_dir)
print(" isdir=", os.path.isdir(data_dir), "data_dir = ", data_dir)
print(" isdir=", os.path.isdir(dom_dir), " dom_dir = ", dom_dir)
print(" isdir=", os.path.isdir(out_dir), " out_dir = ", out_dir)
print("")
for i in range(len(sub_dir)):
    print(" isdir=", os.path.isdir(data_dir + sub_dir[i]), "data_subdir = ", data_dir + sub_dir[i] )
print("")
for i in range(len(HRsubs_dir)):
    print(" isfile=", os.path.isfile(data_dir + sub_dir[0] + HRsubs_dir[i] + "/" +  filename[0] + grid2plt + ext ), \
          LRprnt_lab[0] + " = ", data_dir + sub_dir[0] + HRsubs_dir[i] + "/" + filename[0] + grid2plt + ext  )
print("")
for i in range(len(sub_dir)-1):
    print(" isfile=", os.path.isfile(data_dir + sub_dir[i+1] + "/" +  filename[i+1] + grid2plt + ext ), \
          LRprnt_lab[i+1] + " = ", data_dir + sub_dir[i+1] + "/" + filename[i+1] + grid2plt + ext  )
print("")

LRgrid = Dataset(dom_dir + "/domain_cfg_R3.nc", "r", format="NETCDF4")
HRgrid = Dataset(dom_dir + "/domain_cfg_R27.nc", "r", format="NETCDF4")

HRff_f = HRgrid.variables["ff_f"][0, ::-1, :]
HRny = HRgrid.dimensions["y"].size
HRnx = HRgrid.dimensions["x"].size

LRff_f = LRgrid.variables["ff_f"][0, ::-1, :]
LRny = LRgrid.dimensions["y"].size
LRnx = LRgrid.dimensions["x"].size

maxTime = 720
idk = 1
# Read corresponding to the High resolution models
# ------------------------------------------------
HRdata = Dataset(data_dir + sub_dir[0] + HRsubs_dir[0] + "/" +  filename[0] + grid2plt + ext, 
                      "r", format="NETCDF4")
HRvardims = HRdata.variables[HRvarname].dimensions
HRcmvar = np.zeros([542, 812])
nt = 0
for s in HRsubs_dir:
    print(s)
    if np.size(HRvardims) == 3:
        fin = Dataset(data_dir + sub_dir[0] + HRsubs_dir[i] + "/" +  filename[0] + grid2plt + ext, 
                      "r", format="NETCDF4")
        HRvar = fin.variables[HRvarname][:, ::-1, :]
        nt += fin.dimensions[HRvardims[0]].size
        HRcmvar += np.sum(HRvar, axis=0)
        del HRvar
    elif np.size(HRvardims) == 4:
        fin = Dataset(data_dir + sub_dir[0] + HRsubs_dir[i] + "/" +  filename[0] + grid2plt + ext,
                      "r", format="NETCDF4")
        HRvar = fin.variables[HRvarname][:, idk, ::-1, :]
        HRcmvar[:, :] += np.sum(HRvar, axis=0)
        del HRvar
        nt += fin.dimensions[HRvardims[0]].size
    fin.close
HRvar = np.divide(HRcmvar,HRff_f) / nt 

# Read corresponding to the Low resolution models
# -----------------------------------------------
LRdata = Dataset(data_dir + sub_dir[1] + "/" +  filename[1] + grid2plt + ext,
                 "r", format="NETCDF4")
LUdata = Dataset(data_dir + sub_dir[2] + "/" +  filename[2] + grid2plt + ext,
                 "r", format="NETCDF4")
# Get domain dimension from domain file
# -------------------------------------
LRvardims = LRdata.variables[LRvarname].dimensions

if np.size(LRvardims) == 3:
    LRvar = np.mean(LRdata.variables[LRvarname][1:maxTime, :, :], axis=0)
    LUvar = np.mean(LUdata.variables[LRvarname][1:maxTime, :, :], axis=0)
elif np.size(LRvardims) == 4:
    LRvar = np.mean(LRdata.variables[LRvarname][1:maxTime, idk, ::-1, :], axis=0)
    LUvar = np.mean(LUdata.variables[LRvarname][1:maxTime, idk, ::-1, :], axis=0)
LRvar = np.divide(LRvar,LRff_f) 
LUvar = np.divide(LUvar,LRff_f) 
# %%
maxval = 0.1

# Plotting and plotting settings
# ------------------------------
titles = [r"$\mathrm{R27d}$",
          r"$\mathrm{R3d}$",
          r"$\mathrm{R3LU}$"]
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


fig, (ax1, ax2, ax3) = plt.subplots(1, 3)

im0 = ax1.imshow(np.rot90(HRvar), cmap='RdBu', vmin = -maxval, vmax = maxval)
im1 = ax2.imshow(np.rot90(LRvar), cmap='RdBu', vmin = -maxval, vmax = maxval)
im2 = ax3.imshow(np.rot90(LUvar), cmap='RdBu', vmin = -maxval, vmax = maxval)

ax1.set_title(titles[0], fontsize=10)
ax2.set_title(titles[1], fontsize=10)
ax3.set_title(titles[2], fontsize=10)


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

cb.ax.set_yticklabels([r'-0.1', r'-0.05', r'0.0', r'0.05', r'0.1'],
                      rotation=90, va="center", fontsize=5)

# fig.text(0.5, 0.5, "Preliminary", transform=ax.transAxes,
#         fontsize=20, color='gray', alpha=0.5,
#         ha='center', va='center', rotation=30)

plt.savefig(out_dir + "/" + os.path.basename(__file__)[:-3] + ".pdf", dpi=fig.dpi, bbox_inches='tight')
plt.savefig(out_dir + "/" + os.path.basename(__file__)[:-3] + ".eps", dpi=fig.dpi, bbox_inches='tight')
plt.show()

