#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 15:14:29 2022

@author: ftucciar
"""
# Load modules
# ------------
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt

# Set matplotlib general parameters
# ---------------
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Helvetica"],
    "figure.dpi": 300
    })

# Set parameters
# ---------------
# Output directory
OUTdir = "/home/ftucciar/LU_postprocessing/netCDF_plot/outputs/"

# LR inputs
LRbase_dir = '/home/ftucciar/NEMO_downloads/'
LRsubs_dir = ['proms23/DET00',
              'Schmidt_1000/proms23/EXP00']
LRprnt_lab = [r"R3d",
              r"R3LU"]
LRbasefile = ['GYRE_5d_00010101_00101230_grid_',
              'GYRE_5d_00010101_00101230_grid_']

# HR inputs
HRbase_dir = '/home/ftucciar/dataNEMO/Deter/R27d/'
HRsubs_dir = ['100-102y', '102-104y', '104-106y', '106-108y', '108-110y']
HRbasefile = 'GYRE_5d_00010101_00021230_grid_'
HRprnt_lab = [r"R27d"]
grid2plt = 'T'
ext = '.nc'

LRgrid = '/home/ftucciar/dataNEMO/Deter/domain_cfg_R3.nc'
LRgrid = Dataset(LRgrid, "r", format="NETCDF4")
HRgrid = '/home/ftucciar/dataNEMO/Deter/domain_cfg_R27.nc'
HRgrid = Dataset(HRgrid, "r", format="NETCDF4")

HRvarname = 'votemper'
LRvarname = 'votemper'


# Read corresponding to the Low resolution models
# -----------------------------------------------
DTdata = Dataset(LRbase_dir + LRsubs_dir[0] + "/" + LRbasefile[0] +
                 grid2plt + ext, "r", format="NETCDF4")
LRdata = Dataset(LRbase_dir + LRsubs_dir[1] + "/" + LRbasefile[1] +
                 grid2plt + ext, "r", format="NETCDF4")

# Coordinates for grid-matching
LRlon = LRgrid.variables['nav_lon'][:, :].data
LRlat = LRgrid.variables['nav_lat'][:, :].data
# U-points (for integration) --------------------------------------------
LRe1u = LRgrid.variables['e1u'][0, :, :].data
LRe2u = LRgrid.variables['e2u'][0, :, :].data
# V-points (for integration) --------------------------------------------
LRe1v = LRgrid.variables['e1v'][0, :, :].data
LRe2v = LRgrid.variables['e2v'][0, :, :].data
# W-points (for integration) -------------------------720-------------------
LRe1t = LRgrid.variables['e1t'][0, :, :].data
LRe2t = LRgrid.variables['e2t'][0, :, :].data
LRe3t = LRgrid.variables['e3t_0'][0, :, :].data
LRe3t_1d = LRgrid.variables['e3t_1d'][0, :].data
depth = np.cumsum(np.append([5], LRe3t_1d[:-1]))
# f-points (for integration) --------------------------------------------
LRe1f = LRgrid.variables['e1f'][0, :, :].data
LRe2f = LRgrid.variables['e2f'][0, :, :].data
v3 = np.multiply(LRe1f, LRe2f)
# Mid-point rule area at u and v points ---------------------------------
LRdAu = np.multiply(LRe1u, LRe2u)
LRdAv = np.multiply(LRe1v, LRe2v)
LRdAt = np.multiply(LRe1t, LRe2t)
LRdVt = np.multiply(LRdAt, LRe3t)
# Planetary vorticity
LR_fft = LRgrid.variables['ff_t'][0, :, :].data
# Free memory
del LRe1u, LRe2u, LRe1v, LRe2v, LRe1t, LRe2t


# Read corresponding to the High resolution models
# ------------------------------------------------
HRdata = Dataset(HRbase_dir + HRsubs_dir[0] + "/" +  HRbasefile + grid2plt + ext, 
                      "r", format="NETCDF4")
HRvardims = HRdata.variables[HRvarname].dimensions
# # Coordinates for grid-matching
# HRlon = HRgrid.variables['nav_lon'][:, :].data
# HRlat = HRgrid.variables['nav_lat'][:, :].data
# # U-points (for integration) --------------------------------------------
# HRe1u = HRgrid.variables['e1u'][0, :, :].data
# HRe2u = HRgrid.variables['e2u'][0, :, :].data
# # V-points (for integration) --------------------------------------------
# HRe1v = HRgrid.variables['e1v'][0, :, :].data
# HRe2v = HRgrid.variables['e2v'][0, :, :].data
# # W-points (for integration) --------------------------------------------
HRe1t = HRgrid.variables['e1t'][0, :, :].data
HRe2t = HRgrid.variables['e2t'][0, :, :].data
HRe3t = HRgrid.variables['e3t_0'][0, :, :].data
# # f-points (for integration) --------------------------------------------
# HRe1f = HRgrid.variables['e1f'][0, :, :].data
# HRe2f = HRgrid.variables['e2f'][0, :, :].data
# v3 = np.multiply(LRe1f, LRe2f)
# # Mid-point rule area at u and v points ---------------------------------
# HRdAu = np.multiply(HRe1u, HRe2u)
# HRdAv = np.multiply(HRe1v, HRe2v)
HRdAt = np.multiply(HRe1t, HRe2t)
HRdVt = np.multiply(HRdAt, HRe3t)
# # Planetary vorticity
# HR_fft = HRgrid.variables['ff_t'][0, :, :].data
# # Free memory
# del HRe1u, HRe2u, HRe1v, HRe2v, HRe1t, HRe2t

HRzT = np.empty((5, 30))
for ids, s in enumerate(HRsubs_dir):
    inputfile = HRbase_dir + s + "/" +  HRbasefile + grid2plt + ext
    print(ids, s, inputfile)
    if np.size(HRvardims) == 3:
        fin = Dataset(inputfile, "r", format="NETCDF4")
        HRvar = fin.variables[HRvarname][-1, :, :]
    elif np.size(HRvardims) == 4:
        fin = Dataset(inputfile, "r", format="NETCDF4")
        HRvar = fin.variables[HRvarname][-1, :-1, :, :]
    fin.close
    HRzT[ids, :] = np.sum(np.multiply(HRvar, HRdAt), axis=( -2, -1)
                    )/np.sum(HRdAt, axis=( -2, -1))


# HRvardims = HRdata.variables[HRvarname].dimensions
LRvardims = LRdata.variables[LRvarname].dimensions

# Get domain dimension from domain file
# -------------------------------------
nx = np.zeros(4)
ny = np.zeros(4)
nz = np.zeros(4)
nt = np.zeros(4)

nx[1:] = LRdata.dimensions['x'].size
ny[1:] = LRdata.dimensions['y'].size
nz[1:] = LRdata.dimensions['depth'+grid2plt.lower()].size
nt[1:] = LRdata.dimensions['time_counter'].size

# nx[0] = HRdata.dimensions['x'].size
# ny[0] = HRdata.dimensions['y'].size
# nz[0] = HRdata.dimensions['depthu'].size

s = HRsubs_dir[0]

# nx[0] = HRdata.dimensions[HRvardims[-1]].size
# ny[0] = HRdata.dimensions[HRvardims[-2]].size
nx[1:] = LRdata.dimensions[LRvardims[-1]].size
ny[1:] = LRdata.dimensions[LRvardims[-2]].size
nt[1:] = LRdata.dimensions[LRvardims[0]].size

# if np.size(HRvardims) == 4:
#     nz[0] = HRdata.dimensions[HRvardims[-3]].size
#     nz[1:] = LRdata.dimensions[LRvardims[1]].size

# %% Reading LR temperature
times = [144, 288, 432, 576, 719]
DTzT = np.empty((5, 30))
LRzT = np.empty((5, 30))
for it, t in enumerate(times):
    print(t)
    if np.size(LRvardims) == 3:
        DTvar = DTdata.variables['votemper'][t, :, :]
        LRvar = LRdata.variables['votemper'][t, :, :]
    elif np.size(LRvardims) == 4:
        DTvar = DTdata.variables['votemper'][t, :-1, :, :]
        LRvar = LRdata.variables['votemper'][t, :-1, :, :]
    DTzT[it, :] = np.sum(np.multiply(DTvar, LRdAt), axis=( -2, -1)
                        )/np.sum(LRdAt, axis=( -2, -1))
    LRzT[it, :] = np.sum(np.multiply(LRvar, LRdAt), axis=( -2, -1)
                        )/np.sum(LRdAt, axis=( -2, -1)) 


# %% Plotting routine


# Userdef parameters
sup_title = r"Vertical profile of temperature" 

AxisY = depth[:-1]

positiveY = False
positiveX = True

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Helvetica"],
    "figure.figsize" : (12/2.54, 8/2.54),
    "figure.dpi": 300
    })

height = 1
fig,  (ax1, ax2, ax3, ax4, ax5) = plt.subplots(1,5)

# Plot title
plt.suptitle(sup_title, fontsize=10)

axes_list = [ax1, ax2, ax3, ax4, ax5]
for idx, ax in enumerate(axes_list):
    # Plot command
    ax.plot(HRzT[idx, :], AxisY, color="seagreen")
    ax.plot(DTzT[idx, :], AxisY, color="orange")
    ax.plot(LRzT[idx, :], AxisY, color="maroon")
    
    ax.xaxis.grid(False, which='major')
    ax.xaxis.grid(False, which='minor')

    
    # Y-axis parameters
    # ax.set_ylabel(r"depth [m]")
    
    if (positiveY): 
        ax.set_ylim(AxisY[0], AxisY[-1])
    else: 
        ax.set_ylim(AxisY[-1], AxisY[0])
    
    ax.set_yticks(AxisY, minor = True)
    ax.yaxis.grid(False, which='major')
    ax.yaxis.grid(True, which='minor')
    ax.tick_params(axis='y', which='major', length=7)
    ax.tick_params(axis='y', which='minor', length=0, color='r')
    ax.set_yticks([])
    
    ax.set_title(r"t=" + str((idx+1)*2) + "Y", fontsize=9)
    ax.set_xlim([3, 18])
    ax.set_xlabel(r"[$\null^{\circ}\mathrm{C}$]", fontsize=9)
    
    ax.tick_params(axis='x', which='major', length=2, labelsize=9)


plt.savefig(OUTdir + "Temperature_vertProfile" + ".eps",
            format='eps', dpi='figure', bbox_inches='tight')
plt.show()
