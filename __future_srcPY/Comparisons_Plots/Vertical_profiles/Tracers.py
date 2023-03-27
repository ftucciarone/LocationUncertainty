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
LRbase_dir = '/home/ftucciar/dataNEMO/'
LRsubs_dir = ['debug/Deterministic1y',
              'debug/0w1_36w2',
              'debug/36w1_72w2', 
              'debug/36w1_144w2',
              'debug/36w1_18w2', 
              'debug/72w1_18w2']
LRprnt_lab = [r"Deterministic", 
              r"$(                I - \mathcal{F}_{36} )\boldsymbol{u}$",
              r"$( \mathcal{F}_{36} - \mathcal{F}_{72} )\boldsymbol{u}$", 
              r"$( \mathcal{F}_{36} - \mathcal{F}_{144} )\boldsymbol{u}$",
              r"$( \mathcal{F}_{36} - \mathcal{F}_{18} )\boldsymbol{u}$", 
              r"$( \mathcal{F}_{72} - \mathcal{F}_{18} )\boldsymbol{u}$"]
LRbasefile = ['GYRE_5d_00010101_00051230_grid_',
              'GYRE_5d_00010101_00011230_grid_',
              'GYRE_5d_00010101_00011230_grid_',
              'GYRE_5d_00010101_00011230_grid_',
              'GYRE_5d_00010101_00011230_grid_',
              'GYRE_5d_00010101_00011230_grid_']

# HR inputs
HRbase_dir = '/home/ftucciar/dataNEMO/Deter/R27d/'
HRsubs_dir = ['100-102y', '102-104y', '104-106y', '106-108y', '108-110y']
HRbasefile = 'GYRE_5d_00010101_00021230_grid_'

grid2plt = 'T'
ext = '.nc'

LRgrid = '/home/ftucciar/dataNEMO/Deter/domain_cfg_R3.nc'
LRgrid = Dataset(LRgrid, "r", format="NETCDF4")
HRgrid = '/home/ftucciar/dataNEMO/Deter/domain_cfg_R27.nc'
HRgrid = Dataset(HRgrid, "r", format="NETCDF4")

HRvarname = 'votemper'
LRvarname = 'votemper'
idt = 70

# Read corresponding to the Low resolution models
# -----------------------------------------------
DTdata = Dataset(LRbase_dir + LRsubs_dir[0] + "/" + LRbasefile[0] +
                 grid2plt + ext, "r", format="NETCDF4")
LRdata = Dataset(LRbase_dir + LRsubs_dir[1] + "/" + LRbasefile[1] +
                 grid2plt + ext, "r", format="NETCDF4")
LU1data = Dataset(LRbase_dir + LRsubs_dir[2] + "/" + LRbasefile[2] +
                  grid2plt + ext, "r", format="NETCDF4")
LU2data = Dataset(LRbase_dir + LRsubs_dir[3] + "/" + LRbasefile[3] +
                  grid2plt + ext, "r", format="NETCDF4")
LU3data = Dataset(LRbase_dir + LRsubs_dir[4] + "/" + LRbasefile[4] +
                  grid2plt + ext, "r", format="NETCDF4")
LU4data = Dataset(LRbase_dir + LRsubs_dir[5] + "/" + LRbasefile[5] +
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


# # Read corresponding to the High resolution models
# # ------------------------------------------------
# HRdata = Dataset(HRbase_dir + HRsubs_dir[0] + "/" + HRbasefile +
#                  grid2plt + ext, "r", format="NETCDF4")
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
# HRe1t = HRgrid.variables['e1t'][0, :, :].data
# HRe2t = HRgrid.variables['e2t'][0, :, :].data
# HRe3t = HRgrid.variables['e3t_0'][0, :, :].data
# # f-points (for integration) --------------------------------------------
# HRe1f = HRgrid.variables['e1f'][0, :, :].data
# HRe2f = HRgrid.variables['e2f'][0, :, :].data
# v3 = np.multiply(LRe1f, LRe2f)
# # Mid-point rule area at u and v points ---------------------------------
# HRdAu = np.multiply(HRe1u, HRe2u)
# HRdAv = np.multiply(HRe1v, HRe2v)
# HRdAt = np.multiply(HRe1t, HRe2t)
# HRdVt = np.multiply(HRdAt, HRe3t)
# # Planetary vorticity
# HR_fft = HRgrid.variables['ff_t'][0, :, :].data
# # Free memory
# del HRe1u, HRe2u, HRe1v, HRe2v, HRe1t, HRe2t

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
if np.size(LRvardims) == 3:
    DTvar = DTdata.variables['votemper'][idt, :, :]
    LRvar = LRdata.variables['votemper'][idt, :, :]
    LU1var = LU1data.variables['votemper'][idt, :, :]
    LU2var = LU2data.variables['votemper'][idt, :, :]
    LU3var = LU3data.variables['votemper'][idt, :, :]
    LU4var = LU4data.variables['votemper'][idt, :, :]
elif np.size(LRvardims) == 4:
    DTvar = DTdata.variables['votemper'][idt, :, :, :]
    LRvar = LRdata.variables['votemper'][idt, :, :, :]
    LU1var = LU1data.variables['votemper'][idt, :, :, :]
    LU2var = LU2data.variables['votemper'][idt, :, :, :]
    LU3var = LU3data.variables['votemper'][idt, :, :, :]
    LU4var = LU4data.variables['votemper'][idt, :, :, :]

DTzT = np.sum(np.multiply(DTvar, LRdAt), axis=( -2, -1)
               )/np.sum(LRdAt, axis=( -2, -1))
LRzT = np.sum(np.multiply(LRvar, LRdAt), axis=( -2, -1)
               )/np.sum(LRdAt, axis=( -2, -1))
LU1zT = np.sum(np.multiply(LU1var, LRdAt), axis=( -2, -1)
               )/np.sum(LRdAt, axis=( -2, -1))
LU2zT = np.sum(np.multiply(LU2var, LRdAt), axis=( -2, -1)
                )/np.sum(LRdAt, axis=( -2, -1))
LU3zT = np.sum(np.multiply(LU3var, LRdAt), axis=( -2, -1)
               )/np.sum(LRdAt, axis=( -2, -1))
LU4zT = np.sum(np.multiply(LU4var, LRdAt), axis=( -2, -1)
                )/np.sum(LRdAt, axis=( -2, -1))

# %% Reading LR salinity
if np.size(LRvardims) == 3:
    DTvar = DTdata.variables['vosaline'][idt, :, :]
    LRvar = LRdata.variables['vosaline'][idt, :, :]
    LU1var = LU1data.variables['vosaline'][idt, :, :]
    LU2var = LU2data.variables['vosaline'][idt, :, :]
    LU3var = LU3data.variables['vosaline'][idt, :, :]
    LU4var = LU4data.variables['vosaline'][idt, :, :]
elif np.size(LRvardims) == 4:
    DTvar = DTdata.variables['vosaline'][idt, :, :, :]
    LRvar = LRdata.variables['vosaline'][idt, :, :, :]
    LU1var = LU1data.variables['vosaline'][idt, :, :, :]
    LU2var = LU2data.variables['vosaline'][idt, :, :, :]
    LU3var = LU3data.variables['vosaline'][idt, :, :, :]
    LU4var = LU4data.variables['vosaline'][idt, :, :, :]

DTzS = np.sum(np.multiply(DTvar, LRdAt), axis=( -2, -1)
               )/np.sum(LRdAt, axis=( -2, -1))
LRzS = np.sum(np.multiply(LRvar, LRdAt), axis=( -2, -1)
               )/np.sum(LRdAt, axis=( -2, -1))
LU1zS = np.sum(np.multiply(LU1var, LRdAt), axis=( -2, -1)
               )/np.sum(LRdAt, axis=( -2, -1))
LU2zS = np.sum(np.multiply(LU2var, LRdAt), axis=( -2, -1)
               )/np.sum(LRdAt, axis=( -2, -1))
LU3zS = np.sum(np.multiply(LU3var, LRdAt), axis=( -2, -1)
               )/np.sum(LRdAt, axis=( -2, -1))
LU4zS = np.sum(np.multiply(LU4var, LRdAt), axis=( -2, -1)
               )/np.sum(LRdAt, axis=( -2, -1))
# %% Plotting routine

AxisX = np.empty((2, 6, 30))
# Userdef parameters
sup_title = r"Vertical profile of tracers at instant t=" + str(idt)
sub_title = [r"Temperature", r"Salinity"]

AxisY = depth[:-1]
AxisX[0, :, :] = np.array([DTzT[:-1], LRzT[:-1], LU1zT[:-1], LU2zT[:-1], LU3zT[:-1], LU4zT[:-1]])
AxisX[1, :, :] = np.array([DTzS[:-1], LRzS[:-1], LU1zS[:-1], LU2zS[:-1], LU3zS[:-1], LU4zS[:-1]] )

positiveY = False
positiveX = True

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Helvetica"],
    "figure.dpi": 300
    })

height = 2
fig,  (ax1, ax2) = plt.subplots(1, 2)

# Plot title
plt.suptitle(sup_title)

axes_list = [ax1, ax2]
for idx, ax in enumerate(axes_list):
    # Plot command
    ax.plot(AxisX[idx, 0, :], AxisY, '--', label=LRprnt_lab[0])
    ax.plot(AxisX[idx, 1, :], AxisY, '--', label=LRprnt_lab[1])
    ax.plot(AxisX[idx, 2, :], AxisY, '--', label=LRprnt_lab[2])
    ax.plot(AxisX[idx, 3, :], AxisY, '--', label=LRprnt_lab[3])
    ax.plot(AxisX[idx, 4, :], AxisY, '--', label=LRprnt_lab[4])
    ax.plot(AxisX[idx, 5, :], AxisY, '--', label=LRprnt_lab[5])
    
    # X-axis parameters
    ax.set_xlabel(sub_title[idx])
    
    if (positiveX): 
        ax.set_xlim(0.9*np.min(AxisX[idx, :, :]), 1.1*np.max(AxisX[idx, :, :]))
    else:
        ax.set_xlim(1.1*np.min(AxisX[idx, :, :]), 0.9*np.max(AxisX[idx, :, :]))
    
    ax.xaxis.grid(False, which='major')
    ax.xaxis.grid(False, which='minor')
    ax.tick_params(axis='y', which='major', length=7)
    ax.tick_params(axis='y', which='minor', length=0, color='r')
    
    # Y-axis parameters
    ax.set_ylabel(r"depth [m]")
    
    if (positiveY): 
        ax.set_ylim(AxisY[0], AxisY[-1])
    else: 
        ax.set_ylim(AxisY[-1], AxisY[0])
    
    ax.set_yticks(AxisY, minor = True)
    ax.yaxis.grid(False, which='major')
    ax.yaxis.grid(True, which='minor')
    ax.tick_params(axis='y', which='major', length=7)
    ax.tick_params(axis='y', which='minor', length=0, color='r')

ax1.set_xlim([0, 30])
ax2.set_xlim([25, 35])
ax1.legend(loc="lower right")
ax2.legend(loc="lower left")
ax2.set_ylabel("")
ax2.set_yticks([])

# plt.savefig(OUTdir + "Tracers_VertProfile_t" + str(idt) + ".eps",
#             format='eps', dpi='figure', bbox_inches='tight')
plt.show()
