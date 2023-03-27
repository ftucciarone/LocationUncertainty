#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 15:14:29 2022

@author: ftucciar
"""
# Load modules
# ------------
import pathlib
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt

# Set matplotlib general parameters
# ---------------
#plt.rcParams.update({
#    "text.usetex": True,
#    "font.family": "serif",
#    "font.serif": ["Helvetica"],
#    "figure.dpi": 300
#    })

#plt.rcParams['text.latex.preamble'] = [
#    r'\usepackage{amsmath}',
#    r'\usepackage{amssymb}']

home_dir = str(pathlib.Path.home())

# Set parameters
# ---------------
# Output directory
OUTdir =  home_dir + "/LocationUncertainty/netCDFplot/Outputs/"

# LR inputs
LRbase_dir = home_dir + "/dataNEMO/"
LRsubs_dir = ["Deter",
              "",
              ""]
LRbasefile = ['5*GYRE_5d_00010101_00051230_grid_',
              'pod_unpj_GYRE_5d_00010101_00101230_grid_',
              'pod_proj_GYRE_5d_00010101_00101230_grid_']
# Print labels
LRprnt_lab = [r"Deterministic", 
              r"$( \mathcal{F}_{36} - \mathcal{F}_{144} )\boldsymbol{u}$"]

# HR inputs
#HRbase_dir = '/home/ftucciar/dataNEMO/Deter/R27d/'
#HRsubs_dir = ['100-102y', '102-104y', '104-106y', '106-108y', '108-110y']
#HRprnt_lab = [r"Deterministic"]
#HRbasefile = 'GYRE_5d_00010101_00021230_grid_'

grid2plt = 'T'
ext = '.nc'

LRgrid = home_dir + "/LocationUncertainty/data/domain_cfg_R3.nc"
LRgrid = Dataset(LRgrid, "r", format="NETCDF4")
HRgrid = home_dir + "/LocationUncertainty/data/domain_cfg_R27.nc"
HRgrid = Dataset(HRgrid, "r", format="NETCDF4")

HRvarname = 'votemper'
LRvarname = 'votemper'

# Read corresponding to the Low resolution models
# -----------------------------------------------
LRdata = Dataset(LRbase_dir + LRsubs_dir[0] + "/" + LRbasefile[0] +
                 grid2plt + ext, "r", format="NETCDF4")
LU1data = Dataset(LRbase_dir + LRsubs_dir[1] + "/" + LRbasefile[1] +
                  grid2plt + ext, "r", format="NETCDF4")
LU2data = Dataset(LRbase_dir + LRsubs_dir[2] + "/" + LRbasefile[2] +
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
LRtA = np.sum(LRdAt, axis=( -2, -1))
LRtV = np.sum(LRdVt, axis=( -2, -1))
# Planetary vorticity
LR_fft = LRgrid.variables['ff_t'][0, :, :].data
# Free memory
del LRe1u, LRe2u, LRe1v, LRe2v, LRe1t, LRe2t

# Get domain dimension from domain file
# -------------------------------------
LRvardims = LRdata.variables[LRvarname].dimensions
nx = np.zeros(4)
ny = np.zeros(4)
nz = np.zeros(4)
nt = np.zeros(4)

nx[1:] = LRdata.dimensions['x'].size
ny[1:] = LRdata.dimensions['y'].size
nz[1:] = LRdata.dimensions['depth'+grid2plt.lower()].size
nt[1:] = LRdata.dimensions['time_counter'].size

nx[1:] = LRdata.dimensions[LRvardims[-1]].size
ny[1:] = LRdata.dimensions[LRvardims[-2]].size
nt[1:] = LRdata.dimensions[LRvardims[0]].size

# %% Reading LR salinity

if np.size(LRvardims) == 3:
    LRvar = LRdata.variables[LRvarname][:, :, :]
    LU1var = LU1data.variables[LRvarname][:, :, :]
elif np.size(LRvardims) == 4:
    LRvar = LRdata.variables[LRvarname][:, :-1, :, :]
    LU1var = LU1data.variables[LRvarname][:, :-1, :, :]
    LU2var = LU2data.variables[LRvarname][:, :-1, :, :]

LRzS = np.sum(np.multiply(LRvar, LRdAt), axis=( -2, -1))/LRtA
LU1zS = np.sum(np.multiply(LU1var, LRdAt), axis=( -2, -1))/LRtA
LU2zS = np.sum(np.multiply(LU2var, LRdAt), axis=( -2, -1))/LRtA


# %% Plotting routine

tmax = 720

Salinity = np.stack((#LRzS[:tmax, :], 
                     LU1zS[:tmax, :], 
                     LU2zS[:tmax, :]), axis=0)


DT_Salinity = np.stack((#LRzS[:tmax-1, :]-LRzS[1:tmax, :],
                     	LU1zS[:tmax-1, :]-LU1zS[1:tmax, :], 
                     	LU2zS[:tmax-1, :]-LU2zS[1:tmax, :]), axis=0)



# Userdef parameters
sup_title = r"Vertical profile of tracers at instant t=" + str(1)
sub_title = [r"Temperature", r"Salinity"]


positiveY = False
positiveX = True

var2plt = DT_Salinity


height = 2
fig,  axes_list = plt.subplots(np.size(var2plt, axis=0), 1)

# Plot title
plt.suptitle(sup_title)


for idx, ax in enumerate(axes_list):
    # Plot command
    ax.imshow(np.transpose(var2plt[idx, :, :]),
              cmap='coolwarm',
               # vmin=np.min(var2plt), vmax=np.max(var2plt),
              interpolation='none')  


plt.savefig(OUTdir + "Tracers_Hovemoller.png",
            format='png', dpi='figure', bbox_inches='tight')
plt.show()





