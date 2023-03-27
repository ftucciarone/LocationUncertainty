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


# Set parameters
# ---------------
# Output directory
home_dir = str(pathlib.Path.home())
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
LRvarname = 'vorticity'

# Read corresponding to the Low resolution models
# -----------------------------------------------
LRdata = Dataset(LRbase_dir + LRsubs_dir[0] + "/" + LRbasefile[0] +
                 grid2plt + ext, "r", format="NETCDF4")
LU1data = Dataset(LRbase_dir + LRsubs_dir[1] + LRbasefile[1] +
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

maxTime = 720
idk = 1
# %% Reading LR variables 
if np.size(LRvardims) == 3:
    LRvar = np.mean(LRdata.variables[LRvarname][1:maxTime, :, :], axis=0)
    LU1var = np.mean(LU1data.variables[LRvarname][1:maxTime, :, :], axis=0)
    LU2var = np.mean(LU2data.variables[LRvarname][1:maxTime, :, :], axis=0)
elif np.size(LRvardims) == 4:
    LRvar = np.mean(LRdata.variables[LRvarname][1:maxTime, idk, ::-1, :], axis=0)
    LU1var = np.mean(LU1data.variables[LRvarname][1:maxTime, idk, ::-1, :], axis=0)
    LU2var = np.mean(LU2data.variables[LRvarname][1:maxTime, idk, ::-1, :], axis=0)


plt.rcParams.update({
	"text.usetex": True,
	"font.family": "serif",
	"font.serif": ["Helvetica"],
	"figure.dpi": 300
})


titles = ["No projection", "With projection"]


# %% Vertical plots
time = np.linspace(0, maxTime-1, maxTime-1)

fig, (ax1, ax2) = plt.subplots(1, 2)
fig.set_size_inches(4.8, 3)

im2 = ax1.imshow(np.rot90(LU1var), cmap='RdBu')
ax2.imshow(np.rot90(LU2var), cmap='RdBu')


ax1.set_title(titles[0], fontsize=7.5)
ax2.set_title(titles[1], fontsize=7.5)


box1 = ax1.get_position()
box2 = ax2.get_position()

axes_list = [ax1, ax2]
for idx, ax in enumerate(axes_list):
    ax.tick_params(which='major', label1On=False)
    ax.tick_params(which='major', axis='both', width=0.25,  direction="in")
    ax.tick_params(which='minor', axis='both', width=0.125, direction='in')

    ax.tick_params(which='major', top=True, bottom=True, left=True, right=True)
    ax.tick_params(which='minor', top=True, bottom=True, left=True, right=True)



    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(0.5)

plt.savefig(OUTdir + 'Comparison_P_noP.png',
            format='png', dpi='figure', bbox_inches='tight')
plt.show()
