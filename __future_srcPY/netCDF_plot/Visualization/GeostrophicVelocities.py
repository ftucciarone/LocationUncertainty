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
OUTdir =  home_dir + "/__future_LocationUncertainty/netCDFplot/Outputs/"

# LR inputs
LRbase_dir = home_dir + "/dataNEMO/Deter/R27d/"
LRsubs_dir = ["108-110y"]
LRbasefile = ['GYRE_5d_00010101_00021230_grid_']
# Print labels
LRprnt_lab = [r"Deterministic", 
              r"$( \mathcal{F}_{36} - \mathcal{F}_{144} )\boldsymbol{u}$"]

# HR inputs
HRbase_dir = '/home/ftucciar/dataNEMO/Deter/R27d/'
HRsubs_dir = ["108-110y"]#'100-102y', '102-104y', '104-106y', '106-108y', '108-110y']
HRprnt_lab = [r"Deterministic"]
HRbasefile = 'GYRE_5d_00010101_00021230_grid_'

grid2plt = 'T'
ext = '.nc'

LRgrid = home_dir + "/dataNEMO/Deter/domain_cfg_R3.nc"
LRgrid = Dataset(LRgrid, "r", format="NETCDF4")
HRgrid = home_dir + "/dataNEMO/Deter/domain_cfg_R27.nc"
HRgrid = Dataset(HRgrid, "r", format="NETCDF4")

HRvarname = 'sossheig'
LRvarname = 'sossheig'


# Read corresponding to the Low resolution models
# -----------------------------------------------
HRdata = Dataset(HRbase_dir + HRsubs_dir[0] + "/" + HRbasefile +
                 "T" + ext, "r", format="NETCDF4")


# Coordinates for grid-matching
HRlon = HRgrid.variables['nav_lon'][:, :].data
HRlat = HRgrid.variables['nav_lat'][:, :].data
# U-points (for integration) --------------------------------------------
HRe1u = HRgrid.variables['e1u'][0, :, :].data
HRe2u = HRgrid.variables['e2u'][0, :, :].data
# V-points (for integration) --------------------------------------------
HRe1v = HRgrid.variables['e1v'][0, :, :].data
HRe2v = HRgrid.variables['e2v'][0, :, :].data
# W-points (for integration) -------------------------720-------------------
HRe1t = HRgrid.variables['e1t'][0, :, :].data
HRe2t = HRgrid.variables['e2t'][0, :, :].data
HRe3t = HRgrid.variables['e3t_0'][0, :, :].data
HRe3t_1d = HRgrid.variables['e3t_1d'][0, :].data
depth = np.cumsum(np.append([5], HRe3t_1d[:-1]))
# f-points (for integration) --------------------------------------------
HRe1f = HRgrid.variables['e1f'][0, :, :].data
HRe2f = HRgrid.variables['e2f'][0, :, :].data
v3 = np.multiply(HRe1f, HRe2f)
# Mid-point rule area at u and v points ---------------------------------
HRdAu = np.multiply(HRe1u, HRe2u)
HRdAv = np.multiply(HRe1v, HRe2v)
HRdAt = np.multiply(HRe1t, HRe2t)
HRdVt = np.multiply(HRdAt, HRe3t)
HRtA = np.sum(HRdAt, axis=( -2, -1))
HRtV = np.sum(HRdVt, axis=( -2, -1))
# Planetary vorticity
HR_fff = HRgrid.variables['ff_f'][0, :, :].data
HR_fft = HRgrid.variables['ff_t'][0, :, :].data

# Free memory
#del LRe1u, LRe2u, LRe1v, LRe2v, LRe1t, LRe2t

# Get domain dimension from domain file
# -------------------------------------
HRvardims = HRdata.variables[HRvarname].dimensions
nx = np.zeros(4)
ny = np.zeros(4)
nz = np.zeros(4)
nt = np.zeros(4)

nx[1:] = HRdata.dimensions['x'].size
ny[1:] = HRdata.dimensions['y'].size
nz[1:] = HRdata.dimensions['depth'+grid2plt.lower()].size
nt[1:] = HRdata.dimensions['time_counter'].size

nx[1:] = HRdata.dimensions[HRvardims[-1]].size
ny[1:] = HRdata.dimensions[HRvardims[-2]].size
nt[1:] = HRdata.dimensions[HRvardims[0]].size


# Read corresponding to the Low resolution models
# -----------------------------------------------
idt = 1
idk = 1
# %% Reading LR variables 
if np.size(HRvardims) == 3:
    HRvar = np.mean(HRdata.variables["sossheig"][:idt, :, :], axis=0)
elif np.size(HRvardims) == 4:
    HRvar = np.mean(HRdata.variables["sossheig"][:idt, idk, :, :], axis=0)




ssh_f =  ( HRvar[0:-1,0:-1] + HRvar[1:,0:-1] + HRvar[1:,1:] + HRvar[0:-1,1:] ) * 0.25

f_onU =  ( HR_fft[0:-1,1:] + HR_fft[1:,1:] + HR_fff[1:,0:-1] + HR_fff[1:,1:] ) * 0.25
f_onV =  ( HR_fff[0:-1,1:] + HR_fff[1:,1:] + HR_fft[1:,0:-1] + HR_fft[1:,1:] ) * 0.25

print(np.shape(f_onU))
print(np.shape(f_onV))
print("")

geoU = - ( ssh_f[:,1:] - ssh_f[:, 0:-1] ) / (HRe2u[0:-1,1:-1] * f_onU[:,:-1])
geoV = + ( ssh_f[1:,:] - ssh_f[0:-1,:] ) / (HRe1v[1:-1,0:-1] * f_onV[:-1,:])

geoMod = np.sqrt( (geoU[1:,:] + geoU[0:-1,:])**2 + (geoV[:,1:] + geoV[:, 0:-1])**2 ) 

print(np.shape(HRvar))
print("")
print(np.shape(ssh_f))
print(np.shape(geoU))
print(np.shape(geoV))
print(np.shape(geoMod))





plt.rcParams.update({
	"text.usetex": True,
	"font.family": "serif",
	"font.serif": ["Helvetica"],
	"figure.dpi": 300
})


titles = ["No projection", "With projection"]


# %% Vertical plots

fig, axes = plt.subplots(1, 4)
fig.set_size_inches(4.8, 3)

im2 = axes[0].imshow(np.rot90(HRvar[::-1,:]), cmap='RdBu_r', vmin=-0.5, vmax=0.5)

axes[1].imshow(np.rot90(geoU[::-1,:]), cmap='RdBu_r', vmin=-0.0675, vmax=0.0675)

axes[2].imshow(np.rot90(geoV[::-1,:]), cmap='RdBu_r', vmin=-0.0675, vmax=0.0675)

axes[3].imshow(np.rot90(geoMod[::-1,:]), cmap='RdBu_r', vmin=0.05)


axes_list = axes.flat
for idx, ax in enumerate(axes_list):
    ax.tick_params(which='major', label1On=False)
    ax.tick_params(which='major', axis='both', width=0.25,  direction="in")
    ax.tick_params(which='minor', axis='both', width=0.125, direction='in')

    ax.tick_params(which='major', top=True, bottom=True, left=True, right=True)
    ax.tick_params(which='minor', top=True, bottom=True, left=True, right=True)



    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(0.5)
	
"""
plt.savefig(OUTdir + 'GeostrophicVelocities.png',
            format='png', dpi='figure', bbox_inches='tight')

"""
plt.show()

