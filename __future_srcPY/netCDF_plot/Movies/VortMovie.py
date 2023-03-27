#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 15:00:04 2022

@author: ftucciar
"""



import numpy as np
import scipy.io.netcdf as nc
import seaborn as sns
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import LinearSegmentedColormap
# This is needed to display graphics calculated outside of jupyter notebook
from IPython.display import HTML, display
cmap = LinearSegmentedColormap.from_list('mycmap',
                                         [
                                          [2/255,  25/255, 186/255],
                                          [255/255, 255/255, 255/255],
                                          [180/255,  46/255,  26/255]])



# Set parameters
# ---------------
# Output directory
OUTdir = "/home/ftucciar/LU_postprocessing/netCDF_plot/outputs/"

# LR inputs
LRbase_dir = '/home/ftucciar/dataNEMO/'
LRsubs_dir = ['Deter/R3d', 
              'final/36w1_144w2']
LRprnt_lab = [r"$u_{\mathrm{R27}}$",
              r"$u_{\mathrm{R3}}$",
              r"$( \mathcal{F}_{36} - \mathcal{F}_{144} )u$"]
LRbasefile = ['GYRE_5d_00010101_00301230_grid_',
              'GYRE_5d_00010101_00101230_grid_']

HRbase_dir = '/home/ftucciar/dataNEMO/Deter/R27d/'
HRsubs_dir = ['100-102y','102-104y','104-106y','106-108y','108-110y']

HRbasefile = ['curl_']

grid2plt = 'T'
ext = '.nc'

LRgrid = '/home/ftucciar/dataNEMO/Deter/domain_cfg_R3.nc'
LRgrid = Dataset(LRgrid, "r", format="NETCDF4")
HRgrid = '/home/ftucciar/dataNEMO/Deter/domain_cfg_R27.nc'
HRgrid = Dataset(HRgrid, "r", format="NETCDF4")

HRvarname = 'socurlt'
LRvarname = 'vorticity'

# Read corresponding to the Low resolution models
# -----------------------------------------------
HRdata = Dataset(HRbase_dir + HRsubs_dir[0] + "/" + HRbasefile[0] +
                 grid2plt + ext, "r", format="NETCDF4")
DTdata = Dataset(LRbase_dir + LRsubs_dir[0] + "/" + LRbasefile[0] +
                 grid2plt + ext, "r", format="NETCDF4")
LU1data = Dataset(LRbase_dir + LRsubs_dir[1] + "/" + LRbasefile[1] +
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
# W-points (for integration) --------------------------------------------
LRe1t = LRgrid.variables['e1t'][0, :, :].data
LRe2t = LRgrid.variables['e2t'][0, :, :].data
LRe3t = LRgrid.variables['e3t_0'][0, :, :].data
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

HR_fft = HRgrid.variables['ff_t'][0, :, :].data
# Free memory
del LRe1u, LRe2u, LRe1v, LRe2v, LRe1t, LRe2t

LRvardims = DTdata.variables[LRvarname].dimensions

# %% Reading LR temperature
if np.size(LRvardims) == 3:
    DTvar = DTdata.variables['votemper'][:, :, :]
    LU1var = np.mean(LU1data.variables['vorticity'][:, :, :], axis=0)
elif np.size(LRvardims) == 4:
    DTvar = np.divide(DTdata.variables['vorticity'][:, 0, :, :],
                      LR_fft)
    LU1var = np.divide(LU1data.variables['vorticity'][:, 0, :, :],
                      LR_fft)

HRvar = np.empty([0,542,812])
for s in HRsubs_dir:
    print(s)
    fin = Dataset(HRbase_dir + s + "/" + HRbasefile[0] + 
                  grid2plt + ext, "r", format="NETCDF4")
    HRvar = np.append(HRvar,np.divide(HRdata.variables[HRvarname][:, 0, :, :],
                      HR_fft),axis=0)
    fin.close()



# %%
idz = 1

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Helvetica"],
    "figure.dpi": 300
    })

vbound = np.zeros(2)
vbound[0] = -0.2
vbound[1] =  0.2

# Get a handle on the figure and the axes
fig, (ax0, ax1, ax2) = plt.subplots(1,3, figsize=(12,6))





# Plot the initial frame. 
ax0.imshow(np.flipud(np.rot90(HRvar[0, :, :], k=-1)),
           interpolation='none', cmap=cmap, vmin=vbound[0], vmax=vbound[1])
ax1.imshow(np.flipud(np.rot90(DTvar[0, :, :], k=-1)),
           interpolation='none', cmap=cmap, vmin=vbound[0], vmax=vbound[1])
ax2.imshow(np.flipud(np.rot90(LU1var[0, :, :], k=-1)),
           interpolation='none', cmap=cmap, vmin=vbound[0], vmax=vbound[1])

ax0.tick_params(which='major', label1On=False)
ax0.tick_params(which='major', axis='both', width=0.125,  direction="in")
ax0.tick_params(which='minor', axis='both', width=0.0675, direction='in')

ax0.tick_params(which='major', top=True, bottom=True, left=True, right=True)
ax0.tick_params(which='minor', top=True, bottom=True, left=True, right=True)


for axis in ['top', 'bottom', 'left', 'right']:
    ax0.spines[axis].set_linewidth(0.5)
    
ax0.set_title(LRprnt_lab[0], fontsize=20)

ax1.tick_params(which='major', label1On=False)
ax1.tick_params(which='major', axis='both', width=0.125,  direction="in")
ax1.tick_params(which='minor', axis='both', width=0.0675, direction='in')

ax1.tick_params(which='major', top=True, bottom=True, left=True, right=True)
ax1.tick_params(which='minor', top=True, bottom=True, left=True, right=True)


for axis in ['top', 'bottom', 'left', 'right']:
    ax1.spines[axis].set_linewidth(0.5)
    
ax1.set_title(LRprnt_lab[1], fontsize=20)

ax2.tick_params(which='major', label1On=False)
ax2.tick_params(which='major', axis='both', width=0.125,  direction="in")
ax2.tick_params(which='minor', axis='both', width=0.0675, direction='in')

ax2.tick_params(which='major', top=True, bottom=True, left=True, right=True)
ax2.tick_params(which='minor', top=True, bottom=True, left=True, right=True)


for axis in ['top', 'bottom', 'left', 'right']:
    ax2.spines[axis].set_linewidth(0.5)
    
ax2.set_title(LRprnt_lab[2], fontsize=20)
        

# Next we need to create a function that updates the values for the colormesh, as well as the title.
def animate(frame):
    print("frame " + str(frame))
    ax0.imshow(np.flipud(np.rot90(HRvar[frame, :, :], k=-1)),
           interpolation='none', cmap=cmap, vmin=vbound[0], vmax=vbound[1])
    ax1.imshow(np.flipud(np.rot90(DTvar[frame, :, :], k=-1)),
               interpolation='none', cmap=cmap, vmin=vbound[0], vmax=vbound[1])
    ax2.imshow(np.flipud(np.rot90(LU1var[frame, :, :], k=-1)),
               interpolation='none', cmap=cmap, vmin=vbound[0], vmax=vbound[1])


# Finally, we use the animation module to create the animation.
ani = FuncAnimation(
    fig,             # figure
    animate,         # name of the function above
    frames=720,       # Could also be iterable or list
    interval=100     # ms between frames
)

ani.save('Python_Animation_04.mp4')
display(HTML("<video controls><source src='Python_Animation_04.mp4' type='video/mp4'></video>"))





