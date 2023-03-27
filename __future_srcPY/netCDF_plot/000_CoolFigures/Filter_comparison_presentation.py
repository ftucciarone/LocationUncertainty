#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 11:16:37 2022

@author: ftucciar
"""
# Load modules
# ------------
import numpy as np
import scipy.io.netcdf as nc
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.patches
from matplotlib.colors import LinearSegmentedColormap


cmap = LinearSegmentedColormap.from_list('mycmap',
                                         [  # [204/255, 253/255, 196/255],
                                          [  2/255,  25/255, 186/255],
                                          [255/255, 255/255, 255/255],
                                          [180/255,  46/255,  26/255]])
                                          #  [253/255, 249/255,  83/255]])
plt.rcParams['text.latex.preamble'] = [
    r'\usepackage{amsmath}',
    r'\usepackage{amssymb}']

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Helvetica"],
    "figure.dpi": 300
    })


                           
# Set parammeters
# ---------------

# HR inputs
base_dir = "/home/ftucciar/LocationUncertainty/CoarseGraining/data/"
subs_dir = ['100-102y']  
basefile = 'ocref_r3'
datagrid = ["u", "v", "w"]
suffix = ['_BtBcSpt_36w1_18w2', '_BtBcSpt_72w1_18w2'] 




grid2plt = ''
ext = '.nc'

ingrid = '/home/ftucciar/dataNEMO/Deter/domain_cfg_R3.nc'
grid = nc.netcdf_file(ingrid, 'r')

idx = 0
# Read corresponding to the High resolution models
# ------------------------------------------------
data = Dataset(base_dir + subs_dir[0] + "/" + basefile + suffix[idx] +
               ext, "r", format="NETCDF4")


# Get domain dimension from domain file
# -------------------------------------
nx = grid.dimensions['x']
ny = grid.dimensions['y']
nt = grid.dimensions['time_counter']


u = np.zeros([len(suffix), ny-2, nx-2])
v = np.zeros([len(suffix), ny-2, nx-2])
w = np.zeros([len(suffix), ny-2, nx-2])
vmag = np.zeros([len(suffix), ny-2, nx-2])


magbound = np.zeros(2)
magbound[0] = 0
magbound[1] = 0.2

vbound = np.zeros(2)
vbound[0] = -5e-05
vbound[1] = 5e-05


idtHR = 1

interp = 'none'

for i, idx in enumerate(suffix):
    print(i, idx)
    data = Dataset(base_dir + subs_dir[0] + "/" + basefile + idx + ext,
                   "r", format="NETCDF4")
    u[i, :, :] = data.variables["uband"][idtHR, 0, 1:-1, 1:-1]
    v[i, :, :] = data.variables["vband"][idtHR, 0, 1:-1, 1:-1]
    w[i, :, :] = data.variables["wband"][idtHR, 0, 1:-1, 1:-1]
    vmag[i, :, :] = np.sqrt(u[i, :, :]**2 + v[i, :, :]**2)

    
fig, (ax1, ax2) = plt.subplots(2, 1)
fig.set_size_inches(5, 5)



axes_list = [ax1, ax2] 
for idx, ax in enumerate(axes_list):
    ax.imshow(np.flipud(vmag[idx, :, :]), interpolation=interp,
           cmap=cmap, vmin=magbound[0], vmax=magbound[1])
    
    
axes_list = [ax1, ax2]
for idx, ax in enumerate(axes_list):
    ax.tick_params(which='major', label1On=False)
    ax.tick_params(which='major', axis='both', width=0.125,  direction="in")
    ax.tick_params(which='minor', axis='both', width=0.0675, direction='in')

    ax.tick_params(which='major', top=True, bottom=True, left=True, right=True)
    ax.tick_params(which='minor', top=True, bottom=True, left=True, right=True)

    # Major ticks
    ax.set_xticks(np.arange(10, nx-1, 10))
    ax.set_yticks(np.arange(10, ny-1, 10))

    # Minor ticks
    ax.set_xticks(np.arange(5, nx-1, 5), minor=True)
    ax.set_yticks(np.arange(5, ny-1, 5), minor=True)

    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(0.5)






plt.savefig('filter_comparison_presentation_3.png', format='png', dpi=600, bbox_inches='tight')










