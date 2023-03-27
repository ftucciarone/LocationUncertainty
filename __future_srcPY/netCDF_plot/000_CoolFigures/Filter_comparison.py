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
suffix = ['_BtBcSpt_0w1_36w2', '_BtBcSpt_36w1_72w2', '_BtBcSpt_36w1_144w2'] 




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


u = np.zeros([len(suffix), nx-2, ny-2])
v = np.zeros([len(suffix), nx-2, ny-2])
w = np.zeros([len(suffix), nx-2, ny-2])
vmag = np.zeros([len(suffix), nx-2, ny-2])

phiu = np.zeros([len(suffix), nx-2, ny-2])
phiv = np.zeros([len(suffix), nx-2, ny-2])
phiw = np.zeros([len(suffix), nx-2, ny-2])
phimag = np.zeros([len(suffix), nx-2, ny-2])

magbound = np.zeros(2)
magbound[0] = 0
magbound[1] = 0.2

vbound = np.zeros(2)
vbound[0] = -5e-05
vbound[1] = 5e-05

phimagbound = np.zeros(2)
phimagbound[0] = 0
phimagbound[1] = 0.1

phivbound = np.zeros(2)
phivbound[0] = -1e-05
phivbound[1] = 1e-05
idtHR = 1

interp = 'none'

for i, idx in enumerate(suffix):
    print(i, idx)
    data = Dataset(base_dir + subs_dir[0] + "/" + basefile + idx + ext,
                   "r", format="NETCDF4")
    u[i, :, :] = np.rot90(data.variables["uband"][idtHR, 0, 1:-1, 1:-1], k=-1)
    v[i, :, :] = np.rot90(data.variables["vband"][idtHR, 0, 1:-1, 1:-1], k=-1)
    w[i, :, :] = np.rot90(data.variables["wband"][idtHR, 0, 1:-1, 1:-1], k=-1)
    vmag[i, :, :] = np.sqrt(u[i, :, :]**2 + v[i, :, :]**2)

    phidata = Dataset(base_dir + "3Dpod" + "/" + "spmu" + idx +
                   "_weighted" + ext, "r", format="NETCDF4")
    phiu[i, :, :] = np.rot90(phidata.variables["spat_basis_u_001"][0, 1:-1, 1:-1], k=-1)
    phidata = Dataset(base_dir + "3Dpod" + "/" + "spmv" + idx +
                   "_weighted" + ext, "r", format="NETCDF4")
    phiv[i, :, :] = np.rot90(phidata.variables["spat_basis_v_001"][0, 1:-1, 1:-1], k=-1)
    phidata = Dataset(base_dir + "3Dpod" + "/" + "spmw" + idx +
                   "_weighted" + ext, "r", format="NETCDF4")
    phiw[i, :, :] = np.rot90(phidata.variables["spat_basis_w_001"][0, 1:-1, 1:-1], k=-1)
    phimag[i, :, :] = np.sqrt(phiu[i, :, :]**2 + phiv[i, :, :]**2)
    
fig, ([ax1, ax2, ax7, ax8], [ax3, ax4, ax9, ax10], [ax5, ax6, ax11, ax12]) = plt.subplots(3, 4)
fig.set_size_inches(5, 5)



axes_list = [ax1, ax3, ax5] 
for idx, ax in enumerate(axes_list):
    ax.imshow(np.flipud(vmag[idx, :, :]), interpolation=interp,
           cmap=cmap, vmin=magbound[0], vmax=magbound[1])
    
axes_list = [ax2, ax4, ax6]
for idx, ax in enumerate(axes_list):
    ax.imshow(np.flipud(w[idx, :, :]), interpolation=interp,
           cmap=cmap, vmin=vbound[0], vmax=vbound[1])   

axes_list = [ax7, ax9, ax11] 
for idx, ax in enumerate(axes_list):
    ax.imshow(np.flipud(phimag[idx, :, :]), interpolation=interp,
           cmap=cmap, vmin=phimagbound[0], vmax=phimagbound[1])
    
axes_list = [ax8, ax10, ax12]
for idx, ax in enumerate(axes_list):
    ax.imshow(np.flipud(phiw[idx, :, :]), interpolation=interp,
           cmap=cmap, vmin=phivbound[0], vmax=phivbound[1])   
    
axes_list = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12]
for idx, ax in enumerate(axes_list):
    ax.tick_params(which='major', label1On=False)
    ax.tick_params(which='major', axis='both', width=0.125,  direction="in")
    ax.tick_params(which='minor', axis='both', width=0.0675, direction='in')

    ax.tick_params(which='major', top=True, bottom=True, left=True, right=True)
    ax.tick_params(which='minor', top=True, bottom=True, left=True, right=True)

    # Major ticks
    ax.set_xticks(np.arange(10, ny-1, 10))
    ax.set_yticks(np.arange(10, nx-1, 10))

    # Minor ticks
    ax.set_xticks(np.arange(5, ny-1, 5), minor=True)
    ax.set_yticks(np.arange(5, nx-1, 5), minor=True)

    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(0.5)


box1 = ax1.get_position()
box2 = ax2.get_position()
box3 = ax3.get_position()
box4 = ax4.get_position()
box5 = ax5.get_position()
box6 = ax6.get_position()
box7 = ax7.get_position()
box8 = ax8.get_position()

# Upper Descriptions
plt.figtext(box1.x0 + box1.width/2, box1.y0 + box1.height + 0.02, 
            r"$\sqrt{u^{2} + v^{2}}$",
            fontsize=5,
            horizontalalignment ="center", 
            verticalalignment ="center")
plt.figtext(box2.x0 + box2.width/2, box2.y0 + box2.height + 0.02, 
            r"$w$",
            fontsize=5,
            horizontalalignment ="center", 
            verticalalignment ="center")
plt.figtext(box7.x0 + box7.width/2, box7.y0 + box7.height + 0.02, 
            r"$\sqrt{\Phi^{2}_{x} + \Phi^{2}_{y}}$",
            fontsize=5,
            horizontalalignment ="center", 
            verticalalignment ="center")
plt.figtext(box8.x0 + box8.width/2, box8.y0 + box8.height + 0.02, 
            r"$\Phi_{z}$",
            fontsize=5,
            horizontalalignment ="center", 
            verticalalignment ="center")
plt.figtext(box1.x0 - 0.05, box1.y0 + box1.height/2, 
            r"$\left( I - \mathcal{F}_{36} \right) u^{\downarrow}$",
            fontsize=5,
            horizontalalignment ="center", 
            verticalalignment ="center",
            rotation=90)

plt.figtext(box3.x0 - 0.05, box3.y0 + box3.height/2, 
            r"$\left( \mathcal{F}_{36} - \mathcal{F}_{72} \right) u^{\downarrow}$",
            fontsize=5,
            horizontalalignment ="center", 
            verticalalignment ="center",
            rotation=90)

plt.figtext(box5.x0 - 0.05, box5.y0 + box5.height/2, 
            r"$\left( \mathcal{F}_{36} - \mathcal{F}_{144} \right) u^{\downarrow}$",
            fontsize=5,
            horizontalalignment ="center", 
            verticalalignment ="center",
            rotation=90)




plt.savefig('filter_comparison.png', format='png', dpi=600, bbox_inches='tight')










