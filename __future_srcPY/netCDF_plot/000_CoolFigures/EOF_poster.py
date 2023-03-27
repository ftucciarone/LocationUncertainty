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
# Set parammeters
# ---------------

# HR inputs
HRbase_dir = '/home/ftucciar/dataNEMO/Deter/R27d/'
HRsubs_dir = ['100-102y']  # ['bplus100', 'debminus100j', 'nobia100']
HRbasefile = 'GYRE_5d_00010101_00021230_grid_U'
HRname = 'vozocrtx'

# LR inputs
CGbase_dir = '/home/ftucciar/LocationUncertainty/CoarseGraining/data/' 
CGsubs_dir = ['100-102y']
CGbasefile = ['ocref_r3']
CGname = 'ur'

MODbase_dir = '/home/ftucciar/dataNEMO/PROMS/R3lu/EXP00/' 
MODsubs_dir = ['']
MODbasefile = ['spmu']
MODname = 'spat_basis_u_001'

NOIbase_dir = '/home/ftucciar/dataNEMO/PROMS/R3lu/EXP00/' 
NOIsubs_dir = ['']
NOIbasefile = ['GYRE_5d_00010101_00301230_grid_U']
NOIname = 'tluvelu'



grid2plt = ''
ext = '.nc'

ingrid = '/home/ftucciar/dataNEMO/Deter/R3d/domain_cfg_R3.nc'
grid = nc.netcdf_file(ingrid, 'r')

# Read corresponding to the High resolution models
# ------------------------------------------------
HRdata = Dataset(HRbase_dir + HRsubs_dir[0] + "/" + HRbasefile +
                 grid2plt + ext, "r", format="NETCDF4")

# Read corresponding to the Low resolution models
# -----------------------------------------------
CGdata = Dataset(CGbase_dir + CGsubs_dir[0] + "/" + CGbasefile[0] +
                 grid2plt + ext, "r", format="NETCDF4")

MODdata = Dataset(MODbase_dir + MODsubs_dir[0] + "/" + MODbasefile[0] +
                 grid2plt + ext, "r", format="NETCDF4")

NOIdata = Dataset(NOIbase_dir + NOIsubs_dir[0] + "/" + NOIbasefile[0] +
                 grid2plt + ext, "r", format="NETCDF4")



# Get domain dimension from domain file
# -------------------------------------
nx = np.zeros(8)
ny = np.zeros(8)
nz = np.zeros(8)
nt = np.zeros(8)

nx[1:] = grid.dimensions['x']
ny[1:] = grid.dimensions['y']
# nz[1:] = grid.dimensions['depth'+grid2plt.lower()]
nt[1:] = grid.dimensions['time_counter']

nx[0] = HRdata.dimensions['x'].size
ny[0] = HRdata.dimensions['y'].size
nz[0] = HRdata.dimensions['depthu'+grid2plt.lower()].size
nt[0] = HRdata.dimensions['time_counter'].size

rx = np.zeros(8)
ry = np.zeros(8)
rx[1:] = 1
ry[1:] = 1
rx[0] = (nx[0] - 2)/(nx[1] - 2)
ry[0] = (ny[0] - 2)/(nx[1] - 2)

idtLR = 0
idtHR = idtLR #- 3*nt[0]



varR27 = np.rot90(HRdata.variables[HRname][idtHR, 0, 1:-1, 1:-1], k=-1)

varCG = np.rot90(CGdata.variables[CGname][idtLR, 0, 1:-1, 1:-1], k=-1)

varMOD = np.rot90(MODdata.variables[MODname][ 0, 1:-1, 1:-1], k=-1)

varNOI = np.rot90(NOIdata.variables[NOIname][idtLR, 0, 1:-1, 1:-1], k=-1) 



# %%

titles = [r"$\mathrm{R27d}\: [\mathrm{m}]$",
          r"$\mathrm{R3d}\: [\mathrm{m}]$",
          r"$\mathrm{R3LU}\: [\mathrm{m}]$"]
plt.rcParams['text.latex.preamble'] = [
    r'\usepackage{amsmath}',
    r'\usepackage{amssymb}']

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Helvetica"],
    "figure.dpi": 300
    })

fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8) = plt.subplots(1, 8)
fig.set_size_inches(6, 3)

vbound = np.zeros(2)
vbound[0] = -0.75
vbound[1] = 0.75

mbound = np.zeros(2)
mbound[0] = -0.05
mbound[1] = 0.05

nbound = np.zeros(2)
nbound[0] = -0.03
nbound[1] = 0.03


im2 = ax1.imshow(np.flipud(varR27), interpolation='spline16',
           cmap=cmap, vmin=vbound[0], vmax=vbound[1])
ax2.imshow(np.flipud(varCG), interpolation='spline16',
           cmap=cmap, vmin=vbound[0], vmax=vbound[1])
ax3.imshow(np.flipud(varMOD), interpolation='spline16',
           cmap=cmap, vmin=mbound[0], vmax=mbound[1], alpha=0.2)
ax4.imshow(np.flipud(varMOD), interpolation='spline16',
           cmap=cmap, vmin=mbound[0], vmax=mbound[1], alpha=0.4)
ax5.imshow(np.flipud(varMOD), interpolation='spline16',
           cmap=cmap, vmin=mbound[0], vmax=mbound[1], alpha=0.6)
ax6.imshow(np.flipud(varMOD), interpolation='spline16',
           cmap=cmap, vmin=mbound[0], vmax=mbound[1], alpha=0.8)
ax7.imshow(np.flipud(varMOD), interpolation='spline16',
           cmap=cmap, vmin=mbound[0], vmax=mbound[1])
ax8.imshow(np.flipud(varNOI), interpolation='spline16',
           cmap=cmap, vmin=nbound[0], vmax=nbound[1])

box1 = ax1.get_position()
box2 = ax2.get_position()
box3 = ax3.get_position()
box4 = ax4.get_position()
box5 = ax5.get_position()
box6 = ax6.get_position()
box7 = ax7.get_position()
box8 = ax8.get_position()

pos1 = [box1.x0 - 0.050, box1.y0,         box1.width, box1.height]
pos5 = ax4.get_position() # get the original position 
pos3 = [pos5.x0 - 0.030, pos5.y0 + 0.06,  pos5.width, pos5.height]
pos4 = [pos5.x0 - 0.015, pos5.y0 + 0.03,  pos5.width, pos5.height]
pos6 = [pos5.x0 + 0.015, pos5.y0 - 0.03,  pos5.width, pos5.height]
pos7 = [pos5.x0 + 0.030, pos5.y0 - 0.06,  pos5.width, pos5.height]
ax1.set_position(pos1)
ax3.set_position(pos3)
ax4.set_position(pos4)
ax5.set_position(pos5)
ax6.set_position(pos6)
ax7.set_position(pos7)
ax8.set_position(box6)

axes_list = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]
for idx, ax in enumerate(axes_list):
    ax.tick_params(which='major', label1On=False)
    ax.tick_params(which='major', axis='both', width=0.125,  direction="in")
    ax.tick_params(which='minor', axis='both', width=0.0675, direction='in')

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

# Upper Descriptions
plt.figtext(box1.x0 - 0.050 + box1.width/2, box1.y0 + box1.height + 0.1, 
            r"High Resolution",
            fontsize=5,
            horizontalalignment ="center", 
            verticalalignment ="center")
plt.figtext(box2.x0 + box2.width/2,  box1.y0 + box1.height + 0.1, 
            r"High-pass filtered",
            fontsize=5,
            horizontalalignment ="center", 
            verticalalignment ="center")
plt.figtext(box2.x0 + box2.width/2,  box1.y0 + box1.height + 0.075, 
            r"velocity",
            fontsize=5,
            horizontalalignment ="center", 
            verticalalignment ="center")
plt.figtext(box4.x0 + box4.width/2, box1.y0 + box1.height + 0.1, 
            r"Proper Orthogonal", 
            fontsize=5,
            horizontalalignment ="center", 
            verticalalignment ="center")
plt.figtext(box4.x0 + box4.width/2, box1.y0 + box1.height + 0.075, 
            r"Decomposition", 
            fontsize=5,
            horizontalalignment ="center", 
            verticalalignment ="center")
plt.figtext(box6.x0 + box6.width/2, box1.y0 + box1.height + 0.1, 
            r"Noise Forcing",
            fontsize=5,
            horizontalalignment ="center", 
            verticalalignment ="center")

# Lower descriptions
plt.figtext(box1.x0 - 0.050 + box1.width/2, box1.y0-0.1, 
            r"$\mathbf{u}_{_{\mathrm{HR}}}$",
            fontsize=5,
            horizontalalignment ="center", 
            verticalalignment ="center")
plt.figtext(box2.x0 + box2.width/2, box1.y0-0.1, 
            r"$\mathbf{u}_{_{\mathrm{LR}}} = \left[\left( 1 - \mathcal{G}\right)\mathbf{u}_{_{\mathrm{HR}}}\right]^{\downarrow}_{_{\mathrm{LR}}}$",
            fontsize=5,
            horizontalalignment ="center", 
            verticalalignment ="center")
plt.figtext(box4.x0 + box4.width/2, box1.y0-0.1, 
            r"$\left\lbrace \boldsymbol{\Phi}_{k}(\boldsymbol{x}),\: k\in \left[1, N\right]\right\rbrace$", 
            fontsize=5,
            horizontalalignment ="center", 
            verticalalignment ="center")
plt.figtext(box6.x0 + box6.width/2, box1.y0-0.1, 
            r"$\boldsymbol{\dot{\eta}}(\boldsymbol{x}) = \displaystyle\sum_{\tiny k=1}^{\tiny N} \sqrt{\lambda_k} \boldsymbol{\Phi}_k(\boldsymbol{x}) \xi_k$",
            fontsize=5,
            horizontalalignment ="center", 
            verticalalignment ="center")

# Create the arrow
# 1. Get transformation operators for axis and figure
ax0tr = ax1.transData # Axis 0 -> Display
ax1tr = ax2.transData # Axis 1 -> Display
figtr = fig.transFigure.inverted() # Display -> Figure
# 2. Transform arrow start point from axis 0 to figure coordinates
ptB = figtr.transform(ax0tr.transform((540., 405.)))
# 3. Transform arrow end point from axis 1 to figure coordinates
ptE = figtr.transform(ax1tr.transform((0., 45.)))
# 4. Create the patch
arrow = matplotlib.patches.FancyArrowPatch(
    ptB, ptE, transform=fig.transFigure,  # Place arrow in figure coord system
    fc = "w", connectionstyle="arc3,rad=0.0", arrowstyle='simple', alpha = 0.5,
    mutation_scale = 20.
)
# 5. Add patch to list of objects to draw onto the figure
fig.patches.append(arrow)


# Create the arrow
# 1. Get transformation operators for axis and figure
ax0tr = ax2.transData # Axis 0 -> Display
ax1tr = ax5.transData # Axis 1 -> Display
figtr = fig.transFigure.inverted() # Display -> Figure
# 2. Transform arrow start point from axis 0 to figure coordinates
ptB = figtr.transform(ax0tr.transform((60., 45.)))
# 3. Transform arrow end point from axis 1 to figure coordinates
ptE = figtr.transform(ax1tr.transform((-22., 45.)))
# 4. Create the patch
arrow = matplotlib.patches.FancyArrowPatch(
    ptB, ptE, transform=fig.transFigure,  # Place arrow in figure coord system
    fc = "w", connectionstyle="arc3,rad=0.0", arrowstyle='simple', alpha = 0.5,
    mutation_scale = 20.
)
# 5. Add patch to list of objects to draw onto the figure
fig.patches.append(arrow)


# Create the arrow
# 1. Get transformation operators for axis and figure
ax0tr = ax5.transData # Axis 0 -> Display
ax1tr = ax8.transData # Axis 1 -> Display
figtr = fig.transFigure.inverted() # Display -> Figure
# 2. Transform arrow start point from axis 0 to figure coordinates
ptB = figtr.transform(ax0tr.transform((82., 45.)))
# 3. Transform arrow end point from axis 1 to figure coordinates
ptE = figtr.transform(ax1tr.transform((0., 45.)))
# 4. Create the patch
arrow = matplotlib.patches.FancyArrowPatch(
    ptB, ptE, transform=fig.transFigure,  # Place arrow in figure coord system
    fc = "w", connectionstyle="arc3,rad=0.0", arrowstyle='simple', alpha = 0.5,
    mutation_scale = 20.
)
# 5. Add patch to list of objects to draw onto the figure
fig.patches.append(arrow)


fig.patch.set_alpha(1.0)
#fig.text(1.75, 0.5,
#         'Preliminary version', transform=ax1.transAxes,
#         fontsize=30, color='gray', alpha=0.75,
#         ha='center', va='center', rotation='20')


plt.savefig('poster.png', format='png', dpi=600, bbox_inches='tight')

