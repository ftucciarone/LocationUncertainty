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
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.transforms import Affine2D

cmap = LinearSegmentedColormap.from_list('mycmap',
                                         [  # [204/255, 253/255, 196/255],
                                          [  2/255,  25/255, 186/255],
                                          [255/255, 255/255, 255/255],
                                          [180/255,  46/255,  26/255]])
                                          #  [253/255, 249/255,  83/255]])
cold = (  2/255,  25/255, 186/255)
warm = (180/255,  46/255,  26/255)
# Set parammeters
# ---------------

# HR inputs
HRbase_dir = '/Users/Francesco/dataNEMO/'
HRsubs_dir = ['R27']  # ['bplus100', 'debminus100j', 'nobia100']
HRbasefile = 'GYRE_5d_00010101_00021230_grid_'

grid2plt = 'T'
ext = '.nc'

ingrid = '/Users/Francesco/dataNEMO/R27/domain_cfg_R27.nc'
grid = nc.netcdf_file(ingrid, 'r')

# Read corresponding to the High resolution models
# ------------------------------------------------
HRdataT = Dataset(HRbase_dir + HRsubs_dir[0] + "/" + HRbasefile +
                 'T' + ext, "r", format="NETCDF4")

HRgrid = Dataset(ingrid, "r", format="NETCDF4")

# Get domain dimension from domain file
# -------------------------------------
nx = HRdataT.dimensions['x'].size
ny = HRdataT.dimensions['y'].size
nz = HRdataT.dimensions['depth'+grid2plt.lower()].size
nt = HRdataT.dimensions['time_counter'].size

# Get geographical informations
# -----------------------------
lat = HRgrid.variables['nav_lat'][1:-1, 1:-1]
lon = HRgrid.variables['nav_lon'][1:-1, 1:-1]

maxLat = lat.max()
maxLon = lon.max()
minLat = lat.min()
minLon = lon.min()


idtLR = 376
idtHR = idtLR - 3*nt

varnameR27 = 'sossheig'
varR27 = np.rot90(HRdataT.variables[varnameR27][idtHR, 1:-1, 1:-1], k=-1)








Lat = np.arange(minLat, maxLat, 0.1)



# -----------------------------------------------------------------------------
#     Wind Stress Anlitical Definition (With plot bounds)
# -----------------------------------------------------------------------------
maxtau = 0.105
mintau = maxtau - 0.015
windmax = maxtau * np.sin( np.pi * (Lat - 15.) / (29.-15.) )
windmin = mintau * np.sin( np.pi * (Lat - 15.) / (29.-15.) )
maxW = maxtau + 0.05 
minW = - maxW
windtick = [-maxtau, 0, maxtau]
windticklabels = [-0.1, 0, 0.1]
windUnits =r"$\tau [ \mathrm{N}/\mathrm{m} ]$"

# -----------------------------------------------------------------------------
#     Freshwater flux Anlitical Definition (With plot bounds)
# -----------------------------------------------------------------------------
maxemp = 0.105 
minemp = maxtau - 0.015
empmax = -maxtau * np.sin( np.pi * (Lat - 15.) / (29.-15.) )
empmin = -mintau * np.sin( np.pi * (Lat - 15.) / (29.-15.) )
maxE = maxemp + 0.05 
minE = - maxE
emptick = [-maxtau, 0, maxtau]
empticklabels = [-0.1, 0, 0.1]
empUnits =r"$F [ \mathrm{m}/\mathrm{y} ]$" 


# -----------------------------------------------------------------------------
#     Apparent Air-Temperature Anlitical Definition (With plot bounds)
# -----------------------------------------------------------------------------
maxT = 28.3 
minT = 0.0
summer = 1
winter = 0
Tmpmax = maxT * (1 + 1/50 * summer) * np.cos( np.pi * (Lat - 5.) / (53 * (1 + 11/53 * summer) * 2) )
Tmpmin = maxT * (1 + 1/50 * winter) * np.cos( np.pi * (Lat - 5.) / (53 * (1 + 11/53 * winter) * 2) )
maxT = 32 
minT = 2
Tmptick = [5, 30]
Tmpticklabels = [5, 30]
tempUnits = r"$T^{\star} [ ^{\circ}\mathrm{C} ]$" 


# -----------------------------------------------------------------------------
#     Radiation Forcing Anlitical Definition (With plot bounds)
# -----------------------------------------------------------------------------
maxRad = 230
minRad = 0.0
summer = 1
winter = 0
Radmax = maxRad * np.cos( np.pi * (Lat - 23.5 * summer) / (0.9 - 180.) )
Radmin = maxRad * np.cos( np.pi * (Lat - 23.5 * winter) / (0.9 - 180.) )
maxR = 235
minR = 135
Radtick = [140, 230]
Radticklabels = [140, 230]
radUnits = r"$Q [ \mathrm{W}/\mathrm{m}^{2} ]$"





























# Add 2D affine transformation
t1 = Affine2D().rotate_deg(45)


start = 0
step = 60 
stopX = 540
stopY = 810



# axes_list = [ax1, ax2]
# for idx, ax in enumerate(axes_list):
#     # ax.tick_params(which='major', label1On=False)
#     # ax.tick_params(which='major', axis='both', width=0.125,  direction="in")
#     # ax.tick_params(which='minor', axis='both', width=0.0675, direction='in')

#     # ax.tick_params(which='major', top=True, bottom=True, left=True, right=True)
#     # ax.tick_params(which='minor', top=True, bottom=True, left=True, right=True)

#     for axis in ['top', 'bottom', 'left', 'right']:
#         ax.spines[axis].set_linewidth(0.5)
        


# %

minX = -575
maxX = 389
minY = 955
maxY = -5

# Construct geographical Ticks
spanLonTicks = np.linspace(minX, maxX, 5)
LonTicksLabels = np.linspace(np.abs(minLon), np.abs(maxLon), 5).astype(int)
LonTicksLabels = list(map(str, LonTicksLabels))
LonTicksLabels = [tick_label + "W" for tick_label in LonTicksLabels]

spanLatTicks = np.linspace(minY, maxY, 5)
LatTicksLabels = np.linspace(minLat, maxLat, 5).astype(int)
LatTicksLabels = list(map(str, LatTicksLabels))
LatTicksLabels = [tick_label + "N" for tick_label in LatTicksLabels]




plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Helvetica"],
    "figure.figsize" : (12/2.54, 8/2.54),
    "figure.dpi": 300,
    "axes.labelsize": 6,
    "axes.titlesize": 6,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6
    })


fig=plt.figure(1)

ax1 = plt.subplot2grid((4, 7), (0, 0), rowspan=3, colspan=3)
ax1.imshow(np.flipud(varR27),
            interpolation='spline16', cmap=cmap,  transform=t1 + ax1.transData)
ax1.text(-325, 275, '3180 km', rotation=45, 
        horizontalalignment='center',
        verticalalignment='center', fontsize=9)
ax1.text(225, 175, '2120 km' , rotation=-45, 
        horizontalalignment='center',
        verticalalignment='center', fontsize=9)
minX = -575
maxX = 389
minY = 955
maxY = -5

ax1.set_xlim([minX,maxX])
ax1.set_ylim([minY,maxY])
ax1.tick_params(top=True, bottom=False, 
                left=True, right=False,
                labelbottom=False,  labeltop=True,
                labelleft=True, labelright=False)
ax1.set_xticks(spanLonTicks)
ax1.set_yticks(spanLatTicks)
ax1.set_xticklabels(LonTicksLabels, rotation = -90, ha='right', va='bottom')
ax1.set_yticklabels(LatTicksLabels, rotation =   0, ha='right', va='bottom')

box1 = ax1.get_position()


ax2 = plt.subplot2grid((4, 7), (0, 3), rowspan=3, colspan=1)
ax2.plot(windmax, Lat, color=cold, linewidth=0.75)
ax2.plot(windmin, Lat, color=warm, linewidth=0.75)
ax2.set_xlim([minW,maxW])
ax2.set_ylim([minY,maxY])
ax2.set_yticks(spanLatTicks)
ax2.set_ylim(minLat, maxLat)
ax2.tick_params(top=False, bottom=True, 
                left=True, right=False,
                labelbottom=False,  labeltop=True,
                labelleft=False, labelright=False)
ax2.set_xticks(windtick)
ax2.set_xticklabels(windticklabels, rotation = -90, ha='center', va='bottom')
box2 = ax2.get_position()
ax2.set_position([box2.x0, box1.y0, box2.width, box1.height])
ax2.set_xlabel(windUnits)
ax2.xaxis.grid(True, linewidth=0.5)


ax3 = plt.subplot2grid((4, 7), (0, 4), rowspan=3, colspan=1)
ax3.plot(empmin, Lat, color=cold, linewidth=0.75)
ax3.plot(empmax, Lat, color=warm, linewidth=0.75)
ax3.set_xlim([minE,maxE])
ax3.set_ylim([minY,maxY])
ax3.set_yticks(spanLatTicks)
ax3.set_ylim(minLat, maxLat)
ax3.tick_params(top=False, bottom=True, 
                left=True, right=False,
                labelbottom=False,  labeltop=True,
                labelleft=False, labelright=False)
ax3.set_xticks(emptick)
ax3.set_xticklabels(empticklabels, rotation = -90, ha='center', va='bottom')
box3 = ax3.get_position()
ax3.set_position([box3.x0, box1.y0, box3.width, box1.height])
ax3.set_xlabel(empUnits)
ax3.xaxis.grid(True, linewidth=0.5)

ax4 = plt.subplot2grid((4, 7), (0, 5), rowspan=3, colspan=1)
ax4.plot(Tmpmin, Lat, color=cold, linewidth=0.75)
ax4.plot(Tmpmax, Lat, color=warm, linewidth=0.75)
ax4.set_xlim([minT,maxT])
ax4.set_ylim([minY,maxY])
ax4.set_yticks(spanLatTicks)
ax4.set_ylim(minLat, maxLat)
ax4.tick_params(top=False, bottom=True, 
                left=True, right=False,
                labelbottom=False,  labeltop=True,
                labelleft=False, labelright=False)
ax4.set_xticks(Tmptick)
ax4.set_xticklabels(Tmpticklabels, rotation = -90, ha='center', va='bottom')
box4 = ax4.get_position()
ax4.set_position([box4.x0, box1.y0, box4.width, box1.height])
ax4.set_xlabel(tempUnits)
ax4.xaxis.grid(True, linewidth=0.5)

ax5 = plt.subplot2grid((4, 7), (0, 6), rowspan=3, colspan=1)
ax5.plot(Radmin, Lat, color=cold, linewidth=0.75)
ax5.plot(Radmax, Lat, color=warm, linewidth=0.75)
ax5.set_xlim([minR,maxR])
ax5.set_ylim([minY,maxY])
ax5.set_yticks(spanLatTicks)
ax5.set_ylim(minLat, maxLat)
ax5.tick_params(top=False, bottom=True, 
                left=True, right=False,
                labelbottom=False,  labeltop=True,
                labelleft=False, labelright=False)
ax5.set_xticks(Radtick)
ax5.set_xticklabels(Radticklabels, rotation = -90, ha='center', va='bottom')
box5 = ax5.get_position()
ax5.set_position([box5.x0, box1.y0, box5.width, box1.height])
ax5.set_xlabel(radUnits)
ax5.xaxis.grid(True, linewidth=0.5)

ax_list = fig.axes
for ax in ax_list:
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(0.5)
    ax.tick_params(axis="y",direction="in", width=0.5, length=2)
    ax.tick_params(axis="x",direction="in", width=0.5, length=2)


plt.savefig("Experimental_setup.pdf", dpi=fig.dpi, bbox_inches='tight')
# plt.savefig("Comparison_RMS.eps", dpi=fig.dpi, bbox_inches='tight')
plt.show()




















