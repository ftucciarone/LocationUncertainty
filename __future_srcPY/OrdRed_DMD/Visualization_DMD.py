#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 15:14:29 2022

@author: ftucciar
"""
# Load modules
# ------------
import os
import sys
import time
import pathlib
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

plt.ion()

def curl(u_in, v_in):
    dx = 35000
    dy = 35000
    curl_out = ( v_in[:-1,1:]- v_in[:-1,:-1])/dx - ( u_in[1:,:-1] - u_in[:-1,:-1] )/dy 
    return curl_out


def smooth(array_in, wd, mu, sigma):
    
    array_out = np.pad(array_in, ((wd,wd),(wd,wd)), constant_values=0)
    x, y = np.meshgrid(np.linspace(-1, 1, 2*wd + 1),
                       np.linspace(-1, 1, 2*wd + 1))
    dst = np.sqrt(x**2 + y**2)
    normal = 1. / (2.0 * np.pi * sigma**2)
    gauss = np.exp(-((dst-mu)**2 / (2.0 * sigma**2 ))) * normal
    
    array_out = gaussian_filter(array_out, sigma=sigma)
    return array_out[wd:-wd, wd:-wd]



vorticity = True
#vorticity = False

# Set directories
# ---------------
# home_dir: depending on the system
# work_dir: folder containing the LU procedures
# data_dir: containing the i/o data
#  dom_dir: containing the domain files from NEMO
#  out_dir: output parent directory
home_dir = str(pathlib.Path.home())
work_dir = home_dir + "/__future_LocationUncertainty"
data_dir = work_dir + "/data/DMD_Output"
dom_dir = work_dir + "/data"
out_dir = work_dir + "/data/DMD_Output"
print(" isdir=", os.path.isdir(home_dir), "home_dir = ", home_dir)
print(" isdir=", os.path.isdir(work_dir), "work_dir = ", work_dir)
print(" isdir=", os.path.isdir(data_dir), "data_dir = ", data_dir)
print(" isdir=", os.path.isdir(dom_dir), " dom_dir = ", dom_dir)
print(" isdir=", os.path.isdir(out_dir), " out_dir = ", out_dir)
print("")

if len(sys.argv[1:])<4:
        print("                 ERROR: not enough arguments to run the script")
        print("Correct syntax is:",sys.argv[0], "-prexix PR -weight XXw1_XXw2")
        sys.exit(0)
# Set name of the input and output files
# --------------------------------------
# prefix: choice of the operation on the coarse-grained fluctuation (ke for rescaled)
prefix = sys.argv[2] #"keLR" #"keLR"   # Argument of first input
suffix = sys.argv[4] #"36w1_144w2"     # Argument of second input
infile = data_dir + "/" + prefix + "_ocludat_r9_" + suffix + ".nc"
ingrid = dom_dir + "/domain_cfg_R3.nc"
print("isfile=", os.path.isfile(ingrid), "          Opening domain file:", ingrid)
print("isfile=", os.path.isfile(infile), "Opening coarse-grained fields:", infile)
print("")

idk = 0
# %% LR inputs
fin = Dataset(infile, "r", format="NETCDF4")

omega_r = fin.variables['omega_r'][::2].data
omega_c = fin.variables['omega_c'][::2].data

# Velocity time average
umean = fin.variables['uco'][idk, :, :].data
vmean = fin.variables['vco'][idk, :, :].data
avg_vort = curl(umean, vmean)
avg_vort[:, -1] = 1 # Just to remind yourself what color is positive 
# wmean = fin.variables['wco'][:, :, :].data
#avg_vort[:, :] = 0 
# Correlated part
umod_rp_cor = fin.variables['umode_real_c'][::2, idk, :, :].data
vmod_rp_cor = fin.variables['vmode_real_c'][::2, idk, :, :].data
# wmod_rp_cor = fin.variables['wmode_real_c'][::2, :, :, :].data
umod_ip_cor = fin.variables['umode_imag_c'][::2, idk, :, :].data
vmod_ip_cor = fin.variables['vmode_imag_c'][::2, idk, :, :].data
# wmod_ip_cor = fin.variables['wmode_imag_c'][::2, :, :, :].data

# Random part
umod_rp_rnd = fin.variables['umode_real_r'][::2, idk, :, :].data
vmod_rp_rnd = fin.variables['vmode_real_r'][::2, idk, :, :].data
# wmod_rp_rnd = fin.variables['wmode_real_r'][::2, :, :, :].data
umod_ip_rnd = fin.variables['umode_imag_r'][::2, idk, :, :].data
vmod_ip_rnd = fin.variables['vmode_imag_r'][::2, idk, :, :].data
# wmod_ip_rnd = fin.variables['wmode_imag_r'][::2, :, :, :].data

# %% Generate the noise
start = 86400*360*10
stop = 86400*360*30
step = 86400*5


# avg_vort = smooth(avg_vort, 1, 0, 1)

k = 1
if vorticity: 
    fig, axes = plt.subplots(1,4)
    fig.set_size_inches(9.6, 6)
    maxval = 0.000005
else:
    fig, axes = plt.subplots(3, 4)
    fig.set_size_inches(15.2, 9)
    maxval = 0.5*np.max(np.sqrt(umean**2 + vmean**2))

axes_list = axes.flat

for t in range(start,stop,step):
    for i in range(np.shape(umod_ip_cor)[0]):
        u_cor = 2*np.cos(omega_c[i]*t) * umod_rp_cor[i,:,:] - \
                2*np.sin(omega_c[i]*t) * umod_ip_cor[i,:,:]
        v_cor = 2*np.cos(omega_c[i]*t) * vmod_rp_cor[i,:,:] - \
                2*np.sin(omega_c[i]*t) * vmod_ip_cor[i,:,:]
    for i in range(np.shape(umod_ip_rnd)[0]):
        u_rnd = 2*np.cos(omega_r[i]*t) * umod_rp_rnd[i,:,:] - \
                2*np.sin(omega_r[i]*t) * umod_ip_rnd[i,:,:]
        v_rnd = 2*np.cos(omega_r[i]*t) * vmod_rp_rnd[i,:,:] - \
                2*np.sin(omega_r[i]*t) * vmod_ip_rnd[i,:,:]
    u = umean + u_cor + u_rnd
    v = vmean + v_cor + v_rnd

    if ( k % 1 == 0):
        k = 0
        if vorticity:
            cor_vort = curl(u_cor, v_cor)
            rnd_vort = curl(u_rnd, v_rnd)
            vort =  cor_vort + rnd_vort
           # vort =  avg_vort + cor_vort + rnd_vort
            axes[0].imshow(np.rot90(avg_vort[::-1,:], k=1), interpolation="none", cmap='RdBu_r', vmin=-maxval/2, vmax=maxval/2)
            axes[1].imshow(np.rot90(cor_vort[::-1,:], k=1), interpolation="none", cmap='RdBu_r', vmin=-maxval/2, vmax=maxval/2)
            axes[2].imshow(np.rot90(rnd_vort[::-1,:], k=1), interpolation="none", cmap='RdBu_r', vmin=-maxval/2, vmax=maxval/2)
            axes[3].imshow(np.rot90(vort[::-1,:], k=1), interpolation="none", cmap='RdBu_r', vmin=-maxval/2, vmax=maxval/2)
        else:
            axes[0,0].imshow(umean[::-1,:], cmap='RdBu_r', vmin=-maxval/2, vmax=maxval/2)
            axes[1,0].imshow(vmean[::-1,:], cmap='RdBu_r', vmin=-maxval/2, vmax=maxval/2)
            axes[2,0].imshow(0.5*(umean[::-1,:]**2+vmean[::-1,:]**2), cmap='RdBu_r', vmin=0.0, vmax=(maxval/4)**2)
            axes[0,1].imshow(u_cor[::-1,:], cmap='RdBu_r', vmin=-maxval/2, vmax=maxval/2)
            axes[1,1].imshow(v_cor[::-1,:], cmap='RdBu_r', vmin=-maxval/2, vmax=maxval/2)
            axes[2,1].imshow(0.5*(u_cor[::-1,:]**2+v_cor[::-1,:]**2), cmap='RdBu_r', vmin=0.0, vmax=(maxval/4)**2)
            axes[0,2].imshow(u_rnd[::-1,:], cmap='RdBu_r', vmin=-maxval/2, vmax=maxval/2)
            axes[1,2].imshow(v_rnd[::-1,:], cmap='RdBu_r', vmin=-maxval/2, vmax=maxval/2)
            axes[2,2].imshow(0.5*(u_rnd[::-1,:]**2+v_rnd[::-1,:]**2), cmap='RdBu_r', vmin=0.0, vmax=(maxval/4)**2)
            axes[0,3].imshow(u[::-1,:], cmap='RdBu_r', vmin=-maxval/2, vmax=maxval/2)
            axes[1,3].imshow(v[::-1,:], cmap='RdBu_r', vmin=-maxval/2, vmax=maxval/2)
            axes[2,3].imshow(0.5*(u[::-1,:]**2+v[::-1,:]**2), cmap='RdBu_r', vmin=0.0, vmax=(maxval/4)**2)
    
        plt.title(t/86400)
        for idx, ax in enumerate(axes_list):
            ax.tick_params(which='major', label1On=False)
            ax.tick_params(which='major', axis='both', width=0.25,  direction="in")
            ax.tick_params(which='minor', axis='both', width=0.125, direction='in')

            ax.tick_params(which='major', top=True, bottom=True, left=True, right=True)
            ax.tick_params(which='minor', top=True, bottom=True, left=True, right=True)

            for axis in ['top', 'bottom', 'left', 'right']:
                ax.spines[axis].set_linewidth(0.5)
        fig.canvas.draw()
        fig.canvas.flush_events()
        k = k + 1
	
"""
plt.savefig(OUTdir + 'GeostrophicVelocities.png',
            format='png', dpi='figure', bbox_inches='tight')

"""

