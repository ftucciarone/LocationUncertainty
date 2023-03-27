#!/usr/bin/env python

# Load modules
# ------------
import sys
import pathlib
import numpy as np
from numpy import linalg as LA
import scipy.io.netcdf as nc
from netCDF4 import Dataset
import scipy.io as sio
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.ndimage import gaussian_filter
import os
plt.ion()
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
infile = data_dir + "/" + prefix + "_DMD_r9_" + suffix + ".nc"
ingrid = dom_dir + "/domain_cfg_R9.nc"

print("isfile=", os.path.isfile(ingrid), "          Opening domain file:", ingrid)
print("isfile=", os.path.isfile(infile), "Opening coarse-grained fields:", infile)
print("")

# %% LR_DMD_r3_ZeroMean_36w1_72w2
outfile = prefix + "_ocludat_r9_" + suffix + ".nc"


# Set parammeters
# ---------------
# Inputs
home_dir = str(pathlib.Path.home())


init_file =  home_dir + "/__future_LocationUncertainty/data/CsGrain_Output/" + prefix + "_CsGrained_r9_" + suffix + ".nc"

filt_wd=1
# Cutloffs of physical times
freq_time = np.zeros((10))
freq_time[0] = np.pi/(5.*86400.) # cutoff correlated/uncorrelated (s)
freq_time[1] = np.pi/(10.*86400.) # cutoff correlated/uncorrelated (s) 
freq_time[2] = np.pi/(15.*86400.) # cutoff correlated/uncorrelated (s) 
freq_time[3] = np.pi/(20.*86400.) # cutoff correlated/uncorrelated (s) 
freq_time[4] = np.pi/(25.*86400.) # cutoff correlated/uncorrelated (s) 
freq_time[5] = np.pi/(30.*86400.) # cutoff correlated/uncorrelated (s) 
freq_time[6] = np.pi/(45.*86400.) # cutoff correlated/uncorrelated (s)
freq_time[7] = np.pi/(60.*86400.) # cutoff correlated/uncorrelated (s)
freq_time[8] = np.pi/(90.*86400.) # cutoff correlated/uncorrelated (s)
freq_time[9] = np.pi/(120.*86400.) # cutoff correlated/uncorrelated (s)
freq_time_label = ["  5 days", " 10 days", " 15 days", " 20 days", \
                   " 25 days", " 30 days", " 45 days", " 60 days", \
                   " 90 days", "120 days", "120 days" ] 
radius = 1.1

# Read DMD data
f = Dataset(infile,'r')
x = f.variables['x'][:].copy() 
y = f.variables['y'][:].copy()
z = f.variables['z'][:].copy()
tyrs = f.variables['time'][:].copy()

re = f.variables['lambda_real'][:].copy() 
im = f.variables['lambda_imag'][:].copy()
lam = re + 1j*im

re = f.variables['ct_lamb_real'][:].copy() 
im = f.variables['ct_lamb_imag'][:].copy() 
ct_lam = re + 1j*im

re = f.variables['b_real'][:].copy() 
im = f.variables['b_imag'][:].copy() 
amp = re + 1j*im

uc = f.variables['umean'][...].copy()
vc = f.variables['vmean'][...].copy()
wc = f.variables['wmean'][...].copy()

re = f.variables['umode_real'][...].copy() 
im = f.variables['umode_imag'][...].copy() 
um = re + 1j*im

re = f.variables['vmode_real'][...].copy() 
im = f.variables['vmode_imag'][...].copy() 
vm = re + 1j*im

re = f.variables['wmode_real'][...].copy() 
im = f.variables['wmode_imag'][...].copy() 
wm = re + 1j*im


del re, im
f.close()
dt = (tyrs[1]-tyrs[0])*360.*86400.
nx = len(x)
ny = len(y)
nl = len(z)

# Plot (unfiltered) eigenvalues and amptitudes 
plt.figure(figsize=(5,5))
plt.plot(lam.real, lam.imag, '.')
plt.grid(which='both', axis='both')
plt.xlabel(r'Re ($\lambda$)')
plt.ylabel(r'Im ($\lambda$)')
plt.savefig(out_dir + "/" + prefix + suffix + '_circ.eps', dpi=200, bbox_inches='tight', pad_inches=0)

amp_mod = np.abs(amp)
plt.figure(figsize=(5,5))
plt.plot(1.0e5*ct_lam.imag, amp_mod, '+')
plt.grid(which='both', axis='both')
plt.xlabel(r'Im ($\omega$) $\times 10^{-5}$')
plt.ylabel(r'amplitude')
plt.savefig(out_dir + "/" + prefix + suffix + '_amp.eps', dpi=200, bbox_inches='tight', pad_inches=0)

# Filtering eigenvalues (kill those not on unit cercle)
lam_mod = np.abs(lam)
idm = np.where((lam_mod > 1.0 - 1.0e-2) & (lam_mod <= 1.0 + 1.0e-3))
mask = np.full(len(lam_mod), False, dtype=bool)
mask[idm] = True
lam = lam[mask]
ct_lam = ct_lam[mask]
amp = amp[mask]
amp_mod = amp_mod[mask]
um = um[mask,...]
vm = vm[mask,...]
wm = wm[mask,...]
del lam_mod

plt.figure(figsize=(5,5))
plt.plot(lam.real, lam.imag, '.')
plt.grid(which='both', axis='both')
plt.xlabel(r'Re ($\lambda$)')
plt.ylabel(r'Im ($\lambda$)')
plt.savefig(out_dir + "/" + prefix + suffix + '_circ_filt.eps', dpi=200, bbox_inches='tight', pad_inches=0)

plt.figure(figsize=(5,5))
plt.plot(1.0e5*ct_lam.imag, amp_mod, '+')
plt.grid(which='both', axis='both')
plt.xlabel(r'Im ($\omega$) $\times 10^{-5}$')
plt.ylabel(r'amplitude')
plt.savefig(out_dir + "/" + prefix + suffix + '_amp_filt.eps', dpi=200, bbox_inches='tight', pad_inches=0)

# Seperate temporal scales
plt.figure(figsize=(5,5))
for i in range(len(freq_time)):
    idm_temp = np.where(np.abs(ct_lam.imag) <= freq_time[i])
    mask_temp = np.full(len(ct_lam.imag), False, dtype=bool)
    mask_temp[idm_temp] = True
    plt.plot(lam.real[mask_temp], lam.imag[mask_temp], '.',label=freq_time_label[i])
    del mask_temp, idm_temp
plt.grid(which='both', axis='both')
plt.legend(loc='best')
plt.xlabel(r'Re ($\lambda$)')
plt.ylabel(r'Im ($\lambda$)')
plt.savefig(out_dir + "/" + prefix + suffix + '_circ_temp.eps', dpi=200, bbox_inches='tight', pad_inches=0)


input_days = float(input("Correlated/Random treshold:\n"))
freq_cut = np.pi/(input_days*86400.) # cutoff correlated/uncorrelated (s) 
idm = np.where(np.abs(ct_lam.imag) <= freq_cut)
mask = np.full(len(ct_lam.imag), False, dtype=bool)
mask[idm] = True
um_c = um[mask,...]
vm_c = vm[mask,...]
wm_c = wm[mask,...]
um_nc = um[~mask,...]
vm_nc = vm[~mask,...]
wm_nc = wm[~mask,...]
del um, vm, wm

lam_c = lam[mask]
lam_nc = lam[~mask]
del lam

omg_c = ct_lam[mask]
omg_nc = ct_lam[~mask]
del ct_lam

amp_c = amp[mask]
amp_mod_c = amp_mod[mask]
amp_nc = amp[~mask]
amp_mod_nc = amp_mod[~mask]
del amp, amp_mod

plt.figure(figsize=(5,5))
plt.plot(lam_c.real, lam_c.imag, '.', label=r'correlated')
plt.plot(lam_nc.real, lam_nc.imag, '.', label='uncorrelated')
plt.grid(which='both', axis='both')
plt.xlabel(r'Re ($\lambda$)')
plt.ylabel(r'Im ($\lambda$)')
plt.legend(loc='best')
plt.savefig(out_dir + "/" + prefix + suffix + '_circ_filt_sep.eps', dpi=200, bbox_inches='tight', pad_inches=0)

plt.figure(figsize=(5,5))
plt.plot(1.0e5*omg_c.imag, amp_mod_c, '+', label=r'correlated')
plt.plot(1.0e5*omg_nc.imag, amp_mod_nc, '+', label=r'uncorrelated')
plt.grid(which='both', axis='both')
plt.xlabel(r'Im ($\omega$) $\times 10^{-5}$')
plt.ylabel(r'amplitude')
plt.legend(loc='best')
plt.savefig(out_dir + "/" + prefix + suffix + '_amp_filt_sep.eps', dpi=200, bbox_inches='tight', pad_inches=0)

# Filtering modes by energy
amp_max = float(input("Set trucation amplitude for CORRELATED modes:\n"))
idm = np.where( (amp_mod_c >= amp_max) & (np.abs(omg_c.imag) > 2.5e-7) )
mask = np.full(len(amp_mod_c), False, dtype=bool)
mask[idm] = True
lam_c = lam_c[mask]
omg_c = omg_c[mask]
amp_c = amp_c[mask]
amp_mod_c = amp_mod_c[mask]
um_c = um_c[mask]
vm_c = vm_c[mask]
wm_c = wm_c[mask]
nm_c = len(amp_mod_c)
print('Number of modes used for correlated drift = ',nm_c)

amp_max = float(input("Set trucation amplitude for RANDOM modes:\n"))
idm = np.where( (amp_mod_nc >= amp_max) & (np.abs(omg_nc.imag) < 7.e-6) )
mask = np.full(len(amp_mod_nc), False, dtype=bool)
mask[idm] = True
lam_nc = lam_nc[mask]
omg_nc = omg_nc[mask]
amp_nc = amp_nc[mask]
amp_mod_nc = amp_mod_nc[mask]
um_nc = um_nc[mask]
vm_nc = vm_nc[mask]
wm_nc = wm_nc[mask]
nm_nc = len(amp_mod_nc)
print('Number of modes used for uncorrelated noise = ',nm_nc)

plt.figure(figsize=(5,5))
plt.plot(lam_c.real, lam_c.imag, '.', label=r'correlated')
plt.plot(lam_nc.real, lam_nc.imag, '.', label='uncorrelated')
plt.xlim([-1.1,1.1])
plt.ylim([-1.1,1.1])
plt.grid(which='both', axis='both')
plt.xlabel(r'Re ($\lambda$)')
plt.ylabel(r'Im ($\lambda$)')
plt.legend(loc='best')
plt.savefig(out_dir + "/" + prefix + suffix + '_circ_filt_sep_cut.eps', dpi=200, bbox_inches='tight', pad_inches=0)

plt.figure(figsize=(5,5))
plt.plot(1.0e5*omg_c.imag, amp_mod_c, '+', label=r'correlated')
plt.plot(1.0e5*omg_nc.imag, amp_mod_nc, '+', label=r'uncorrelated')
plt.xlim([-0.8,0.8])
plt.grid(which='both', axis='both')
plt.xlabel(r'Im ($\omega$) $\times 10^{-5}$')
plt.ylabel(r'amplitude')
plt.legend(loc='best')
plt.savefig(out_dir + "/" + prefix + suffix + '_amp_filt_sep_cut.eps', dpi=200, bbox_inches='tight', pad_inches=0)


# Gramian procedure
# Rescale truncated amplitudes
fin = Dataset(init_file,'r')
u = fin.variables['uband'][0,...].copy() 
v = fin.variables['vband'][0,...].copy() 
w = fin.variables['wband'][0,...].copy() 
fin.close()

u -= uc  # This depends on "if i have already removed the mean or not from my coarse grained velocity"
v -= vc
w -= wc

X0 = u.flatten() 
X0 = np.append(X0, v.flatten())
X0 = np.append(X0, w.flatten())
del u, v, w

nu = nl*nx*ny
nv = nl*nx*ny
nw = nl*nx*ny

re = np.zeros((nu+nv+nw,nm_c), dtype=np.float64)
im = np.zeros((nu+nv+nw,nm_c), dtype=np.float64)
re[:nu,:] = um_c.real.reshape((nm_c,nu)).transpose() 
im[:nu,:] = um_c.imag.reshape((nm_c,nu)).transpose()  
re[nu:nu+nv,:] = vm_c.real.reshape((nm_c,nv)).transpose() 
im[nu:nu+nv,:] = vm_c.imag.reshape((nm_c,nv)).transpose() 
re[nu+nv:,:] = wm_c.real.reshape((nm_c,nw)).transpose()
im[nu+nv:,:] = wm_c.imag.reshape((nm_c,nw)).transpose() 
Phi = re + 1j*im
del re, im

Psi = np.matmul( Phi, LA.pinv(np.matmul(Phi.conj().transpose(), Phi)) )
amp_c = np.matmul(Psi.conj().transpose(), X0)
amp_mod_c = np.abs(amp_c)

re = np.zeros((nu+nv+nw,nm_nc), dtype=np.float64)
im = np.zeros((nu+nv+nw,nm_nc), dtype=np.float64)
re[:nu,:] = um_nc.real.reshape((nm_nc,nu)).transpose()  
im[:nu,:] = um_nc.imag.reshape((nm_nc,nu)).transpose() 
re[nu:nu+nv,:] = vm_nc.real.reshape((nm_nc,nv)).transpose() 
im[nu:nu+nv,:] = vm_nc.imag.reshape((nm_nc,nv)).transpose()
re[nu+nv:,:] = wm_nc.real.reshape((nm_nc,nw)).transpose() 
im[nu+nv:,:] = wm_nc.imag.reshape((nm_nc,nw)).transpose() 
Phi = re + 1j*im
del re, im

Psi = np.matmul( Phi, LA.pinv(np.matmul(Phi.conj().transpose(), Phi)) )
amp_nc = np.matmul(Psi.conj().transpose(), X0)
amp_mod_nc = np.abs(amp_nc)
del Phi, Psi, X0

plt.figure(figsize=(5,5))
plt.plot(1.0e5*omg_c.imag, amp_mod_c, '+', label=r'correlated')
plt.plot(1.0e5*omg_nc.imag, amp_mod_nc, '+', label=r'uncorrelated')
plt.xlim([-0.8,0.8])
plt.grid(which='both', axis='both')
plt.xlabel(r'Im ($\omega$) $\times 10^{-5}$')
plt.ylabel(r'amplitude')
plt.legend(loc='best')
plt.savefig(out_dir + "/" + prefix + suffix + '_amp_filt_sep_cut_res.eps', dpi=200, bbox_inches='tight', pad_inches=0)

# Sorting the arrays
omg_c_idx = np.abs(omg_c.imag).argsort()
omg_c = omg_c[omg_c_idx[::-1]]
amp_c = amp_c[omg_c_idx[::-1]]
um_c = um_c[omg_c_idx[::-1],:,:,:]
vm_c = vm_c[omg_c_idx[::-1],:,:,:]
wm_c = wm_c[omg_c_idx[::-1],:,:,:]

omg_nc_idx = np.abs(omg_nc.imag).argsort()
omg_nc = omg_nc[omg_nc_idx[::-1]]
amp_nc = amp_nc[omg_nc_idx[::-1]]
um_nc = um_nc[omg_nc_idx[::-1],:,:,:]
vm_nc = vm_nc[omg_nc_idx[::-1],:,:,:]
wm_nc = wm_nc[omg_nc_idx[::-1],:,:,:]


# Rescale DMD modes
for i in range(nm_c):
    um_c[i,...] *= amp_c[i]
    vm_c[i,...] *= amp_c[i]
    wm_c[i,...] *= amp_c[i]

for i in range(nm_nc):
    um_nc[i,...] *= amp_nc[i]
    vm_nc[i,...] *= amp_nc[i]
    wm_nc[i,...] *= amp_nc[i]

# Save outputs
fout = Dataset(out_dir + "/" + outfile, 'w', format='NETCDF4')

fout.createDimension('mode_c', nm_c)
fout.createDimension('mode_r', nm_nc)

fout.createDimension('y', ny)
fout.createDimension('x', nx)
fout.createDimension('z', nl)

yid = fout.createVariable('y', 'f8', ('y',))
xid = fout.createVariable('x', 'f8', ('x',))
zid = fout.createVariable('z', 'f8', ('z',))

mcid = fout.createVariable('mode_c', 'i4', ('mode_c',))
mrid = fout.createVariable('mode_r', 'i4', ('mode_r',))

umid = fout.createVariable('uco', 'f8', ('z','y','x',))
vmid = fout.createVariable('vco', 'f8', ('z','y','x',))
wmid = fout.createVariable('wco', 'f8', ('z','y','x',))

wcid = fout.createVariable('omega_c', 'f8', ('mode_c',))
wrid = fout.createVariable('omega_r', 'f8', ('mode_r',))

ucreid = fout.createVariable('umode_real_c', 'f8', ('mode_c','z','y','x',))
ucimid = fout.createVariable('umode_imag_c', 'f8', ('mode_c','z','y','x',))
vcreid = fout.createVariable('vmode_real_c', 'f8', ('mode_c','z','y','x',))
vcimid = fout.createVariable('vmode_imag_c', 'f8', ('mode_c','z','y','x',))
wcreid = fout.createVariable('wmode_real_c', 'f8', ('mode_c','z','y','x',))
wcimid = fout.createVariable('wmode_imag_c', 'f8', ('mode_c','z','y','x',))

urreid = fout.createVariable('umode_real_r', 'f8', ('mode_r','z','y','x',))
urimid = fout.createVariable('umode_imag_r', 'f8', ('mode_r','z','y','x',))
vrreid = fout.createVariable('vmode_real_r', 'f8', ('mode_r','z','y','x',))
vrimid = fout.createVariable('vmode_imag_r', 'f8', ('mode_r','z','y','x',))
wrreid = fout.createVariable('wmode_real_r', 'f8', ('mode_r','z','y','x',))
wrimid = fout.createVariable('wmode_imag_r', 'f8', ('mode_r','z','y','x',))

# Add attributes
mcid.long_name = 'Mode index (for correction drift)'
mrid.long_name = 'Mode index (for random noise)'

xid.long_name = 'Ocean X axis'
yid.long_name = 'Ocean Y axis'
zid.long_name = 'Ocean Z axis'
xid.units = 'km'
yid.units = 'km'
zid.units = 'km'

wcid.long_name = 'Continuous eigenvalues for correction drift'
wrid.long_name = 'Continuous eigenvalues for random noise'
wcid.units = 's^-1'
wrid.units = 's^-1'

umid.long_name = 'Mean of zonal correction drift'
vmid.long_name = 'Mean of meridional correction drift'
wmid.long_name = 'Mean of vertical correction drift'
umid.units = 'm/s'
vmid.units = 'm/s'
wmid.units = 'm/s'

ucreid.long_name = 'Zonal DMD modes for correction drift (real part)'
vcreid.long_name = 'Meridional DMD modes for correction drift (real part)'
wcreid.long_name = 'Vertical DMD modes for correction drift (real part)'
ucreid.units = 'm/s'
vcreid.units = 'm/s'
wcreid.units = 'm/s'

ucimid.long_name = 'Zonal DMD modes for correction drift (imag part)'
vcimid.long_name = 'Meridional DMD modes for correction drift (imag part)'
wcimid.long_name = 'Vertical DMD modes for correction drift (imag part)'
ucimid.units = 'm/s'
vcimid.units = 'm/s'
wcimid.units = 'm/s'

urreid.long_name = 'Zonal DMD modes for random noise (real part)'
vrreid.long_name = 'Meridional DMD modes for random noise (real part)'
wrreid.long_name = 'Vertical DMD modes for random noise (real part)'
urreid.units = 'm/s'
vrreid.units = 'm/s'
wrreid.units = 'm/s'

urimid.long_name = 'Zonal DMD modes for random noise (imag part)'
vrimid.long_name = 'Meridional DMD modes for random noise (imag part)'
wrimid.long_name = 'Vertical DMD modes for random noise (imag part)'
vrimid.units = 'm/s'
urimid.units = 'm/s'
wrimid.units = 'm/s'

# Write data
mcid[:] = range(nm_c)
mrid[:] = range(nm_nc)

yid[:] = y
xid[:] = x
zid[:] = z

wcid[:] = omg_c.imag
wrid[:] = omg_nc.imag

umid[...] = gaussian_filter(uc, sigma=filt_wd)
vmid[...] = gaussian_filter(vc, sigma=filt_wd)
wmid[...] = gaussian_filter(wc, sigma=filt_wd)

ucreid[...] = gaussian_filter(um_c.real, sigma=filt_wd)
vcreid[...] = gaussian_filter(vm_c.real, sigma=filt_wd)
wcreid[...] = gaussian_filter(wm_c.real, sigma=filt_wd)

ucimid[...] = gaussian_filter(um_c.imag, sigma=filt_wd)
vcimid[...] = gaussian_filter(vm_c.imag, sigma=filt_wd)
wcimid[...] = gaussian_filter(wm_c.imag, sigma=filt_wd)

urreid[...] = gaussian_filter(um_nc.real, sigma=filt_wd)
vrreid[...] = gaussian_filter(vm_nc.real, sigma=filt_wd)
wrreid[...] = gaussian_filter(wm_nc.real, sigma=filt_wd)

urimid[...] = gaussian_filter(um_nc.imag, sigma=filt_wd)
vrimid[...] = gaussian_filter(vm_nc.imag, sigma=filt_wd)
wrimid[...] = gaussian_filter(wm_nc.imag, sigma=filt_wd)

# Close file
fout.close()

print("")
print('Program terminates')


