# Load modules
# ------------
import os
import sys
import glob
import pathlib
import numpy as np
from numpy import linalg as LA
import scipy.io as sio
import scipy.io.netcdf as nc
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from netCDF4 import Dataset

usrname = "ftucciar"
rootdir = "/home/" + usrname + "def_JAMES/"

infile = "/home/ftucciar/LU_JAMES/srcF90/Output/keLR_CsGrained_r3_WithMean_108d1_216d2.nc"
ingrid = "/home/ftucciar/dataNEMO/Deter/domain_cfg_R3.nc"
outfile1 = "/home/ftucciar/LU_JAMES/srcF90/Output/LR_aDMD_r3_WithMean_108d1_216d2.nc"
outfile2 = "/home/ftucciar/LU_JAMES/srcF90/Output/LR_diag_r3_WithMean_108d1_216d2.nc"


# infile = "/data/Stockage_12T_A/ftucciar/CsGrain/LR_CsGrained_r9_WithMean_6d1_24d2.nc"
# ingrid = "/home/ftucciar/dataNEMO/Deter/domain_cfg_R9.nc"
# outfile1 = "/data/Stockage_12T_A/ftucciar/CsGrain/LR_aDMD_r9_WithMean_6d1_24d2.nc"
# outfile2 = "/data/Stockage_12T_A/ftucciar/CsGrain/LR_diag_r9_WithMean_6d1_24d2.nc"



# Outputs
sfx = ['u', 'v', 'w']
drc = ['zonal', 'meridional', 'vertical']

dt = 5*86400

# Opening files
# -------------
fin = Dataset(infile, "r")
grid = Dataset(ingrid, "r")

# Get domain dimension from domain file
# -----------------------------------------------------------------------
nx = grid.variables['jpiglo'][:]
ny = grid.variables['jpjglo'][:]
nz = grid.variables['jpkglo'][:]

# Evaluation of the number of modes (i.e. of the times instants)
# -----------------------------------------------------------------------
nt = fin.variables["time_counter"][:].size
nt = 360




nm = nt - 1

pc = np.empty((nm, nz))

umean = np.empty((nz, ny, nx))
vmean = np.empty((nz, ny, nx))
wmean = np.empty((nz, ny, nx))

Uavg = np.empty((nz, ny, nx))
Vavg = np.empty((nz, ny, nx))
Wavg = np.empty((nz, ny, nx))

Uavg = fin.variables["umean"][:, :, :].data
Vavg = fin.variables["vmean"][:, :, :].data
Wavg = fin.variables["wmean"][:, :, :].data
Uavg[-1, :, :] = 0
Vavg[-1, :, :] = 0
Wavg[-1, :, :] = 0

umode = np.empty((nm, nz, ny, nx))
vmode = np.empty((nm, nz, ny, nx))
wmode = np.empty((nm, nz, ny, nx))

tmode = np.empty((nm, nt, nz))  # eigenvectors
tmeco = np.empty((nm, nz))      #

# %% Get Velocities
U = np.empty((0, nz, ny, nx))
V = np.empty((0, nz, ny, nx))
 
# Append velocity data from every .nc file
ids = "f1"
Uv = fin.variables["u" + ids][:nt, :, :, :].data
Vv = fin.variables["v" + ids][:nt, :, :, :].data
Wv = fin.variables["w" + ids][:nt, :, :, :].data
SSHv = fin.variables["ssh" + ids][:nt, :, :].data
Uv[:, -1, :, :] = 0
Vv[:, -1, :, :] = 0
Wv[:, -1, :, :] = 0

Um = np.mean(Uv, axis=0, dtype=np.float64)
Vm = np.mean(Vv, axis=0, dtype=np.float64)
Wm = np.mean(Wv, axis=0, dtype=np.float64)
Uv -= Um[None] 
Vv -= Vm[None] 
Wv -= Wm[None] 

# %% Create X1
X1 = np.concatenate((Uv.reshape((nt,nz*ny*nx)).T,
                     Vv.reshape((nt,nz*ny*nx)).T,
                     Wv.reshape((nt,nz*ny*nx)).T), axis=0)


sshX = SSHv.reshape((nt,ny*nx)).T

#X1 = np.concatenate((X1,X1), axis=0)

def DMD_modes( Xb, Xa, r, dt ):

    """    
       Name    : DMD
       Purpose : Performs the standard Dynamical Mode Decomposition

    Parameters
    ----------
            Xb : double,      First matrix of snapshots, "before".
            Xa : double,      Second matrix of snapshots, "after".
             r : double/int,  Desired rank of SVD
            dt : double,      Sampling step of the data

       Returns
       -------
       dmd_mds : double,      DMD modes
        ct_eig : double,      Continuous-time DMD eigenvalues
        dt_eig : double,      Discrete-time DMD eigenvalues
    """
    
    nd = np.size(Xb,0) # spatial dimension
    ns = np.size(Xb,1) # snapshots number
    
    # Compute SVD
    U, S, V = np.linalg.svd(Xb, full_matrices=False)

    # Dimensionality reduction parameter
    if isinstance(r, int):
        r = min(r, ns)     
    elif isinstance(r, float):
        ric = np.cumsum(S, dtype=np.float64)/np.sum(S, dtype=np.float64)
        r = np.abs(ric - r).argmin() + 1
    else:
        r = ns
    
    # Reduce dimensionality
    U = U[:,:r]
    S = S[:r]
    V = V[:,:r]
    
    # Low rank dynamics
    Atilde = np.matmul(np.matmul(np.matmul(U.transpose(), Xa), V), np.diag(1./S))

    # Discrete-time eigenvalues (\tilde{A}W=W\Lambda)
    dt_eig, W = np.linalg.eig(Atilde) 

    # DMD modes
    dmd_mds = np.matmul(np.matmul(np.matmul(Xa, V), np.diag(1./S)), W)

    return dmd_mds, dt_eig, U, S, V, W



def DMD_amplitudes( dmd_mds, X0 ):

    """    
       Name    : DMD_dynamics
       Purpose : Compute the dynamics encoded in the DMD modes

    Parameters
    ----------
       dmd_mds : double,  DMD modes
            X0 : double,  Initial snapshot

       Returns
       -------
           amp : double,      Vector of amplitudes of modes dmd_mds
    """    
    
    # DMD amplitudes (through Moore-Penrose inverse)
    amp = np.matmul(np.linalg.pinv(dmd_mds), X0) 

    return amp



def DMD_dynamics( dmd_mds, dt_eig, amp, dt, nt ):
    
    """    
       Name    : DMD_dynamics
       Purpose : Compute the dynamics encoded in the DMD modes

    Parameters
    ----------
       dmd_mds : double,  DMD modes
        dt_eig : double,  discrete-time DMD eigenvalues
           amp : double,  vector of amplitudes of modes dmd_mds
            dt : double,  Sampling step of the data
            nt : int,     Number of snapshots

       Returns
       -------
          recX : double,  data matrix reconstructed from dmd_mds, ct_eig, amp
    """    
    
    ns = np.size(dmd_mds,1) # modes number
    
    # Continuous-time eigenvalues
    ct_eig = np.log(dt_eig)/dt
    
    # Time vector
    t = np.linspace(0, nt, nt, endpoint=True) * dt
    
    # Time dynamics
    time_dynamics = np.zeros((ns, nt), dtype=complex)

    for iter in range(nt):
        time_dynamics[:,iter] = amp * np.exp( ct_eig * t[iter] )

    recX = np.matmul(dmd_mds, time_dynamics) 

    return recX.real, ct_eig


# POD




dmd_mds, dt_eig, U, S, V, W = DMD_modes( X1[:,:-1], X1[:,1:], 2000, dt )
amp = DMD_amplitudes( dmd_mds, X1[:,0] )
recX, ct_eig = DMD_dynamics( dmd_mds, dt_eig, amp, dt, np.size(X1,1) )



ssh_dmd_mds, ssh_dt_eig, ssh_U, ssh_S, ssh_V, ssh_W = DMD_modes( sshX[:,:-1], sshX[:,1:], 2000, dt )
ssh_amp = DMD_amplitudes( ssh_dmd_mds, sshX[:,0] )
ssh_recX, ssh_ct_eig = DMD_dynamics( ssh_dmd_mds, ssh_dt_eig, ssh_amp, dt, np.size(sshX,1) )

"""
U1, S1, V1 = LA.svd(X1[:,:-1], full_matrices=False)
plt.figure()
plt.plot(ric)
plt.show()

print(r'%d EOFs is used to capture %4.3f of total energy'%(idm,enpp))
# Low-rank
U1 = U1[:,:idm]
S1 = S1[:idm]
V1 = V1[:,:idm]
nm = len(S1)

# DMD
A1 = np.matmul(np.matmul(np.matmul(U1.transpose(), X1[:,1:]), V1), np.diag(1./S1))
lamb, W1 = LA.eig(A1) # 'lamb' is discrete-time eigenvalues
Phi = np.matmul(np.matmul(np.matmul(X1[:,1:], V1), np.diag(1./S1)), W1) # DMD modes
ct_lamb = np.log(lamb)/dt # continuous-time eigenvalues
b = np.matmul(LA.pinv(Phi), X1[:,0]) # DMD amptitude
"""


#ureid = (recX[:nx*ny*nz,:].real.transpose()).reshape((nt-1,nz,ny,nx))

#dmd_mds = dmd_mds[:2*nx*ny*nz,:]

print("done")




osc_idx = np.where(np.abs(np.abs(dt_eig) - 1) > 0.2 )
osc_mds = np.delete(dmd_mds, osc_idx, axis=1)
osc_amp = np.delete(amp, osc_idx)
osc_dt_eig = np.delete(dt_eig, osc_idx)
osc_ct_eig = np.delete(ct_eig, osc_idx)



grw_idx = np.where(np.abs(osc_dt_eig.imag) < 0.002 )
osc_mds = np.delete(osc_mds, grw_idx, axis=1)
osc_amp = np.delete(osc_amp, grw_idx)
osc_dt_eig = np.delete(osc_dt_eig, grw_idx)
osc_ct_eig = np.delete(osc_ct_eig, grw_idx)

nm = np.size(osc_dt_eig)
nt = np.size(osc_dt_eig)

recX, ct_eig = DMD_dynamics( osc_mds, osc_dt_eig, osc_amp, dt, nt )

"""


"""


ordered_idx = np.argsort(osc_ct_eig.real)


# %% Create ncfile
Phi = osc_mds[:,ordered_idx]*osc_amp[ordered_idx]
lamb = osc_dt_eig[ordered_idx]
ct_lamb = osc_ct_eig[ordered_idx]
b = osc_amp[ordered_idx]

ssh_osc_idx = np.where(np.abs(np.abs(ssh_dt_eig) - 1) > 0.2 )
ssh_osc_mds = np.delete(ssh_dmd_mds, ssh_osc_idx, axis=1)
ssh_osc_amp = np.delete(ssh_amp, ssh_osc_idx)
ssh_osc_dt_eig = np.delete(ssh_dt_eig, ssh_osc_idx)
ssh_osc_ct_eig = np.delete(ssh_ct_eig, ssh_osc_idx)


ssh_grw_idx = np.where(np.abs(ssh_osc_dt_eig.imag) < 0.002 )
ssh_osc_mds = np.delete(ssh_osc_mds, ssh_grw_idx, axis=1)
ssh_osc_amp = np.delete(ssh_osc_amp, ssh_grw_idx)
ssh_osc_dt_eig = np.delete(ssh_osc_dt_eig, ssh_grw_idx)
ssh_osc_ct_eig = np.delete(ssh_osc_ct_eig, ssh_grw_idx)

ssh_ordered_idx = np.argsort(ssh_osc_ct_eig.real)

ssh_Phi = ssh_osc_mds[:,ssh_ordered_idx]*ssh_osc_amp[ssh_ordered_idx]
ssh_lamb = ssh_osc_dt_eig[ssh_ordered_idx]
ssh_ct_lamb = ssh_osc_ct_eig[ssh_ordered_idx]
ssh_b = ssh_osc_amp[ssh_ordered_idx]

fout = Dataset(outfile1 , 'w', format='NETCDF4')

# Create dimensions
fout.createDimension('time', nt)
fout.createDimension('mode', nm)
fout.createDimension('y', ny)
fout.createDimension('x', nx)
fout.createDimension('z', nz)
# Create variables
tid = fout.createVariable('time', 'f8', ('time',))
mid = fout.createVariable('mode', 'i4', ('mode',))
yid = fout.createVariable('y', 'f8', ('y',))
xid = fout.createVariable('x', 'f8', ('x',))
zid = fout.createVariable('z', 'f8', ('z',))

lreid = fout.createVariable('lambda_real', 'f8', ('mode',))
limid = fout.createVariable('lambda_imag', 'f8', ('mode',))

ctlreid = fout.createVariable('ct_lamb_real', 'f8', ('mode',))
ctlimid = fout.createVariable('ct_lamb_imag', 'f8', ('mode',))

breid = fout.createVariable('b_real', 'f8', ('mode',))
bimid = fout.createVariable('b_imag', 'f8', ('mode',))

umid = fout.createVariable('umean', 'f8', ('z','y','x',))
vmid = fout.createVariable('vmean', 'f8', ('z','y','x',))
wmid = fout.createVariable('wmean', 'f8', ('z','y','x',))

ureid = fout.createVariable('umode_real', 'f8', ('mode','z','y','x',))
uimid = fout.createVariable('umode_imag', 'f8', ('mode','z','y','x',))

vreid = fout.createVariable('vmode_real', 'f8', ('mode','z','y','x',))
vimid = fout.createVariable('vmode_imag', 'f8', ('mode','z','y','x',))

wreid = fout.createVariable('wmode_real', 'f8', ('mode','z','y','x',))
wimid = fout.createVariable('wmode_imag', 'f8', ('mode','z','y','x',))

ssh_reid = fout.createVariable('ssh_mode_real', 'f8', ('mode','y','x',))
ssh_imid = fout.createVariable('ssh_mode_imag', 'f8', ('mode','y','x',))

ssh_lreid = fout.createVariable('ssh_lambda_real', 'f8', ('mode',))
ssh_limid = fout.createVariable('ssh_lambda_imag', 'f8', ('mode',))

ssh_ctlreid = fout.createVariable('ssh_ct_lamb_real', 'f8', ('mode',))
ssh_ctlimid = fout.createVariable('ssh_ct_lamb_imag', 'f8', ('mode',))

ssh_breid = fout.createVariable('ssh_b_real', 'f8', ('mode',))
ssh_bimid = fout.createVariable('ssh_b_imag', 'f8', ('mode',))

# Add attributes
tid.long_name = 'Time axis'
tid.units = 'years'
mid.long_name = 'Mode index'
yid.long_name = 'Ocean Y axis (T-grid)'
yid.units = 'km'
xid.long_name = 'Ocean X axis (T-grid)'
xid.units = 'km'
zid.long_name = 'Ocean mid-layer axis'
zid.units = 'km'
umid.long_name = 'Mean of zonal velocity'
umid.units = 'm/s'
vmid.long_name = 'Mean of meridional velocity'
vmid.units = 'm/s'
wmid.long_name = 'Mean of vertical velocity'
wmid.units = 'm/s'
lreid.long_name = 'DMD eigenvalues (real part)'
limid.long_name = 'DMD eigenvalues (imag part)'
ctlreid.long_name = 'Continuous time eigenvalues (real part)'
ctlimid.long_name = 'Continuous time eigenvalues (imag part)'
breid.long_name = 'DMD amplitudes (real part)'
bimid.long_name = 'DMD amplitudes (imag part)'
ureid.long_name = 'Zonal DMD modes (real part)'
vreid.long_name = 'Meridional DMD modes (real part)'
wreid.long_name = 'Vertical DMD modes (real part)'
uimid.long_name = 'Zonal DMD modes (imag part)'
vimid.long_name = 'Meridional DMD modes (imag part)'
wimid.long_name = 'Vertical DMD modes (imag part)'


ssh_reid.long_name = 'SSH DMD modes (real part)'
ssh_imid.long_name = 'SSH DMD modes (imag part)'

# Write data
tid[:] = range(nt)
mid[:] = range(nm)
yid[:] = range(ny)
xid[:] = range(nx)
zid[:] = range(nz)
lreid[:] = lamb.real
limid[:] = lamb.imag
ctlreid[:] = ct_lamb.real
ctlimid[:] = ct_lamb.imag
breid[:] = b.real
bimid[:] = b.imag
umid[:,:,:] = Uavg
vmid[:,:,:] = Vavg
wmid[:,:,:] = Wavg
ureid[:,:,:,:] = (Phi[:nx*ny*nz,:].real.transpose()).reshape((nm,nz,ny,nx))
vreid[:,:,:,:] = (Phi[nx*ny*nz:2*nx*ny*nz,:].real.transpose()).reshape((nm,nz,ny,nx))
wreid[:,:,:,:] = (Phi[2*nx*ny*nz:,:].real.transpose()).reshape((nm,nz,ny,nx))
uimid[:,:,:,:] = (Phi[:nx*ny*nz,:].imag.transpose()).reshape((nm,nz,ny,nx))
vimid[:,:,:,:] = (Phi[nx*ny*nz:2*nx*ny*nz,:].imag.transpose()).reshape((nm,nz,ny,nx))
wimid[:,:,:,:] = (Phi[2*nx*ny*nz:,:].imag.transpose()).reshape((nm,nz,ny,nx))


ssh_lreid[:] = ssh_lamb.real
ssh_limid[:] = ssh_lamb.imag
ssh_ctlreid[:] = ssh_ct_lamb.real
ssh_ctlimid[:] = ssh_ct_lamb.imag
ssh_breid[:] = ssh_b.real
ssh_bimid[:] = ssh_b.imag
ssh_reid[:,:,:] = (ssh_Phi[:nx*ny,:].real.transpose()).reshape((nm,ny,nx))
ssh_imid[:,:,:] = (ssh_Phi[:nx*ny,:].imag.transpose()).reshape((nm,ny,nx))
# Close file
fout.close()

print("done")
# Reshape outputs

u_mdre = (Phi[:nx*ny*nz,:].real.transpose()).reshape((nm,nz,ny,nx))
u_mdim = (Phi[:nx*ny*nz,:].imag.transpose()).reshape((nm,nz,ny,nx))
v_mdre = (Phi[nx*ny*nz:2*nx*ny*nz,:].real.transpose()).reshape((nm,nz,ny,nx))
v_mdim = (Phi[nx*ny*nz:2*nx*ny*nz,:].imag.transpose()).reshape((nm,nz,ny,nx))

u_mds = (Phi[:nx*ny*nz,:].transpose()).reshape((nm,nz,ny,nx))
v_mds = (Phi[nx*ny*nz:2*nx*ny*nz,:].transpose()).reshape((nm,nz,ny,nx))
w_mds = (Phi[2*nx*ny*nz:3*nx*ny*nz,:].transpose()).reshape((nm,nz,ny,nx))

# Save outputs

nm_c = 100
nm_r = nm - nm_c

print(nm, nm_c, nm_r)
print(np.shape(lamb.imag),np.shape(lamb[nm_c:].imag),np.shape(lamb[:nm_c].imag))

omg_c = lamb[:nm_c].imag
omg_nc = lamb[nm_c:].imag

uc = Um
vc = Vm
wc = Wm
um_c = u_mds[:nm_c,...]
vm_c = v_mds[:nm_c,...]
wm_c = w_mds[:nm_c,...]
um_nc = u_mds[nm_c:,...]
vm_nc = v_mds[nm_c:,...]
wm_nc = w_mds[nm_c:,...]


print(np.shape(omg_c),np.shape(omg_nc))
print(np.shape(um_c),np.shape(vm_c),np.shape(um_nc),np.shape(vm_nc))







fout = Dataset(outfile2 , 'w', format='NETCDF4')

fout.createDimension('mode_c', nm_c)
fout.createDimension('mode_r', nm_r)

fout.createDimension('y', ny)
fout.createDimension('x', nx)
fout.createDimension('z', nz)

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
mrid[:] = range(nm_r)

yid[:] = range(ny)
xid[:] = range(nx)
zid[:] = range(nz)


wcid[:] = omg_c
wrid[:] = omg_nc
umid[...] = uc
vmid[...] = vc

ucreid[...] = um_c.real
vcreid[...] = vm_c.real

ucimid[...] = um_c.imag
vcimid[...] = vm_c.imag

urreid[...] = um_nc.real
vrreid[...] = vm_nc.real
wrreid[...] = wm_nc.real

urimid[...] = um_nc.imag
vrimid[...] = vm_nc.imag
wrimid[...] = wm_nc.imag

# Close file
fout.close()

print("")
print('Program terminates')

























