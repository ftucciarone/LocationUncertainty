#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 23 16:55:00 2023


@author: ftucciar
"""
import copy
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.colors as colors

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


def DMD_continuous_eig( dt_eig, dt ):
    """    
       Name    : DMD_continuous_eig
       Purpose : Compute the DMD continuous eigenvalues

    Parameters
    ----------
        dt_eig : double,  discrete-time DMD eigenvalues
            dt : double,  Sampling step of the data

       Returns
       -------
        ct_eig : double,  DMD continuous eigenvalues
    """    
    
    # Continuous-time eigenvalues
    ct_eig = np.log(dt_eig)/dt
    
    return ct_eig


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

    return recX.real


def DefineXandXp( S ):
    """    
       Name    : DefineXandXp
       Purpose : Defines X and Xp from the snapshot matrix

    Parameters
    ----------
             S : double,  Matrix of snapshots

       Returns
       -------
             X : double,  First matrix of snapshots
            Xp : double,  Second matrix of snapshots
    """   
    
    X = S[:,:-1]
    Xp = S[:,1:]
    
    return  X,  Xp
    

def ShiftStack( S, m, s ):
    """    
       Name    : ShiftStack
       Purpose : Stacks 's' times the snapshot matrix Snap

    Parameters
    ----------
             S : double,  Matrix of snapshots
             m : int,     Number of snapshots
             s : int,     Number of stacks performed

       Returns
       -------
           ssS : double,  DMD modes
    """    

    ssS = S[...,:m-s]
    for i in range(s):
        ssS = np.concatenate((ssS, S[...,i+1:m-s+i+1]), axis=0)

    return ssS


def CleanUnitCircle( dmd_mds, ct_eig, dt_eig, amp, tol ):    
    """    
       Name    : CleanUnitCircle
       Purpose : Removes the modes too far from the unit circle 

    Parameters
    ----------
       dmd_mds : double,  DMD modes
        ct_eig : double,  continuous-time DMD eigenvalues
        dt_eig : double,  discrete-time DMD eigenvalues
           amp : double,  vector of amplitudes of modes dmd_mds
           tol : double,  tolerance to drop 

    """   
    
    non_osc_idx = np.where(np.abs(np.abs(dt_eig) - 1) > tol )

    circ_amp = np.delete(amp, non_osc_idx)
    circ_mds = np.delete(dmd_mds, non_osc_idx, axis=1)
    circ_ct_eig = np.delete(ct_eig, non_osc_idx)
    circ_dt_eig = np.delete(dt_eig, non_osc_idx)

    return circ_mds, circ_ct_eig, circ_dt_eig, circ_amp
 

def PlotAmplitudes(ct_eig, amp, title=None, split=None, freq_time=None, saveout=None, filename=None):
    """    
       Name    : PlotAmplitudes
       Purpose : Plot the amplitudes and saves the plot with a given path/filename

    Parameters
    ----------
        ct_eig : double,    continuous-time DMD eigenvalues
           amp : double,    vector of amplitudes of modes dmd_mds
         title : character, optional title 
       saveout : logical,   flag to enable saving
      filename : character, filename and path to save the plot

    """   
    plt.figure(figsize=(5,5))

    if split is not None:
        tresh = 1.0e5*np.pi/(split*86400.)
        plt.axvline(x=tresh, color="orange" )
        plt.axvline(x=-tresh, color="orange" )
    plt.plot(1.0e5*ct_eig.imag, np.abs(amp), '+')
     
    if title is not None:
        plt.title(r'Amplidute distribution '+ title)
    else:
        plt.title(r'Amplidute distribution')
    plt.grid(which='both', axis='both')
    plt.grid(which='both', axis='both')
 
    if freq_time is None:
        plt.xlabel(r'Im ($\omega$) $\times 10^{-5}$')
    else:   
        freqs =  np.append(freq_time, -freq_time[::-1])
        ticks = 1.0e5*np.pi / (freqs*86400.)
        labels = [ str(abs(x)) for x in freqs ]
        plt.xticks(ticks, labels, rotation=90)
        plt.xlabel(r'Im ($\omega$) (days)')
    plt.ylabel(r'Amplitude')
    plt.ylim([0, 60])
    if saveout is not None: plt.savefig(filename, dpi=200, bbox_inches='tight', pad_inches=0)


def PlotCircle(ct_eig, dt_eig, title=None, split=None, freq_time=None, saveout=None, filename=None):
    """    
       Name    : PlotCircle
       Purpose : Plot the eigencircle and saves the plot with a given path/filename

    Parameters
    ----------
        ct_eig : double,    continuous-time DMD eigenvalues
        dt_eig : double,    discrete-time DMD eigenvalues
         title : character, optional title 
     freq_time : double,    frequency (in days) array for scale separation
       saveout : logical,   flag to enable saving
      filename : character, filename and path to save the plot

    """   
    plt.figure(figsize=(5,5))
    if split is not None:
        id_split = np.argmin(np.abs(np.abs(ct_eig.imag) - np.pi/(split*86400.)))
        plt.plot([0, dt_eig[id_split].real], [0, dt_eig[id_split].imag], color="orange" )
        plt.plot([0, dt_eig[id_split].real], [0, -dt_eig[id_split].imag], color="orange")

    if freq_time is None: 
        plt.plot(dt_eig.real, dt_eig.imag, ".")
    else:
        for i in range(len(freq_time)):
            idm_temp = np.where(np.abs(ct_eig.imag) <= np.pi/(freq_time[i] * 86400.))
            mask_temp = np.full(len(ct_eig.imag), False, dtype=bool)
            mask_temp[idm_temp] = True
            plt.plot(dt_eig.real[mask_temp], dt_eig.imag[mask_temp], '.',label=str(freq_time[i]) + " days")
            del mask_temp, idm_temp
            plt.legend(loc='best')
    plt.grid(which='both', axis='both')

    if title is not None:
        plt.title(r'Eigenvalue $\lambda$ disposition '+ title)
    else:
        plt.title(r'Eigenvalue $\lambda$ disposition')

    plt.xlim([-1.05, 1.05])
    plt.ylim([-1.05, 1.05])
    plt.xlabel(r'Re ($\lambda$)')
    plt.ylabel(r'Im ($\lambda$)')
    if saveout is not None: plt.savefig(filename, dpi=200, bbox_inches='tight', pad_inches=0)
        

def CleanExponentialModes( dmd_mds, ct_eig, dt_eig, amp, tol ):    
    """    
       Name    : CleanExponentialModes
       Purpose : Removes modes with small Imaginary part 

    Parameters
    ----------
       dmd_mds : double,  DMD modes
        ct_eig : double,  continuous-time DMD eigenvalues
        dt_eig : double,  discrete-time DMD eigenvalues
           amp : double,  vector of amplitudes of modes dmd_mds
           tol : double,  tolerance to drop type('CopyOfB', self.__bases__, dict(self.__dict__))

    """   

    exp_idx = np.where(np.abs(dt_eig.imag) < tol )

    nonexp_amp = np.delete(amp, exp_idx)
    nonexp_mds = np.delete(dmd_mds, exp_idx, axis=1)
    nonexp_ct_eig = np.delete(ct_eig, exp_idx)
    nonexp_dt_eig = np.delete(dt_eig, exp_idx)

    return nonexp_mds, nonexp_ct_eig, nonexp_dt_eig, nonexp_amp


def SortByFrequency( dmd_mds, ct_eig, dt_eig, amp ):
    """    
       Name    : SortByFrequency
       Purpose : Order modes from low to high frequency 

    Parameters
    ----------
       dmd_mds : double,  DMD modes
        ct_eig : double,  continuous-time DMD eigenvalues
        dt_eig : double,  discrete-time DMD eigenvalues
           amp : double,  vector of amplitudes of modes dmd_mds

    """   
    
    sort_index = np.argsort(np.abs(ct_eig))

    sorted_amp = amp[sort_index]
    sorted_mds = osc_mds[:, sort_index]
    sorted_ct_eig = ct_eig[sort_index]
    sorted_dt_eig = dt_eig[sort_index]

    return sorted_mds, sorted_ct_eig, sorted_dt_eig, sorted_amp


class DMD_vecmodes: 
    
    nx = None
    ny = None
    nz = None
    nt = None

    attr = []
    xyz_sizes = []
    xyz_shapes = []
    
    def __init__(self, attr, sizes, shapes, SnapMat, SVDrank, delta_T):
        self.attr = attr
        self.xyz_sizes = sizes
        self.xyz_shapes = shapes
        self.modes, self.dt_eig = DMD_modes( SnapMat[:,:-1], SnapMat[:,1:], SVDrank, delta_T )[:2]
        self.amp = DMD_amplitudes( self.modes, SnapMat[:,0] )
        # Continuous-time eigenvalues
        self.ct_eig = np.log(self.dt_eig)/delta_T
        self.dt = delta_T
        self.is_scaled = False
        
    def setDimension(nx=None, ny=None, nz=None):
        if nx is not None: self.nx = nx
        if ny is not None: self.ny = ny 
        if nz is not None: self.nz = nz 
        
    def CleanUnitCirc(self, tol):
        self.modes, self.ct_eig, self.dt_eig, self.amp = CleanUnitCircle( self.modes, self.ct_eig, self.dt_eig, self.amp, tol)
        
    def CleanExpModes(self, tol):
        self.modes, self.ct_eig, self.dt_eig, self.amp = CleanExponentialModes( self.modes, self.ct_eig, self.dt_eig, self.amp, tol )
        
    def SortByFreq(self):
        self.modes, self.ct_eig, self.dt_eig, self.amp = SortByFrequency( self.modes, self.ct_eig, self.dt_eig, self.amp )
        
    def ScaleWithAmplitude(self):
        self.modes = self.modes * self.amp
        self.is_scaled = True

    def SplitByTreshold(self, treshold_days):
        idm = np.where(np.abs(self.ct_eig.imag) <= np.pi/(treshold_days*86400.))
        mask = np.full(len(self.ct_eig.imag), False, dtype=bool)
        mask[idm] = True

        lowTresh = copy.deepcopy(self)
        highTresh = copy.deepcopy(self)
        
        lowTresh.amp = lowTresh.amp[mask]
        lowTresh.ct_eig = lowTresh.ct_eig[mask]
        lowTresh.dt_eig = lowTresh.dt_eig[mask]
        lowTresh.modes = lowTresh.modes[..., mask]
        lowTresh.nt = len(lowTresh.amp)

        highTresh.amp = highTresh.amp[~mask]
        highTresh.ct_eig = highTresh.ct_eig[~mask]
        highTresh.dt_eig = highTresh.dt_eig[~mask]
        highTresh.modes = highTresh.modes[..., ~mask]
        highTresh.nt = len(highTresh.amp)
        return lowTresh, highTresh

        
    def ToSpatial(self):
        spt_mod = DMD_phyModes( self.amp, self.ct_eig, self.dt_eig ) 
        spt_mod.nm = np.shape(self.modes)[-1]
        spt_mod.nx = self.xyz_shapes[0][-1]
        spt_mod.ny = self.xyz_shapes[0][-2]
        spt_mod.nz = self.xyz_shapes[0][-3]
        spt_mod.attr = self.attr
        Phi = self.modes*self.amp
        for i in range(len(self.attr)):
            data = Phi[i*self.xyz_sizes[1]:(i+1)*self.xyz_sizes[1],:].T
            data = data.reshape((spt_mod.nm, self.xyz_shapes[i][-3],self.xyz_shapes[i][-2],self.xyz_shapes[i][-1]))
            setattr(spt_mod, self.attr[i][:-4]+"mod", data)
        del data
        return spt_mod

    def SaveNetCDF(self, filepath, nm, nx, ny, nz=None):
        spatial = self.ToSpatial( nx, ny, nz)
        spatial.SaveNetCDF(filepath, nm, nx, ny, nz)




class DMD_phyModes:
    
    nx = None
    ny = None
    nz = None
    nt = None

    ddim = 0
    
    
    def __init__(self, amp, ct_eig, dt_eig ):
        self.amp = amp
        self.ct_eig = ct_eig
        self.dt_eig = dt_eig
        self.is_scaled = False

    def ScaleWithAmplitude(self):
        self.Umod = self.Umod * self.amp
        self.Vmod = self.Vmod * self.amp
        if Wmod is not None: self.Wmod = self.Wmod * self.amp
        self.is_scaled = True

    def SplitByTreshold(self, treshold_days):
        idm = np.where(np.abs(self.ct_lam.imag) <= np.pi/(treshold_days*86400.))
        mask = np.full(len(self.ct_lam.imag), False, dtype=bool)
        mask[idm] = True
        
        lowTresh = self
        lowTresh.amp = self.amp[mask]
        lowTresh.ct_eig = self.ct_eig[mask]
        lowTresh.dt_eig = self.dt_eig[mask]
        lowTresh.nt = len(lowTresh.amp)
          
        for i in range(len(highTresh.attr)):
            data = getattr(self, self.attr[i][:-4]+"mod") 
            data = data[mask,...]
            setattr(highTresh, self.attr[i][:-4]+"mod", data)
        del data
        
        
        highTresh = self   
        highTresh.amp = self.amp[~mask]
        highTresh.ct_eig = self.ct_eig[~mask]
        highTresh.dt_eig = self.dt_eig[~mask]
        highTresh.nt = len(highTresh.amp)
          
        for i in range(len(highTresh.attr)):
            data = getattr(self, self.attr[i][:-4]+"mod") 
            data = data[~mask,...]
            setattr(highTresh, self.attr[i][:-4]+"mod", data)
        del data
    
    
    def SaveNetCDF(self, filepath, nm, nx, ny, nz=None):
        """
        
           Name    : SaveNetCDF
           Purpose : Saves Real and Imaginary part in a (split) NetCDF dataset
    
        Parameters
        ----------
              self : class,       DMD class
          filepath : char,        Name of the output file
                nm : int,         Number of modes
                nx : int,         Number of x gridpoints
                ny : int,         Number of y gridpoints
                nz : int,         Number of vertical levels
    
           Returns
           -------
              None : Writes external file 
    
        """
        from netCDF4 import Dataset
        
        Phi = dmd_mds*dmd_amp
        
        fout = Dataset(filepath , 'w', format='NETCDF4')
        
        # Create dimensions
        fout.createDimension('mode', nm)
        
        fout.createDimension('y', ny)
        fout.createDimension('x', nx)
        nz_ = 1 if (nz is None) else nz 
        if nz is not None: fout.createDimension('z', nz)
        
        # Create variables
        mid = fout.createVariable('mode', 'i4', ('mode',))
        
        xid = fout.createVariable('x', 'f8', ('x',))
        yid = fout.createVariable('y', 'f8', ('y',))
        if nz is not None: zid = fout.createVariable('z', 'f8', ('z',))
        
        lreid = fout.createVariable('lambda_real', 'f8', ('mode',))
        limid = fout.createVariable('lambda_imag', 'f8', ('mode',))
        
        ctlreid = fout.createVariable('ct_lamb_real', 'f8', ('mode',))
        ctlimid = fout.createVariable('ct_lamb_imag', 'f8', ('mode',))
        
        breid = fout.createVariable('b_real', 'f8', ('mode',))
        bimid = fout.createVariable('b_imag', 'f8', ('mode',))
        
        if nz is not None:
            umid = fout.createVariable('umean', 'f8', ('z','y','x',))
            vmid = fout.createVariable('vmean', 'f8', ('z','y','x',))
            wmid = fout.createVariable('wmean', 'f8', ('z','y','x',))
            
            ureid = fout.createVariable('umode_real', 'f8', ('mode','z','y','x',))
            uimid = fout.createVariable('umode_imag', 'f8', ('mode','z','y','x',))
            
            vreid = fout.createVariable('vmode_real', 'f8', ('mode','z','y','x',))
            vimid = fout.createVariable('vmode_imag', 'f8', ('mode','z','y','x',))
            
            wreid = fout.createVariable('wmode_real', 'f8', ('mode','z','y','x',))
            wimid = fout.createVariable('wmode_imag', 'f8', ('mode','z','y','x',))
        else:
            umid = fout.createVariable('umean', 'f8', ('y','x',))
            vmid = fout.createVariable('vmean', 'f8', ('y','x',))
            
            ureid = fout.createVariable('umode_real', 'f8', ('mode','y','x',))
            uimid = fout.createVariable('umode_imag', 'f8', ('mode','y','x',))
            
            vreid = fout.createVariable('vmode_real', 'f8', ('mode','y','x',))
            vimid = fout.createVariable('vmode_imag', 'f8', ('mode','y','x',))
            
        # Add attributes
        mid.long_name = 'Mode index'
        
        xid.long_name = 'Ocean X axis (T-grid)'
        yid.long_name = 'Ocean Y axis (T-grid)'
        zid.long_name = 'Ocean mid-layer axis'
        xid.units = 'km'
        yid.units = 'km'
        zid.units = 'km'
        
        umid.long_name = 'Mean of zonal velocity'
        vmid.long_name = 'Mean of meridional velocity'
        umid.units = 'm/s'
        vmid.units = 'm/s'
        
        if nz is not None:    
            wmid.long_name = 'Mean of vertical velocity'
            wmid.units = 'm/s'
        
        lreid.long_name = 'DMD eigenvalues (real part)'
        limid.long_name = 'DMD eigenvalues (imag part)'
        ctlreid.long_name = 'Continuous time eigenvalues (real part)'
        ctlimid.long_name = 'Continuous time eigenvalues (imag part)'
        breid.long_name = 'DMD amplitudes (real part)'
        bimid.long_name = 'DMD amplitudes (imag part)'
        
        ureid.long_name = 'Zonal DMD modes (real part)'
        uimid.long_name = 'Zonal DMD modes (imag part)'
        vreid.long_name = 'Meridional DMD modes (real part)'
        vimid.long_name = 'Meridional DMD modes (imag part)'
        
        if nz is not None:    
            wreid.long_name = 'Vertical DMD modes (real part)'
            wimid.long_name = 'Vertical DMD modes (imag part)'
        
    
        # Write data
        mid = range(nm)
        yid = range(ny)
        xid = range(nx)
        if nz is not None: zid[:] = range(nz)
            
        lreid = self.dt_eig.real
        limid = self.dt_eig.imag
        ctlreid = self.ct_eig.real
        ctlimid = self.ct_eig.imag
        breid = self.dmd_amp.real
        bimid = self.dmd_amp.imag
        #umid = Uavg
        #vmid = Vavg
        #if nz is not None: wmid = Wavg
           
        ureid = self.Umod.real
        uimid = self.Umod.imag
        vreid = self.Vmod.real
        vimid = self.Vmod.imag
        
        if nz is not None: 
            wreid = self.Wmod.real
            wimid = self.Wmod.imag
        
        fout.close()
        print("Saving on NetCDF file:")
        print("   " + filepath)



        

def SaveNetCDF(dmd_corr, dmd_uncorr, utrend, vtrend, filepath ):
    # Save outputs
    fout = Dataset(filepath, 'w', format='NETCDF4')
    
    fout.createDimension('mode_c', dmd_corr.nm)
    fout.createDimension('mode_r', dmd_uncorr.nm)
    
    fout.createDimension('y', dmd_corr.ny)
    fout.createDimension('x', dmd_corr.nx)
    fout.createDimension('z', dmd_corr.nz)
    
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
    
    # xid.long_name = 'Ocean X axis'
    # yid.long_name = 'Ocean Y axis'
    # zid.long_name = 'Ocean Z axis'
    # xid.units = 'km'
    # yid.units = 'km'
    # zid.units = 'km'
    
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
    mcid[:] = range(dmd_corr.nm)
    mrid[:] = range(dmd_uncorr.nm)
    
    # yid[:] = dmd_uncorr.y
    # xid[:] = dmd_uncorr.x
    # zid[:] = dmd_uncorr.z
    
    wcid[:] = dmd_corr.ct_eig.imag
    wrid[:] = dmd_uncorr.ct_eig.imag
    
    umid[...] = utrend 
    vmid[...] = vtrend
    #wmid[...] = wc
    
    ucreid[...] = dmd_corr.Umod.real
    vcreid[...] = dmd_corr.Vmod.real
    wcreid[...] = dmd_corr.Wmod.real
    
    ucimid[...] = dmd_corr.Umod.imag 
    vcimid[...] = dmd_corr.Vmod.imag
    wcimid[...] = dmd_corr.Wmod.imag
    
    urreid[...] = dmd_uncorr.Umod.real
    vrreid[...] = dmd_uncorr.Vmod.real
    wrreid[...] = dmd_uncorr.Wmod.real
    
    urimid[...] = dmd_uncorr.Umod.imag
    vrimid[...] = dmd_uncorr.Vmod.imag
    wrimid[...] = dmd_uncorr.Wmod.imag
    
    # Close file
    fout.close()
    
    print("")
    print('Program terminates')

