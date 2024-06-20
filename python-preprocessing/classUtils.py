#!/usr/bin/env python3
#! -*- coding: utf-8 -*-


# Load modules
# ------------
import numpy as np
import copy
from netCDF4 import Dataset




def mean_removal(*fields, move=None):
    """
    
       Name    : mean_removal
       Purpose : Removes the (moving or static) mean of a dataset

    Parameters
    ----------
        fields : double,      Input to detrend.

       Returns
       -------
          trnd : double,      Trend
       detrndd : double,      Detrended input

    """
    trnd = []
    dtrnd = []
    if move is None:
        for idx, field in enumerate(fields):
            trnd.append( np.mean(field, axis=0, dtype=np.float64) )
            dtrnd.append( field - np.mean(field, axis=0, dtype=np.float64) )
    else: 
        if (move % 2) == 0:
            print("Adding 1 to move")
            move +=1 
        for idx, field in enumerate(fields):
            cumsum = np.cumsum(field, axis=0)
            trnd.append( (cumsum[:-move, ...] - cumsum[move:,...])/move)
            dtrnd.append( field[:-move, ...] - (cumsum[:-move, ...] - cumsum[move:,...])/move)
    return trnd, dtrnd

        
def ReadComplexField(field_name):
    """
    
       Name    : ReadComplexField
       Purpose : Reads Real and Imaginary part from a (split) NetCDF dataset

    Parameters
    ----------
    field_name : char,           Name of the input field to read from NetCDF

       Returns
       -------
    field_cplx : double complex, Associated complex field

    """
    re = f.variables[field_name + "_real"][:].copy() 
    im = f.variables[field_name + "_imag"][:].copy()
    field_cplx = re + 1j*im
    return field_cplx



    
    
    
def SaveDMDsplitNetCDF(dmd_mds, dmd_amp, ct_eig, dt_eig, filepath, nm_r, nm_c, nx, ny, nz=None):
    """
    
       Name    : ReadComplexField
       Purpose : Reads Real and Imaginary part from a (split) NetCDF dataset

    Parameters
    ----------
       dmd_mds : double,      DMD modes
       dmd_amp : double,      Vector of amplitudes of modes dmd_mds
        ct_eig : double,      Continuous-time DMD eigenvalues
        dt_eig : double,      Discrete-time DMD eigenvalues
      filepath : char,        Name of the output file
          nm_r : int,         Number of random modes
          nm_c : int,         Number of correlated modes
            nx : int,         Number of x gridpoints
            ny : int,         Number of y gridpoints
            nz : int,         Number of vertical levels

       Returns
       -------
          None : Writes external file 
    

    """   
    # Save outputs
    fout = Dataset(filepath, 'w', format='NETCDF4')
    
    fout.createDimension('mode_r', nm_r)
    fout.createDimension('mode_c', nm_c)
    
    fout.createDimension('y', ny)
    fout.createDimension('x', nx)
    if nz is not None: fout.createDimension('z', nz)
    
    yid = fout.createVariable('y', 'f8', ('y',))
    xid = fout.createVariable('x', 'f8', ('x',))
    if nz is not None: zid = fout.createVariable('z', 'f8', ('z',))
    
    mcid = fout.createVariable('mode_c', 'i4', ('mode_c',))
    mrid = fout.createVariable('mode_r', 'i4', ('mode_r',))
    
    wcid = fout.createVariable('omega_c', 'f8', ('mode_c',))
    wrid = fout.createVariable('omega_r', 'f8', ('mode_r',))
        
    if nz is not None: 
        umid = fout.createVariable('uco', 'f8', ('z','y','x',))
        vmid = fout.createVariable('vco', 'f8', ('z','y','x',))
        wmid = fout.createVariable('wco', 'f8', ('z','y','x',))
        
        ucreid = fout.createVariable('umode_real_c', 'f8', ('mode_c','z','y','x',))
        ucimid = fout.createVariable('umode_imag_c', 'f8', ('mode_c','z','y','x',))
        urreid = fout.createVariable('umode_real_r', 'f8', ('mode_r','z','y','x',))
        urimid = fout.createVariable('umode_imag_r', 'f8', ('mode_r','z','y','x',))
        vcreid = fout.createVariable('vmode_real_c', 'f8', ('mode_c','z','y','x',))
        vcimid = fout.createVariable('vmode_imag_c', 'f8', ('mode_c','z','y','x',))
        vrreid = fout.createVariable('vmode_real_r', 'f8', ('mode_r','z','y','x',))
        vrimid = fout.createVariable('vmode_imag_r', 'f8', ('mode_r','z','y','x',))
        
        wcreid = fout.createVariable('wmode_real_c', 'f8', ('mode_c','z','y','x',))
        wcimid = fout.createVariable('wmode_imag_c', 'f8', ('mode_c','z','y','x',))
        wrreid = fout.createVariable('wmode_real_r', 'f8', ('mode_r','z','y','x',))
        wrimid = fout.createVariable('wmode_imag_r', 'f8', ('mode_r','z','y','x',))
    else:
        umid = fout.createVariable('uco', 'f8', ('y','x',))
        vmid = fout.createVariable('vco', 'f8', ('y','x',))
        wmid = fout.createVariable('wco', 'f8', ('y','x',))
        
        ucreid = fout.createVariable('umode_real_c', 'f8', ('mode_c','y','x',))
        ucimid = fout.createVariable('umode_imag_c', 'f8', ('mode_c','y','x',))
        urreid = fout.createVariable('umode_real_r', 'f8', ('mode_r','y','x',))
        urimid = fout.createVariable('umode_imag_r', 'f8', ('mode_r','y','x',))
        vcreid = fout.createVariable('vmode_real_c', 'f8', ('mode_c','y','x',))
        vcimid = fout.createVariable('vmode_imag_c', 'f8', ('mode_c','y','x',))
        vrreid = fout.createVariable('vmode_real_r', 'f8', ('mode_r','y','x',))
        vrimid = fout.createVariable('vmode_imag_r', 'f8', ('mode_r','y','x',))
        
        wcreid = fout.createVariable('wmode_real_c', 'f8', ('mode_c','y','x',))
        wcimid = fout.createVariable('wmode_imag_c', 'f8', ('mode_c','y','x',))
        wrreid = fout.createVariable('wmode_real_r', 'f8', ('mode_r','y','x',))
        wrimid = fout.createVariable('wmode_imag_r', 'f8', ('mode_r','y','x',))
    
    # Add attributes
    mcid.long_name = 'Mode index (for correction drift)'
    mrid.long_name = 'Mode index (for random noise)'
    
    xid.long_name = 'Ocean X axis'
    yid.long_name = 'Ocean Y axis'
    if nz is not None: zid.long_name = 'Ocean Z axis'
    
    xid.units = 'km'
    yid.units = 'km'
    if nz is not None: zid.units = 'km'
    
    wcid.long_name = 'Continuous eigenvalues for correction drift'
    wrid.long_name = 'Continuous eigenvalues for random noise'
    wcid.units = 's^-1'
    wrid.units = 's^-1'
    
    umid.long_name = 'Mean of zonal correction drift'
    ucreid.long_name = 'Zonal DMD modes for correction drift (real part)'
    ucimid.long_name = 'Zonal DMD modes for correction drift (imag part)'
    urreid.long_name = 'Zonal DMD modes for random noise (real part)'
    urimid.long_name = 'Zonal DMD modes for random noise (imag part)'
    umid.units = 'm/s'
    ucreid.units = 'm/s'
    ucimid.units = 'm/s'
    urreid.units = 'm/s'
    urimid.units = 'm/s'
    
    vmid.long_name = 'Mean of meridional correction drift'
    vcreid.long_name = 'Meridional DMD modes for correction drift (real part)'
    vcimid.long_name = 'Meridional DMD modes for correction drift (imag part)'
    vrreid.long_name = 'Meridional DMD modes for random noise (real part)'
    vrimid.long_name = 'Meridional DMD modes for random noise (imag part)'
    vmid.units = 'm/s'
    vcreid.units = 'm/s'   
    vcimid.units = 'm/s'
    vrreid.units = 'm/s'
    vrimid.units = 'm/s'
    
    if nz is not None:     
        wcreid.long_name = 'Vertical DMD modes for correction drift (real part)'
        wcimid.long_name = 'Vertical DMD modes for correction drift (imag part)'
        wrreid.long_name = 'Vertical DMD modes for random noise (real part)'
        wrimid.long_name = 'Vertical DMD modes for random noise (imag part)'
        wcimid.units = 'm/s'
        wcreid.units = 'm/s'
        wrreid.units = 'm/s'
        wrimid.units = 'm/s'
    
    # Write data
    mcid = range(nm_c)
    mrid = range(nm_r)
    
    yid = y
    xid = x
    if nz is not None: zid = z
    
    wcid = omg_c.imag
    wrid = omg_r.imag
    
    umid = uc
    ucreid = um_c.real
    ucimid = um_c.imag
    urreid = um_r.real
    urimid = um_r.imag
    
    vmid = vc
    vcreid = vm_c.real
    vcimid = vm_c.imag
    vrreid = vm_r.real
    vrimid = vm_r.imag
    
    if nz is not None:
        wmid = wc
        wcreid = wm_c.real
        wcimid = wm_c.imag
        wrreid = wm_r.real
        wrimid = wm_r.imag
    
    # Close file
    fout.close()
    print("Saving on NetCDF file:")
    print("   " + filename)
    
    
    
    
    
    
    
    
    
    
    
    
