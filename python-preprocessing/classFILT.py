#!/usr/bin/env python3

import numpy as np
import scipy.ndimage
from netCDF4 import Dataset
    


def SpatialDownsample(hres, lres, field, fnam=None):
    """   

       Name    : SpatialDownsample
       Purpose : Subsamples a given field

    Parameters
    ----------
          hres : int,     Original resolution (that of field)
          lres : int,     Target resolution (that of flt)
         field : double,  Input field of resolution hres
        

       Returns
       -------
           fds : double,  Downsampled field

          Note : field is not modified
          ----

    """   
    
    # %%
    if fnam == None: 
        fnam = ""
    else:
        fnam += " "
    
    ratio= int(hres/lres)

    print("Downsampling " + fnam)
    fds = field[..., ::ratio, ::ratio]
    
    return fds


def SpatialFiltDown(hres, lres, scal, field, fnam=None):
    """   

       Name    : SpatialFiltDown
       Purpose : Applies given spatial filter to a given field and subsamples it

    Parameters
    ----------
          hres : int,     Original resolution (that of field)
          lres : int,     Target resolution (that of flt)
          scal : double,  Scaling parameter for sigma
         field : double,  Input field of resolution hres
        

       Returns
       -------
           flt : double,  Filtered field
           avg : double,  0th-Axis average
           fct : double,  Fluctuation field

          Note : field is not modified, it assumes Axis=0 is time, always
          ----

    """   
    
    # %%
    if fnam == None: 
        fnam = ""
    else:
        fnam += " "

    ratio= int(hres/lres)
    flt_sigma = scal*ratio

    print("Filtering and downsampling " + fnam)
    flt = scipy.ndimage.gaussian_filter(field, sigma=flt_sigma, mode='constant')[..., ::ratio, ::ratio]
    
    # Set bottom boundary contition 
    if field.ndim == 4:
        flt[..., -1, :, :] = 0.
   
    if field.ndim > 2:
        print("0th-Axis averaging of " + fnam)
        avg = np.mean(flt, axis=0)
        fct = flt - avg
    else:
        print("Array " + fnam + "is two-dimensional in space, avg and fct will be zeros")
        avg = np.zeros_like(flt)
        fct = np.zeros_like(flt)

    return flt, avg, fct




def TemporalFiltering(fcut, fs, field):
    """    
       
       Name    : TemporalFiltering
       Purpose : Applies given temporal filter to a given field 

    Parameters
    ----------
          fcut : double,  Cutting frequency for the filter
            fs : double,  Sampling frequency
         field : double,  discrete-time DMD 
        

       Returns
       -------
           flt : double,  Filtered field
           fct : double,  Fluctuation field

          Note : field is not modified
          ----

    """   
    
    print('Perform eddy-mean decomposition')
    wid = fcut / (0.5/fs)
    b, a = signal.butter(5, wid, 'low') # 5th order Butterworth
    flt = signal.filtfilt(b, a, field, axis=0) # mean
    fct = field - flt # eddy  

    return flt, fct




class sp_filt:
    
    nx = None
    ny = None
    nz = None
    nt = None

    ddim = 0
    
    
    def __init__(self, nt, levs, hres, lres, scal, listofdicts ):
        self.hres = hres
        self.lres = lres
        self.scal = scal
        self.vars = []

        self.nt = nt
        self.levs = levs

        for dict in listofdicts:
            key = next(iter(dict))
            print("Processing " + key)

            vel_flt, avg_ds, vel_fct = SpatialFiltDown(hres, lres, scal, np.squeeze(dict[key]), key + "vel")
            vel_ds = SpatialDownsample(hres, lres, np.squeeze(dict[key]), key + "vel")
            vel_hp = vel_ds - vel_flt

            # Define filtered field
            tempdict = {}
            tempdict["data"] = vel_flt
            tempdict["name"] = key + "flt" 
            tempdict["shape"] = tempdict["data"].shape
            tempdict["units"] = "m/s"
            tempdict["longname"] = "Filtered and downsampled " + key.upper() + " velocity"
            setattr(self, key+"flt", tempdict)
            self.vars.append(key+"flt")

            # Define averaged filtered field
            tempdict = {}
            tempdict["data"] = avg_ds
            tempdict["name"] = key + "avg_ds" 
            tempdict["shape"] = tempdict["data"].shape
            tempdict["units"] = "m/s"
            tempdict["longname"] = "Temporal average of filtered " + key.upper() + " velocity " + key + "flt"
            setattr(self, key+"avg_ds", tempdict)
            self.vars.append(key+"avg_ds")

            # Define fluctuations field
            tempdict = {}
            tempdict["data"] = vel_fct
            tempdict["name"] = key + "fct" 
            tempdict["shape"] = tempdict["data"].shape
            tempdict["units"] = "m/s"
            tempdict["longname"] = "Temporal fluctuations of filtered " + key.upper() + " velocity " + key + "flt"
            setattr(self, key+"fct", tempdict)
            self.vars.append(key+"fct")

            # Define downsampled field
            tempdict = {}
            tempdict["data"] = vel_ds
            tempdict["name"] = key + "_ds" 
            tempdict["shape"] = tempdict["data"].shape
            tempdict["units"] = "m/s"
            tempdict["longname"] = "Downsampled " + key.upper() + " velocity (no filter)"
            setattr(self, key+"_ds", tempdict)
            self.vars.append(key+"_ds")

            # Define high-passed field
            tempdict = {}
            tempdict["data"] = vel_hp
            tempdict["name"] = key + "_hp" 
            tempdict["shape"] = tempdict["data"].shape
            tempdict["units"] = "m/s"
            tempdict["longname"] = "High passed " + key.upper() + " velocity ( " + key + "_ds - " + key + "_flt )"
            setattr(self, key+"_hp", tempdict)
            self.vars.append(key+"_hp")



    def SpatialNetCDFSaveOut(self, outfile):

        print("Writing filtered field in file: ", outfile)

        fout = Dataset(outfile, "w", format="NETCDF4")
   

        nx = 30*self.lres + 2
        ny = 20*self.lres + 2
        nz = self.levs
        nt = self.nt

        # Create dimensions
        fout.createDimension("scalar", 0)
        fout.createDimension("t", nt)
        fout.createDimension("x", nx)
        fout.createDimension("y", ny)
        fout.createDimension("z", nz)

        for var in self.vars:
            var_dict = getattr(self, var)

            if var_dict["data"].ndim == 3: 
                # Check if last variable is nz
                if var_dict["data"].shape[0] == nz: 
                    dims = ("z", "y", "x",)
                    tmp = np.zeros([nz, ny, nx])
                # Check if last variable is nt
                if var_dict["data"].shape[0] == nt: 
                    dims = ("t", "y", "x",)
                    tmp = np.zeros([nt, ny, nx])


            if var_dict["data"].ndim == 4: 
                dims = ("t", "z", "y", "x",)
                tmp = np.zeros([nt, nz, ny, nx])

            var_id = fout.createVariable(var, "f4", dims)
            var_id.units = var_dict["units"]
            var_id.long_name = ["longname"]

            tmp[...,1:-1,1:-1] = var_dict["data"] 
            if var_dict["data"].ndim == 3: var_id[:,:,:] = tmp  
            if var_dict["data"].ndim == 4: var_id[:,:,:,:] = tmp  

        fout.close()
