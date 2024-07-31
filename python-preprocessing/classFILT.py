#!/usr/bin/env python3

import numpy as np
import scipy.ndimage
from scipy import signal
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


def SpatialFiltering(hres, lres, scal, field, fnam=None):
    """   

       Name    : SpatialFiltering
       Purpose : Applies given spatial filter to a given field

    Parameters
    ----------
          hres : int,     Original resolution (that of field)
          lres : int,     Target resolution (that of flt)
          scal : double,  Scaling parameter for sigma
         field : double,  Input field of resolution hres
        

       Returns
       -------
           flp : double,  Filtered field
           fhp : double,  High pass filter ( as I-F )

          Note : field is not modified, fct has the same dimension
          ----

    """   
    
    # %%
    if fnam == None: 
        fnam = ""
    else:
        fnam += " "

    ratio= int(hres/lres)
    flt_sigma = (scal*ratio, scal*ratio)
    flp = np.zeros_like(field)

    print("Spatial filtering " + fnam)
    flp = scipy.ndimage.gaussian_filter(field, sigma=flt_sigma, mode='constant', axes=(-2,-1))
    
    # Set bottom boundary contition 
    if field.ndim == 4:
        flp[..., -1, :, :] = 0.
   
    fhp = field - flp 

    return flp, fhp


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
    flt_sigma = (scal*ratio, scal*ratio)

    print("Downsampling filtered " + fnam)
    fds = scipy.ndimage.gaussian_filter(field, sigma=flt_sigma, mode='constant', axes=(-2,-1))[..., ::ratio, ::ratio]
    
    return fds

class sp_filt:
    
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
            #
            # Filter High resolution in space
            vel_hlp, vel_hhp = SpatialFiltering(hres, lres, scal, np.squeeze(dict[key]), key + "vel")
            #
            # Filter another time and downsample
            vel_flt = SpatialFiltDown(hres, lres, scal, vel_hlp, "low passed high resolution " + key + "vel")
            vel_fhp = SpatialFiltDown(hres, lres, 2*scal, vel_hhp, "high passed high resolution " + key + "vel")

            vel_avg = np.mean(vel_flt, axis=0)
            vel_ds = SpatialDownsample(hres, lres, np.squeeze(dict[key]), key + "vel")

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
            tempdict["data"] = vel_avg
            tempdict["name"] = key + "avg_ds" 
            tempdict["shape"] = tempdict["data"].shape
            tempdict["units"] = "m/s"
            tempdict["longname"] = "Temporal average of filtered " + key.upper() + " velocity " + key + "flt"
            setattr(self, key+"avg_ds", tempdict)
            self.vars.append(key+"avg_ds")

            # Define fluctuations field
            tempdict = {}
            tempdict["data"] = vel_fhp - np.mean(vel_fhp, axis=0)

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
            tempdict["data"] = vel_fhp
            tempdict["name"] = key + "_hp" 
            tempdict["shape"] = tempdict["data"].shape
            tempdict["units"] = "m/s"
            tempdict["longname"] = "Filtered high passed " + key.upper() + " velocity ( " + key + "vel - F(" + key + "vel) )"
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

        dims = ("scalar",)
        hres = fout.createVariable("hres", "i4", dims)
        hres.units = "1/degrees"
        hres.long_name = "High resolution factor"
        hres = self.hres
        
        lres = fout.createVariable("lres", "i4", dims)
        lres.units = "1/degrees"
        lres.long_name = "Low resolution factor"
        lres = self.lres

        scal = fout.createVariable("scal", "i4", dims)
        scal.long_name = "Resolution ratio scale: sigma=scal*Hres/Lres"
        scal = self.scal

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
            var_id.long_name = var_dict["longname"]

            tmp[...,1:-1,1:-1] = var_dict["data"] 
            if var_dict["data"].ndim == 3: var_id[:,:,:] = tmp  
            if var_dict["data"].ndim == 4: var_id[:,:,:,:] = tmp  

        fout.close()




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

class tm_filt:
    
    def __init__(self, nt, levs, hres, lres, scal, fcut, fs, listofdicts, HighPass=False):
        self.hres = hres
        self.lres = lres
        self.scal = scal
        self.fcut = fcut
        self.fs = fs
        self.vars = []

        self.nt = nt
        self.levs = levs
        self.HighPass = HighPass

        for dict in listofdicts:
            key = next(iter(dict))
            print("Processing " + key)
            
            # Temporal filtering
            vel_hlp, vel_hhp = TemporalFiltering(fcut, fs, np.squeeze(dict[key]))

            vel_lp = SpatialFiltDown(hres, lres, scal, vel_hlp, key + "vel")
            vel_hp = SpatialFiltDown(hres, lres, scal, vel_hhp, key + "vel")

            vel_ds = SpatialDownsample(hres, lres, np.squeeze(dict[key]), key + "vel")
            avg_ds = np.mean(vel_lp, axis=0)
            vel_fct = vel_lp - avg_ds


            # Define filtered field
            tempdict = {}
            tempdict["data"] = vel_lp
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

        dims = ("scalar",)
        hres = fout.createVariable("hres", "i4", dims)
        hres.units = "1/degrees"
        hres.long_name = "High resolution factor"
        hres = self.hres
        
        lres = fout.createVariable("lres", "i4", dims)
        lres.units = "1/degrees"
        lres.long_name = "Low resolution factor"
        lres = self.lres

        scal = fout.createVariable("scal", "i4", dims)
        scal.long_name = "Resolution ratio scale: sigma=scal*Hres/Lres"
        scal = self.scal

        fcut = fout.createVariable("fcut", "i4", dims)
        fcut.units = "1/s"
        fcut.long_name = "Cut frequency"
        fcut = self.fcut

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
            var_id.long_name = var_dict["longname"]

            tmp[...,1:-1,1:-1] = var_dict["data"] 
            if var_dict["data"].ndim == 3: var_id[:,:,:] = tmp  
            if var_dict["data"].ndim == 4: var_id[:,:,:,:] = tmp  

        fout.close()


