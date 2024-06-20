#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 23 16:55:00 2023


@author: ftucciar
"""

import numpy as np


class Phys_data:

    nx = None
    ny = None
    nz = None
    nt = None

    attr = []
    sizes = []
    shapes = []
   

    dict = {}

    def __init__(self, nx=None, ny=None, nz=None, nt=None):
        """
    
           Name    : __init__
           Purpose : Constructor

        Parameters
        ----------
            fields : double,      Input to detrend.

           Returns
           -------
              trnd : double,      Trend
           detrndd : double,      Detrended input
    
        """
        self.attr = []
        self.sizes = []
        self.shapes = []
        if nx is not None: self.nx = nx
        if ny is not None: self.ny = ny
        if nz is not None: self.nz = nz
        if nt is not None: self.nt = nt
       

    def append(self, varname, data, move=None):
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
        setattr(self, varname, data)
        self.shapes.append(np.array(np.shape(data)[1:]))
        self.sizes.append(np.size(data[0,...]))
        self.attr.append(varname)
       

    def time_detrend(self, varnam_in, varnam_out, move=None):
        """
    
           Name    : time_detrend
           Purpose : Removes the (moving or static) mean of a dataset

        Parameters
        ----------
            fields : double,      Input to detrend.

           Returns
           -------
              trnd : double,      Trend
           detrndd : double,      Detrended input

        """

        if move is None:
            trnd = np.mean(getattr(self, varnam_in ), axis=0, dtype=np.float64)
            dtrnd = getattr(self, varnam_in ) - trnd
            
            setattr(self, varnam_out + "trnd", trnd )
            setattr(self, varnam_out + "fctn", dtrnd )
        else: 
            if (move % 2) == 0:
                print("Adding 1 to move")
                move +=1 
            
            field = getattr(self, varnam_in )
            cumsum = np.cumsum(field, axis=0)
            trnd = ( cumsum[:-move, ...] - cumsum[move:,...] ) / move 
            dtrnd = field - trnd
            
            setattr(self, varnam_out + "trnd", trnd )
            setattr(self, varnam_out + "fctn", dtrnd )


    def unroll(self, varnam):
        """
    
           Name    : time_detrend
           Purpose : Removes the (moving or static) mean of a dataset

        Parameters
        ----------
            fields : double,      Input to detrend.

           Returns
           -------
              trnd : double,      Trend
           detrndd : double,      Detrended input

        """
        
        nx_ = 1 if (self.nx is None) else self.nx 
        ny_ = 1 if (self.ny is None) else self.ny 
        nz_ = 1 if (self.nz is None) else self.nz
        
        return getattr(self, varnam ).reshape((self.nt,nz_*ny_*nx_))
    
    
    def rollup(self, field_in):
        """
    
           Name    : time_detrend
           Purpose : Removes the (moving or static) mean of a dataset

        Parameters
        ----------
            fields : double,      Input to detrend.

           Returns
           -------
              trnd : double,      Trend
           detrndd : double,      Detrended input

        """
        
        nx_ = 1 if (self.nx is None) else self.nx 
        ny_ = 1 if (self.ny is None) else self.ny 
        nz_ = 1 if (self.nz is None) else self.nz
        
        return field_in.reshape((self.nt,nz_, ny_, nx_))
        
