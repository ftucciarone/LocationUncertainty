#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 11:02:13 2021

         Computes the EOF through an Proper Orthogonal Decomposition
         procedure starting of a 3D velocity field.
         Prints the output in 3 different files

@author: ftucciar
"""

# Load modules
# ------------
import os
import sys
import yaml
import pathlib


def yaml2dict(path):
    print(path)
    # Open YAML file to read paths
    # ----------------------------
    # work_dir: folder containing the LU procedures
    # data_dir: containing the i/o data
    #  dom_dir: containing the domain files from NEMO
    #  out_dir: output parent directory
    with open(path) as file:
        # The FullLoader parameter handles the conversion from YAML
        # scalar values to Python the dictionary format
        IO_dict = yaml.load(file, Loader=yaml.FullLoader)
        HR_dict = IO_dict['HRexps']
        LR_dict = IO_dict['LRexps']
    
    # Set directories
    # ---------------
    # home_dir: depending on the system
    home_dir = str(pathlib.Path.home())
    IO_dict['home_dir'] = str(pathlib.Path.home())
    IO_dict['work_dir'] = home_dir + IO_dict['work_dir']
    IO_dict['data_dir'] = home_dir + IO_dict['data_dir']
    IO_dict['dom_dir'] = home_dir + IO_dict['dom_dir']
    IO_dict['out_dir'] = IO_dict['work_dir'] + IO_dict['out_dir']
    
    print("")
    print("Sanity check on existence of folders and files")
    print("")
    print("  isdir =", os.path.isdir(IO_dict['home_dir']), "    home_dir = ", IO_dict['home_dir'])
    print("  isdir =", os.path.isdir(IO_dict['work_dir']), "    work_dir = ", IO_dict['work_dir'])
    print("  isdir =", os.path.isdir(IO_dict['data_dir']), "    data_dir = ", IO_dict['data_dir'])
    print("  isdir =", os.path.isdir(IO_dict['dom_dir']), "     dom_dir = ", IO_dict['dom_dir'])
    print("  isdir =", os.path.isdir(IO_dict['out_dir']), "     out_dir = ", IO_dict['out_dir'])
    
    print("")
    print("Sanity check on existence of High Resolution folders and files")
    print("")
    HRfiles = []
    for HRexp in range(len(HR_dict)):
        tempdict = {}
        temp =  HR_dict[list(HR_dict)[HRexp]]
        tempdir = IO_dict['data_dir'] + temp['base_dir']
        templist = []
        print("  isdir =", os.path.isdir(tempdir), "data_basedir = ", tempdir)
        print("") 
        for i in range(len(temp['subs_dir'])):
            tempfile = tempdir + list(temp['subs_dir'])[i] + "/" + temp['curlfile'] + ".nc"
            templist.append(tempfile)
            print(" isfile =", os.path.isfile(templist[i]) , "                " + templist[i] ) 
        tempdict["curl_files"] = templist
        templist = []
        print("")
        for i in range(len(temp['subs_dir'])):
            tempfile = tempdir + list(temp['subs_dir'])[i] + "/" + temp['tempfile'] + ".nc"
            templist.append(tempfile)
            print(" isfile =", os.path.isfile(templist[i]) , "                " + templist[i] ) 
        tempdict["temp_files"] = templist
        templist = []
        print("")
        for i in range(len(temp['subs_dir'])):
            tempfile = tempdir + list(temp['subs_dir'])[i] + "/" + temp['uvelfile'] + ".nc"
            templist.append(tempfile)
            print(" isfile =", os.path.isfile(templist[i]) , "                " + templist[i] ) 
        tempdict["uvel_files"] = templist
        templist = []
        print("")
        for i in range(len(temp['subs_dir'])):
            tempfile = tempdir + list(temp['subs_dir'])[i] + "/" + temp['vvelfile'] + ".nc"
            templist.append(tempfile)
            print(" isfile =", os.path.isfile(templist[i]) , "                " + templist[i] ) 
        tempdict["vvel_files"] = templist
        templist = []
        print("")
        for i in range(len(temp['subs_dir'])):
            tempfile = tempdir + list(temp['subs_dir'])[i] + "/" + temp['wvelfile'] + ".nc"
            templist.append(tempfile)
            print(" isfile =", os.path.isfile(templist[i]) , "                " + templist[i] ) 
        tempdict["wvel_files"] = templist
        HRfiles.append(tempdict)
        templist = []
        del temp
        print("")
    
    print("Sanity check on existence of Low Resolution folders and files")
    print("")
    LRfiles = []
    for LRexp in range(len(LR_dict)):
        tempdict = {}
        temp =  LR_dict[list(LR_dict)[LRexp]]
        tempdir = IO_dict['data_dir'] + temp['base_dir']
        templist = []
        print("  isdir =", os.path.isdir(tempdir), "data_basedir = ", tempdir)
        print("")
        for i in range(len(temp['subs_dir'])):
            tempfile = tempdir + list(temp['subs_dir'])[i] + "/" + temp['curlfile'] + ".nc"
            templist.append(tempfile)
            print(" isfile =", os.path.isfile(templist[i]) , "                " + templist[i] ) 
        tempdict["curl_files"] = templist
        templist = []
        print("")
        for i in range(len(temp['subs_dir'])):
            tempfile = tempdir + list(temp['subs_dir'])[i] + "/" + temp['tempfile'] + ".nc"
            templist.append(tempfile)
            print(" isfile =", os.path.isfile(templist[i]) , "                " + templist[i] ) 
        tempdict["temp_files"] = templist
        templist = []
        print("")
        for i in range(len(temp['subs_dir'])):
            tempfile = tempdir + list(temp['subs_dir'])[i] + "/" + temp['uvelfile'] + ".nc"
            templist.append(tempfile)
            print(" isfile =", os.path.isfile(templist[i]) , "                " + templist[i] ) 
        tempdict["uvel_files"] = templist
        templist = []
        print("")
        for i in range(len(temp['subs_dir'])):
            tempfile = tempdir + list(temp['subs_dir'])[i] + "/" + temp['vvelfile'] + ".nc"
            templist.append(tempfile)
            print(" isfile =", os.path.isfile(templist[i]) , "                " + templist[i] ) 
        tempdict["vvel_files"] = templist
        templist = []
        print("")
        for i in range(len(temp['subs_dir'])):
            tempfile = tempdir + list(temp['subs_dir'])[i] + "/" + temp['wvelfile'] + ".nc"
            templist.append(tempfile)
            print(" isfile =", os.path.isfile(templist[i]) , "                " + templist[i] ) 
        tempdict["wvel_files"] = templist
        LRfiles.append(tempdict)
        del temp
        print("")
    
    return HRfiles, LRfiles, IO_dict

