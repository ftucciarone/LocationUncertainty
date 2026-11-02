"""Compare the variables of two NetCDF files. Run: python tests/compare_nc.py new.nc reference.nc"""
import sys

import netCDF4
import numpy as np

with netCDF4.Dataset(sys.argv[1]) as a, netCDF4.Dataset(sys.argv[2]) as b:
    missing = set(b.variables) ^ set(a.variables)
    if missing:
        print("only in one file:", sorted(missing))
    for name in sorted(set(a.variables) & set(b.variables)):
        x, y = np.ma.filled(a[name][:], 0.0), np.ma.filled(b[name][:], 0.0)
        if x.shape != y.shape:
            print(f"{name:10s} shape {x.shape} vs {y.shape}")
            continue
        d = np.abs(x - y).max()
        scale = np.abs(y).max()
        corr = np.corrcoef(x.ravel(), y.ravel())[0, 1] if x.size > 1 and x.std() and y.std() else float("nan")
        print(f"{name:14s} corr {corr:.6f}, max|diff| {d:.3e}, relative {d/scale if scale else 0:.3e}, max|ref| {scale:.3e}")
