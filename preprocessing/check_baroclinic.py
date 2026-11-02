#!/usr/bin/env python3
"""
Check the output of baroclinic.py against what the split must satisfy, for any source:
the baroclinic part has zero e3-weighted depth mean, the two parts rebuild the velocity,
and land is NaN where the mesh says so and only there.

    python check_baroclinic.py --source glorys12 --dst /path/to/out --start 1993-01 --end 1993-01
"""
import argparse
import os
import sys

import numpy as np

import datakit
from datakit import clean


def check_date(src, out, var, date, path, e3, mask, tol):
    """Errors of one variable on one date, and the number of wet points without a valid source value."""
    raw = src.snapshot(path, var, np.float64)
    u = clean(raw, fill=0.0)
    invalid = int((mask & ~(np.abs(np.nan_to_num(raw, nan=1e30)) <= 10)).sum())
    dry_values = int((~mask & np.isfinite(raw) & (raw != 0)).sum())
    written = out.groups(var, start=date, end=date)[0][0][1]
    bc = out.snapshot(written, var, np.float64)
    bt = out.snapshot(written, datakit.BAROTROPIC[var], np.float64)
    wet2d = mask.any(axis=0)
    land_ok = (np.isnan(bc) == ~mask).all() and (np.isnan(bt) == ~wet2d).all()
    weight = e3.astype(np.float64)*mask
    mean = np.sum(np.where(mask, bc, 0.0)*weight, axis=0)/np.where(wet2d, weight.sum(axis=0), 1.0)
    rebuild = np.abs(np.where(mask, bc + np.where(wet2d, bt, 0.0)[None] - u, 0.0)).max()
    print(f"  {date} {var}: |depth mean| {np.abs(mean).max():.2e}, |rebuild| {rebuild:.2e}, "
          f"max |baroclinic| {np.nanmax(np.abs(bc)):.2f}, max |barotropic| {np.nanmax(np.abs(bt)):.2f}, land ok {land_ok}, wet without value {invalid}, dry with value {dry_values}", flush=True)
    return land_ok and np.abs(mean).max() < tol and rebuild < tol


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source", required=True)
    ap.add_argument("--dst", required=True)
    ap.add_argument("--start")
    ap.add_argument("--end")
    ap.add_argument("--tol", type=float, default=1e-4)
    a = ap.parse_args()
    src = datakit.source(a.source)
    out = datakit.derived(os.path.join(a.dst, "split"), src)
    e3u, e3v, umask, vmask = src.vertical(np.float64)
    metric = {"vozocrtx": (e3u, umask),
              "vomecrty": (e3v, vmask)}
    print(f"wet points without thickness: U {(umask & (e3u <= 0)).sum()}, V {(vmask & (e3v <= 0)).sum()}", flush=True)
    good = True
    for date, *paths in src.groups(*metric, start=a.start, end=a.end)[0]:
        for var, path in zip(metric, paths):
            good &= check_date(src, out, var, date, path, *metric[var], a.tol)
    print("ALL OK" if good else "FAILED")
    sys.exit(0 if good else 1)


if __name__ == "__main__":
    main()
