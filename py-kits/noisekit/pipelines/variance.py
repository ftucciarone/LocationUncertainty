#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
V1, V2, V3: constant-variance noise statistics.

Optionally smooths each snapshot, keeps running mean, variance and covariance
of u and v, optionally low-passes them, samples them onto eORCA1 and writes them.

    V3 (GLORYS04 surface):  --dataset glorys04 --level 0 --lowpass 12 --names noise_stat
    V2 (filtered CMEMS):    --dataset cmems --smooth 4
    V1 (CMEMS):             --dataset cmems

--names picks the variable names: 'tlu' (u_tavg, ... as tlu_gss.F90 reads) or
'noise_stat' (umean, ... as build_noise_stat.py wrote).

The flux variant of V2 (download_and_process.py) is not covered.
"""
import argparse
import importlib
import os
import time
from functools import partial

import netCDF4
import numpy as np

from datakit import clean
from noisekit import welford
from orcakit import paths
from orcakit import filters, netcdf, orca

DOMAIN_ORCA1 = paths.get("domain_orca1")
NAMES = {"tlu": ("u_tavg", "v_tavg", "u_tvar", "v_tvar", "uv_tvcov"),
         "noise_stat": ("umean", "vmean", "uvar", "vvar", "uvvar")}


def run(out, dataset="cmems", start=None, end=None, level=None, smooth=0, lowpass=0, names="tlu", map_file=None, domain=DOMAIN_ORCA1):
    r"""
    run: V1 to V3: mean, variance and covariance of u and v, sampled onto eORCA1.

    Snapshots may be smoothed first and the statistics low-passed after. The covariance pairs U(i,j) with
    V(i,j), the points of the same index, as NEMO's noise input expects.

        Parameters:
            out      (str): File to write
            dataset  (str): datakit module name
            start    (str): First date, or None
            end      (str): Last date, or None
            level    (int): One level, or None for all
            smooth (float): Gaussian smoothing of each snapshot, in grid points; 0 for none
            lowpass (float): Gaussian low-pass of the statistics, in grid points; 0 for none
            names    (str): Key of NAMES, the variable names written
            map_file (str): Index map cache, default <out dir>/<dataset>_to_eorca1.nc
            domain   (str): eORCA1 domain_cfg

        Returns:
            None; writes out

        Formula:
            \mu,\ \sigma^2_u,\ \sigma^2_v,\ c_{uv} \text{ of } G_{\sigma_1} * u,\ G_{\sigma_1} * v,\ \text{then } \text{masked\_smooth}_{\sigma_2}

        See also:
            noisekit.welford.Welford_covariance, orcakit.filters.masked_smooth
    """
    t0 = time.time()
    data = importlib.import_module(f"datakit.{dataset}")
    agg = None
    u_land = v_land = None
    n = 0
    for date, u, v in data.stream("vozocrtx", "vomecrty", start=start, end=end, level=level):
        if agg is None:
            u_land, v_land = np.isnan(u), np.isnan(v)
            agg = welford.Welford_aggregate_xy_init(np.zeros(u.shape, dtype=np.float32))
        u, v = clean(u, fill=0.0).astype(np.float32), clean(v, fill=0.0).astype(np.float32)
        if smooth:
            u = orca.smooth_2D(u, smooth, isMeridional=True, point="U")
            v = orca.smooth_2D(v, smooth, isMeridional=True, point="V")
        agg = welford.Welford_covariance(agg, u, v)
        n += 1
        if n % 50 == 0:
            print(f"{date}: {n} snapshots, {time.time()-t0:.0f} s", flush=True)
    umean, vmean, uvar, vvar, cov = welford.Welford_cov_finalize(agg)
    u_wet, v_wet = (~u_land).astype(np.float32), (~v_land).astype(np.float32)

    if lowpass:
        def smoother(point, flip):
            r"""
            smoother: the ORCA Gaussian for one point and fold sign.

                Parameters:
                    point  (str): 'T', 'U', 'V' or 'F'
                    flip  (bool): True for a field that flips sign at the fold

                Returns:
                    smooth (callable): smooth(field, sigma)

                Formula:
                    f \mapsto \text{smooth\_2D}(f, \sigma, \text{point}, \text{flip})

                See also:
                    orcakit.orca.smooth_2D
            """
            return partial(orca.smooth_2D, mask=None, isMeridional=flip, point=point)
        umean = filters.masked_smooth(umean*u_wet, u_wet, lowpass, smoother("U", True), smoother("U", False))
        vmean = filters.masked_smooth(vmean*v_wet, v_wet, lowpass, smoother("V", True), smoother("V", False))
        uvar = filters.masked_smooth(uvar*u_wet, u_wet, lowpass, smoother("U", False))
        vvar = filters.masked_smooth(vvar*v_wet, v_wet, lowpass, smoother("V", False))
        # u'v' pairs a U and a V point, half a cell apart: T is the nearest fold mapping
        cov = filters.masked_smooth(cov*u_wet*v_wet, u_wet*v_wet, lowpass, smoother("T", False))

    with netCDF4.Dataset(domain) as nc:
        lat, lon = nc["nav_lat"][:, :], nc["nav_lon"][:, :]
        mask = nc["top_level"][0, :, :] > 0
    with netCDF4.Dataset(data.files("vozocrtx")[0]) as nc:
        lat_f, lon_f = nc["nav_lat"][:, :], nc["nav_lon"][:, :]
    map_file = map_file or os.path.join(os.path.dirname(os.path.abspath(out)), f"{dataset}_to_eorca1.nc")
    ilat, ilon = orca.load_or_make_map(map_file, lat, lon, lat_f, lon_f, mask=mask)
    fields = {}
    for name, f in zip(NAMES[names], (umean, vmean, uvar, vvar, cov)):
        ds = orca.downsample_field(f, ilat, ilon, mask)
        fields[name] = np.ma.array(ds, mask=np.broadcast_to(~mask, ds.shape))
    netcdf.write_fields(domain, out, fields)
    print(f"wrote {out}, {n} snapshots, {time.time()-t0:.0f} s", flush=True)


def main():
    r"""
    main: command-line entry point of V1 to V3.

        Parameters:
            None; options are read from the command line

        Returns:
            None

        Formula:
            \text{run}(\text{options})

        See also:
            run
    """
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--dataset", default="cmems")
    ap.add_argument("--start")
    ap.add_argument("--end")
    ap.add_argument("--level", type=int, help="one vertical level, default all")
    ap.add_argument("--smooth", type=float, default=0, help="smooth each snapshot, sigma in gridpoints")
    ap.add_argument("--lowpass", type=float, default=0, help="low-pass the statistics, sigma in gridpoints")
    ap.add_argument("--names", choices=sorted(NAMES), default="tlu")
    a = ap.parse_args()
    run(a.out, a.dataset, a.start, a.end, a.level, a.smooth, a.lowpass, a.names)


if __name__ == "__main__":
    main()
