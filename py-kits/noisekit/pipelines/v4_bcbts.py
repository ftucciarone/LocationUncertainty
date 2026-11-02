#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
V4: barotropic/baroclinic split noise statistics from daily GLORYS04.

Keeps running statistics of the barotropic and baroclinic parts of u and v,
low-passes them, samples them onto the eORCA1 grid and writes them. to_tlu()
converts the result to the file tlu_gss.F90 reads.

Ported from the earlier scripts build_noise_BcBtS_stats.py and
build_TLUstat_from_BcBtS.py.

    python -m noisekit.pipelines.v4_bcbts --out stats.nc --tlu tlu_stats.nc
"""
import argparse
import os
import time
from functools import partial

import netCDF4
import numpy as np

from datakit import clean, glorys04
from noisekit import welford
from orcakit import paths
from orcakit import filters, netcdf, orca

DOMAIN_ORCA1 = paths.get("domain_orca1")
SIGMA_LOWPASS = 12.0


def accumulate(start=None, end=None, progress_every=20):
    r"""
    accumulate: running mean and variance of the barotropic and baroclinic parts of daily GLORYS04 u and v.

        Parameters:
            start          (str): First date, or None
            end            (str): Last date, or None
            progress_every (int): Print progress every this many days

        Returns:
            stats (dict): {name: (mean, var)} for ubt, vbt, ubc and vbc
            umask (bool): (nz, ny, nx) U-point wet mask
            vmask (bool): (nz, ny, nx) V-point wet mask

        Formula:
            \mu,\ \sigma^2 \text{ of } \bar u^{z},\ \bar v^{z},\ u',\ v' \text{ by Welford over all days}

        See also:
            orcakit.orca.split_barotropic, noisekit.welford.Welford_variance
    """
    e3t_0, e3t_ps, tmask, umask, vmask = glorys04.mesh("e3t_0", "e3t_ps", "tmask", "umask", "vmask")
    tmask, umask, vmask = tmask.astype(bool), umask.astype(bool), vmask.astype(bool)
    e3u, e3v = orca.e3uv_from_partial_steps(e3t_0, e3t_ps, tmask, umask, vmask)
    e3u, e3v = e3u.astype(np.float32), e3v.astype(np.float32)
    umask_f, vmask_f = umask.astype(np.float32), vmask.astype(np.float32)

    agg = {"ubt": welford.Welford_aggregate_init(np.zeros(umask.shape[1:], dtype=np.float32)),
           "vbt": welford.Welford_aggregate_init(np.zeros(vmask.shape[1:], dtype=np.float32)),
           "ubc": welford.Welford_aggregate_init(np.zeros(umask.shape, dtype=np.float32)),
           "vbc": welford.Welford_aggregate_init(np.zeros(vmask.shape, dtype=np.float32))}
    t0 = time.time()
    n = 0
    for date, u, v in glorys04.stream("vozocrtx", "vomecrty", start=start, end=end, level=None):
        u = clean(u, limit=np.inf, fill=0.0).astype(np.float32)
        v = clean(v, limit=np.inf, fill=0.0).astype(np.float32)
        ubt, ubc = orca.split_barotropic(u, umask_f, e3=e3u)
        vbt, vbc = orca.split_barotropic(v, vmask_f, e3=e3v)
        for name, part in (("ubt", ubt), ("vbt", vbt), ("ubc", ubc), ("vbc", vbc)):
            agg[name] = welford.Welford_variance(agg[name], part)
        n += 1
        if n % progress_every == 0:
            print(f"{date}: {n} days, {time.time()-t0:.0f} s", flush=True)
    return {k: welford.Welford_finalize(a)[:2] for k, a in agg.items()}, umask, vmask


def lowpass(stats, umask, vmask, sigma=SIGMA_LOWPASS):
    r"""
    lowpass: mask-aware low-pass of the statistics.

    Means flip sign at the north fold; variances and masks do not.

        Parameters:
            stats (dict): From accumulate
            umask (bool): U-point wet mask
            vmask (bool): V-point wet mask
            sigma (float): Gaussian lengthscale, in grid points

        Returns:
            fields (dict): {name_mean, name_var: field}, 0 on land

        Formula:
            \tilde s = \text{masked\_smooth}(s),\ \text{sign at the fold} = -1 \text{ for means},\ +1 \text{ for variances}

        See also:
            orcakit.filters.masked_smooth
    """
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
    wet = {"ubt": umask.any(axis=0).astype(np.float32),
           "vbt": vmask.any(axis=0).astype(np.float32),
           "ubc": umask.astype(np.float32),
           "vbc": vmask.astype(np.float32)}
    out = {}
    for name, (mean, var) in stats.items():
        w = wet[name]
        point = name[0].upper()
        plain = smoother(point, False)
        out[f"{name}_mean"] = filters.masked_smooth(mean*w, w, sigma, smoother(point, True), plain)
        out[f"{name}_var"] = filters.masked_smooth(var*w, w, sigma, plain)
    return out


def to_orca1(fields, map_file, domain=DOMAIN_ORCA1):
    r"""
    to_orca1: sample fields onto eORCA1 at the nearest GLORYS04 point.

        Parameters:
            fields  (dict): {name: (nz, ny, nx) or (ny, nx) field}
            map_file (str): Index map cache, see orca.load_or_make_map
            domain   (str): eORCA1 domain_cfg

        Returns:
            out (dict): {name: MaskedArray on eORCA1}, masked on land

        Formula:
            f^c_{j,i} = f^f_{\,\text{ilat}_{j,i},\,\text{ilon}_{j,i}}

        See also:
            orcakit.orca.downsample_field
    """
    with netCDF4.Dataset(domain) as nc:
        lat, lon = nc["nav_lat"][:, :], nc["nav_lon"][:, :]
        mask = nc["top_level"][0, :, :] > 0
    with netCDF4.Dataset(glorys04.files("vozocrtx")[0]) as nc:
        lat_f, lon_f = nc["nav_lat"][:, :], nc["nav_lon"][:, :]
    ilat, ilon = orca.load_or_make_map(map_file, lat, lon, lat_f, lon_f, mask=mask)
    out = {}
    for name, f in fields.items():
        nz = f.shape[0] if f.ndim == 3 else None
        ds = orca.downsample_field(f, ilat, ilon, mask, nz=nz)
        out[name] = np.ma.array(ds, mask=np.broadcast_to(~mask, ds.shape))
    return out


def to_tlu(stats_file, tlu_file):
    r"""
    to_tlu: write NEMO's noise statistics file from the V4 statistics, baroclinic part only.

    The file has u_tavg, v_tavg, u_tvar, v_tvar and uv_tvcov. uv_tvcov is 0: it is not computed here, and
    tlu_gss.F90 does not use it.

        Parameters:
            stats_file (str): Statistics written by main
            tlu_file   (str): File to write

        Returns:
            None

        Formula:
            (u, v)_{tavg} = \mu_{bc},\quad (u, v)_{tvar} = \sigma^2_{bc},\quad uv_{tvcov} = 0

        See also:
            main
    """
    with netCDF4.Dataset(stats_file) as src, netCDF4.Dataset(tlu_file, "w") as dst:
        shape = src["ubc_mean"].shape
        for name, size in zip(("z", "y", "x"), shape):
            dst.createDimension(name, size)
        dst.source = f"noisekit.pipelines.v4_bcbts, baroclinic half of {os.path.basename(stats_file)}"
        for new, old in (("u_tavg", "ubc_mean"), ("v_tavg", "vbc_mean"), ("u_tvar", "ubc_var"), ("v_tvar", "vbc_var")):
            dst.createVariable(new, "f8", ("z", "y", "x"), fill_value=0.0)[:] = src[old][:]
        dst.createVariable("uv_tvcov", "f8", ("z", "y", "x"), fill_value=0.0)[:] = np.zeros(shape)


def main():
    r"""
    main: command-line entry point of V4: statistics of the barotropic and baroclinic parts, on eORCA1.

        Parameters:
            None; options are read from the command line (--out, --tlu, --map, --start, --end)

        Returns:
            None; writes --out and, with --tlu, the tlu_gss.F90 file

        Formula:
            \text{to\_orca1}(\text{lowpass}(\text{accumulate}()))

        See also:
            accumulate, lowpass, to_orca1, to_tlu
    """
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tlu", help="also write the tlu_gss.F90 file here")
    ap.add_argument("--map", help="index map cache, default <out dir>/glorys04_to_eorca1.nc")
    ap.add_argument("--start")
    ap.add_argument("--end")
    a = ap.parse_args()
    map_file = a.map or os.path.join(os.path.dirname(os.path.abspath(a.out)), "glorys04_to_eorca1.nc")
    stats, umask, vmask = accumulate(a.start, a.end)
    fields = to_orca1(lowpass(stats, umask, vmask), map_file)
    netcdf.write_fields(DOMAIN_ORCA1, a.out, fields)
    if a.tlu:
        to_tlu(a.out, a.tlu)
    print("wrote", a.out, flush=True)


if __name__ == "__main__":
    main()
