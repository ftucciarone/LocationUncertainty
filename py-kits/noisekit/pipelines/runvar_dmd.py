#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
D1 and D2: time-varying noise amplitude from DMD of a running variance.

Smooths velocity, takes its variance (or standard deviation) over a sliding window of
snapshots, decomposes the anomalies with exact or kernel DMD, keeps one mode per
conjugate pair, samples onto eORCA1 and writes noise_DMDstatistics.nc for tlu_gss.F90.

Ported from the earlier CMEMS running-variance and kernel DMD notebooks.
Differences from the notebooks:
  - no grid flip, and the whole grid is read, not the [-2::-1, 1:-1] crop;
  - snapshots stay in time order. The notebooks start with i = -whs, which
    stores the first window in the last column;
  - the tapers are computed here, not read from CMECS_to_eOrca1.nc;
  - the default statistic is variance, as in the backup notebook, which
    produced CMECS_varDMD_2016-2017.nc. The later notebook switched to std.

    python -m noisekit.pipelines.runvar_dmd --out noise_DMDstatistics.nc --start 2016-01 --end 2016-12
"""
import argparse
import importlib
import os
import time

import netCDF4
import numpy as np

from datakit import clean
from modredkit import dmd, modes
from noisekit.running import rolling_stats
from orcakit import paths
from orcakit import netcdf, orca

DOMAIN_ORCA1 = paths.get("domain_orca1")
DT_MONTH = 31 * 86400
POINT = {"vozocrtx": "U",
         "vomecrty": "V"}


def _months_before(start, n):
    r"""
    _months_before: the month n months before a given one.

        Parameters:
            start (str): 'YYYY-MM' or 'YYYY-MM-DD'
            n     (int): Months back

        Returns:
            month (str): 'YYYY-MM'

        Formula:
            (y, m) - n \text{ months}

        See also:
            window_stat
    """
    y, m = map(int, start.split("-")[:2])
    total = y*12 + (m-1) - n
    return f"{total//12:04d}-{total%12+1:02d}"


def window_stat(data, var, start, end, window, sigma, taper, stat="var", level=None):
    r"""
    window_stat: running variance or standard deviation of a smoothed, tapered velocity.

    How strongly a field varies itself changes in time, with the seasons for example. Each window's
    statistic follows that change. The record starts window - 1 months before start, so the first window
    ends at start.

        Parameters:
            data   (module): datakit module, e.g. datakit.cmems
            var       (str): 'vozocrtx' or 'vomecrty'
            start     (str): Month of the first window's end
            end       (str): Month of the last window's end
            window    (int): Snapshots per window
            sigma   (float): Gaussian smoothing of each snapshot, in grid points; 0 for none
            taper   (float): Coastal taper the snapshots are multiplied by
            stat      (str): 'var' or 'std'
            level     (int): One level, or None for all

        Returns:
            stat (float): (nt, [nz,] ny, nx), one per window

        Formula:
            s_t = \mathrm{var}\big(\tau\, G_\sigma * u\big) \text{ over } t - W + 1, \dots, t

        See also:
            noisekit.running.rolling_stats, run
    """
    first = _months_before(start, window-1)

    def snapshots():
        r"""
        snapshots: the cleaned, smoothed and tapered snapshots, one at a time.

            Parameters:
                None

            Returns:
                generator: (ny, nx) or (nz, ny, nx) arrays, land as 0

            Formula:
                \tau\, G_\sigma * u_t,\ \text{land as } 0

            See also:
                window_stat
        """
        for _, f in data.stream(var, start=first, end=end, level=level):
            f = clean(f, fill=0.0).astype(np.float64)
            if sigma:
                f = orca.smooth_2D(f, sigma, isMeridional=True, point=POINT.get(var, "T"))
            yield f * taper
    return np.stack([std**2 if stat == "var" else std for _, _, std in rolling_stats(snapshots(), window)])


def decompose(X, rank, kernel, dt):
    r"""
    decompose: exact or kernel DMD of a matrix of anomalies.

        Parameters:
            X    (float): (n, m) anomalies, one snapshot per column
            rank (float): Rank, or share of information to keep
            kernel (str): None for exact DMD, or 'linear' or 'poly' for kernel DMD (see modes.gram)
            dt   (float): Time between snapshots, in seconds

        Returns:
            mds (complex): Modes
            eig (complex): Discrete-time eigenvalues
            ct  (complex): Continuous-time eigenvalues

        Formula:
            \text{exactDMD}(A_{:,:-1}, A_{:,1:}) \text{ or } \text{kernelDMD}(\text{gram}(A))

        See also:
            modredkit.dmd.exactDMD, modredkit.dmd.kernelDMD
    """
    if kernel is None:
        mds, eig, _ = dmd.exactDMD(X[:, :-1], X[:, 1:], rank)
        return mds, eig, dmd.DMD_continuous_eig(eig, dt)
    S2 = modes.gram(X, kernel)
    mds, ct, eig, *_ = dmd.kernelDMD(S2, X[:, :-1].T, dt, rank)
    return mds, eig, ct


def run(out, start, end, dataset="cmems", window=3, rank=100.0, kernel=None, sigma=4, map_file=None, domain=DOMAIN_ORCA1, nmodes=None, stat="var", level=None):
    r"""
    run: D1 and D2: modes of the running variance of u and v, as NEMO's noise input.

    The running variance follows how variability changes in time. Its DMD gives a few oscillating spatial
    patterns; one mode per conjugate pair is kept, and the mean and modes are tapered at the coasts,
    sampled onto eORCA1 and written in the format of tlu_gss.F90. D1 uses exact DMD, D2 kernel DMD.

        Parameters:
            out       (str): File to write
            start     (str): Month of the first window's end
            end       (str): Month of the last window's end
            dataset   (str): datakit module name
            window    (int): Snapshots per window
            rank    (float): Rank, or share of information to keep
            kernel    (str): None for D1, 'linear' or 'poly' for D2
            sigma   (float): Gaussian smoothing of each snapshot, in grid points
            map_file  (str): Index map cache, default <out dir>/<dataset>_to_eorca1.nc
            domain    (str): eORCA1 domain_cfg
            nmodes    (int): Keep this many modes with the largest energy, or None for all
            stat      (str): 'var' or 'std'
            level     (int): One level, or None for all

        Returns:
            None; writes out

        Formula:
            s(t) \approx \bar s + \mathrm{Re}\sum_k 2\, b_k \phi_k\, e^{\lambda^c_k t}

        See also:
            window_stat, decompose, orcakit.netcdf.write_modes
    """
    t0 = time.time()
    data = importlib.import_module(f"datakit.{dataset}")
    tmask = data.mesh("tmask").astype(np.float32)
    taper_f = orca.coastal_taper(tmask, 3, 8)
    with netCDF4.Dataset(domain) as nc:
        lat, lon = nc["nav_lat"][:, :], nc["nav_lon"][:, :]
        top, bot = nc["top_level"][0].data, nc["bottom_level"][0].data
    k = np.arange(1, tmask.shape[0]+1)[:, None, None]
    tmask_c = ((k >= top[None]) & (k <= bot[None])).astype(np.float32)
    taper_c = orca.coastal_taper(tmask_c, 1.5, 4)
    mask_c = tmask_c[0] > 0
    lat_f, lon_f = data.mesh("nav_lat", "nav_lon")
    map_file = map_file or os.path.join(os.path.dirname(os.path.abspath(out)), f"{dataset}_to_eorca1.nc")
    ilat, ilon = orca.load_or_make_map(map_file, lat, lon, lat_f, lon_f, mask=mask_c)
    print(f"grids ready, {time.time()-t0:.0f} s", flush=True)

    taper = taper_f if level is None else taper_f[level]
    ustd = window_stat(data, "vozocrtx", start, end, window, sigma, taper, stat, level)
    vstd = window_stat(data, "vomecrty", start, end, window, sigma, taper, stat, level)
    print(f"{ustd.shape[0]} windows, {time.time()-t0:.0f} s", flush=True)
    wet = np.ones(ustd.shape[1:], dtype=bool)
    X = modes.pack([ustd, vstd], [wet, wet])
    mean = X.mean(axis=1, keepdims=True)
    mds, eig, ct = decompose(X - mean, rank, kernel, DT_MONTH)
    amp = dmd.DMD_amplitudes_pinv(mds, (X - mean)[:, 0])
    print(f"{eig.size} modes, {time.time()-t0:.0f} s", flush=True)

    # Lowest frequency first: tlu_gss.F90 reads only the first modes of the file
    mds, eig, ct, amp = modes.apply_mask(modes.mask_one_per_pair(ct), mds, eig, ct, amp)
    order = np.argsort(np.abs(ct.imag))
    if nmodes:
        order = order[:nmodes]
    mds, ct, amp = mds[:, order], ct[order], amp[order]

    (um, vm) = modes.unpack(mean, [wet, wet])
    (ur, vr) = modes.unpack(mds*amp, [wet, wet])

    def coarse(f):
        r"""
        coarse: taper, sample onto eORCA1 and taper again.

            Parameters:
                f (float): Fine-grid field, all levels, or one level when run's level is set

            Returns:
                field (float): On eORCA1, all jpk levels; with one level, only that level is filled

            Formula:
                \tau_c \cdot \text{downsample}(\tau_f\, f)

            See also:
                orcakit.orca.downsample_field
        """
        ds = orca.downsample_field(f*taper, ilat, ilon, mask_c)
        if level is None:
            return ds * taper_c
        full = np.zeros(ds.shape[:-2] + (tmask_c.shape[0],) + ds.shape[-2:], dtype=ds.dtype)
        full[..., level, :, :] = ds * taper_c[level]
        return full

    netcdf.write_modes(out, {"uavg": coarse(um[0]),
                             "vavg": coarse(vm[0])},
                       {"umode": coarse(ur),
                        "vmode": coarse(vr)}, ct.imag)
    print(f"wrote {out}, {time.time()-t0:.0f} s", flush=True)


def main():
    r"""
    main: command-line entry point of D1 and D2.

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
    ap.add_argument("--start", required=True, help="YYYY-MM of the first window's end")
    ap.add_argument("--end", required=True)
    ap.add_argument("--dataset", default="cmems")
    ap.add_argument("--window", type=int, default=3)
    ap.add_argument("--rank", type=float, default=100.0, help="an information percentage")
    ap.add_argument("--kernel", choices=("linear", "poly"), help="use kernel DMD (D2)")
    ap.add_argument("--sigma", type=float, default=4)
    ap.add_argument("--nmodes", type=int)
    ap.add_argument("--stat", choices=("var", "std"), default="var")
    ap.add_argument("--level", type=int, help="one vertical level, default all")
    a = ap.parse_args()
    run(a.out, a.start, a.end, a.dataset, a.window, a.rank, a.kernel, a.sigma, nmodes=a.nmodes, stat=a.stat, level=a.level)


if __name__ == "__main__":
    main()
