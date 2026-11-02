#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_noise_sigma.py

MPI-tiled production script: read GLORYS04's barotropic-removed daily u/v
(daily-baroclinic/, from build_baroclinic_daily.py), high-pass filter each
day IN PLACE per rank (spatial tile + halo read straight off disk -- no MPI
communication, see mpi_orca.read_slab_with_halo), fold the high-pass
residual into a running (Welford) joint mean/variance/covariance, and at the
end downsample the result onto the eORCA1 grid and write it as
noise_sigma.nc -- the file tlu_gss.F90 already
knows how to read (ln_tlu_stat / cn_tlu_stat_file), with EXACTLY the schema
it expects: u_tavg, v_tavg, u_tvar, v_tvar, uv_tvcov, all (z, y, x) on the
eORCA1 grid domain_cfg.nc describes.

Why a high-pass instead of the coarse-grain-then-subtract residual
--------------------------------------------------------------------
build_noise_BcResid_stats.py computes the same *kind* of quantity (the
sub-eORCA1-scale part of the flow) a different way: low-pass, bin-average
onto eORCA1, reconstruct on the fine grid, subtract. That extra
bin-average/reconstruct round trip mainly buys "the resolved part sits
exactly on eORCA1's own discretisation" -- a small conceptual refinement,
not different physics. A direct spatial high-pass (raw - low-pass, same
sigma) removes the same band of scales more simply, is the same operation
the SBC pipeline itself already uses elsewhere (see the noise-generation design notes
sec 4.2.3's H operator), and -- the reason this script exists -- factorises
cleanly into a per-tile, per-day, disk-halo-read operation with no
inter-rank communication at all, unlike the coarse-grain step (which needs
a global bin-average across tile boundaries). See the noise-generation design notes
sec 3.1 for the shared design background.

Two halo modes -- read-from-disk (day loop) vs exchange (post-processing)
---------------------------------------------------------------------------
Every rank owns one (nz, jl, il) horizontal tile, ALL 75 levels ("z-pencil",
same shape as mpi_filter.py), for the ENTIRE 2-year day loop -- the tiling is
decided once up front from the (static) grid size and never changes. Two
DIFFERENT spatial filters run in this script, and they get their halo two
different ways:

  1. Per-day HIGH-PASS (day loop, `highpass_tile` / `process_one_day`).
     Because the source files are read-only and every rank can open any one
     of them independently, the halo this Gaussian filter needs is read
     DIRECTLY OFF DISK (mpi_orca.read_slab_with_halo) instead of exchanged
     with neighbouring ranks after reading just the interior. This removes
     ALL per-day MPI traffic -- the day loop itself does zero communication.

  2. Post-processing anti-alias LOW-PASS (`exchange_lowpass`, once, after
     the day loop, before the eORCA1 downsample -- same purpose as
     build_noise_stat.py's own pre-downsample low-pass). This filter's
     input is the FINISHED Welford mean/variance/covariance -- a result
     that lives only in this rank's memory, computed by folding 365+ days
     together. There is no file to re-read a wider window from, so this
     step genuinely needs a halo EXCHANGE with neighbouring ranks
     (mpi_orca.add_halo / exchange_halo, the same primitives `mpi_filter.py`
     uses) -- the one and only place this script does per-tile MPI
     communication. Because it runs once (not once per day), the cost is
     irrelevant next to the day loop's.

Besides that one exchange (and the trailing gathers below), the only other
communication in the whole run is a handful of ONE-TIME gathers after the
day loop, to assemble the finished (nz, ny_g, nx_g) statistics on rank 0
for the downsample + write.

Memory: O(1) in the number of days (Welford), same discipline as
build_noise_BcResid_stats.py / build_noise_BcBtS_stats.py. Per rank, working
set is a handful of (nz, jl+2hw, il+2hw) tiles plus the Welford aggregate of
the same shape -- with, say, an 8x8 tiling (jl,il ~ 128, 180) and hw=32,
that is a few hundred MB, not the several GB a single full-domain rank
would need.

Run
---
    # one-time correctness check (serial vs an 8x8 tiling vs the whole-domain
    # reference smoother, one day only -- see check_tiling() below):
    mpirun -np 64 .../python build_noise_sigma.py --check

    # production: full 2-year archive, 8x8 = 64 ranks
    mpirun -np 64 .../python build_noise_sigma.py

    # production: one year only
    mpirun -np 64 .../python build_noise_sigma.py --year 2015

    # quick smoke test: a handful of days, any rank count (incl. serial)
    python build_noise_sigma.py --limit 3

Choosing the tile count: any perfect square (mpi_orca.Decomp.from_comm's
only requirement); pick P*P <= (free cores), and hw <= min(jl, il) is NOT
required here (read_slab_with_halo has no neighbour to run out of, unlike
add_halo) so a narrow tile with a wide halo is fine, just wasteful (a tile
narrower than hw re-reads a lot of overlap). 8x8 = 64 ranks against a
(1021, 1442) grid gives tiles ~ (128, 180) against hw = ceil(4*8) = 32 --
comfortable.

Runtime: high-pass filtering here is the SAME per-point cost as
build_noise_BcResid_stats.py's smoothing step (~9 s/component/day
full-domain, serial) -- but split across P*P tiles with no communication
overhead per day, so wall time should drop roughly by that factor (64x
tiling => minutes rather than hours for the full 2-year archive, network/
disk I/O rather than compute becomes the likely bottleneck at that point).
Not yet benchmarked at full scale; the --check run measures one day first.
"""

import argparse
import glob
import os
import re
import sys
import time

import numpy as np
import netCDF4 as netcdf

HERE = os.path.dirname(os.path.abspath(__file__))

from orcakit import paths

from orcakit import orca
from noisekit import welford
from orcakit.netcdf import nc_copy

from orcakit import mpi_orca as mo

try:
    from mpi4py import MPI
    _HAVE_MPI = True
except Exception:
    _HAVE_MPI = False

# -- Configuration ----------------------------------------------------------
acc_dtype = np.float32
database_dir = os.path.join(paths.get("glorys04"), "daily-baroclinic")
orca1_domcfg = paths.get("domain_orca1")
output_file = os.path.join(paths.get("outputs"), "noise_sigma.nc")

sigma_hp_default = 8.0     # Gaussian high-pass sigma, GLORYS04 gridpoints (matches
                            # build_noise_BcResid_stats.py's sigma_cg -- same scale cut)
sigma_lowpass_default = 12.0  # anti-alias low-pass before the eORCA1 downsample
                               # (matches build_noise_stat.py's own sigma_lowpass)

udate_re = re.compile(r"vozocrtx_glorys04_(\d{4})-(\d{2})-(\d{2})\.nc$")
vdate_re = re.compile(r"vomecrty_glorys04_(\d{4})-(\d{2})-(\d{2})\.nc$")


def discover_files(years=None, limit=None):
    r"""
    discover_files: the daily baroclinic u and v files, paired by date.

        Parameters:
            years (list): Years to keep, e.g. [2015], or None for all
            limit  (int): Keep only the first limit days, for tests; None for all

        Returns:
            ufiles (list): Sorted u files
            vfiles (list): v files of the same dates

        Formula:
            \{(u_d, v_d)\},\ \text{sorted, the first limit days}

        See also:
            main
    """
    ufiles = sorted(glob.glob(os.path.join(database_dir, "vozocrtx", "vozocrtx_glorys04_*.nc")))
    if years is not None:
        ufiles = [f for f in ufiles if int(udate_re.search(os.path.basename(f)).group(1)) in years]
    vfiles = []
    for uf in ufiles:
        y, m, d = udate_re.search(os.path.basename(uf)).groups()
        vf = os.path.join(database_dir, "vomecrty", f"vomecrty_glorys04_{y}-{m}-{d}.nc")
        if not os.path.exists(vf):
            raise RuntimeError(f"missing v-component partner for {uf}: {vf}")
        vfiles.append(vf)
    if limit is not None:
        ufiles, vfiles = ufiles[:limit], vfiles[:limit]
    return ufiles, vfiles


# ============================================================================
#  The per-tile, per-day operation: high-pass a field already read WITH its
#  real halo (mpi_orca.read_slab_with_halo) -- no exchange needed first.
# ============================================================================

def _masked_gaussian_padded(a, w, sigma):
    r"""
    _masked_gaussian_padded: land-aware Gaussian low-pass of an array that already carries its halo.

    The field times the mask and the mask are smoothed with the same kernel and divided, a weighted
    average over wet points; where the smoothed mask is below 0.01 the result is 0. orcakit.filters.lowpass
    with smooth=filters.gaussian does the same.

        Parameters:
            a   (float): (nz, jl+2hw, il+2hw) padded field
            w   (float): Same shape, wet mask
            sigma (float): Gaussian lengthscale, in grid points; levels are not mixed

        Returns:
            low (float): Low-pass, same padded shape

        Formula:
            \tilde f = \frac{G_\sigma * (f m)}{G_\sigma * m} \text{ where } G_\sigma * m > 0.01,\ \text{else } 0

        See also:
            orcakit.filters.lowpass
    """
    from scipy.ndimage import gaussian_filter
    sig = (0.0, sigma, sigma)
    num = gaussian_filter(a * w, sigma=sig, mode="nearest")
    den = gaussian_filter(w, sigma=sig, mode="nearest")
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0.01, num / den, 0.0)


def highpass_tile(a, w, sigma, hw, margin=0):
    r"""
    highpass_tile: land-aware Gaussian high-pass of a tile whose halo was read from disk.

    Pairs with mpi_orca.read_slab_with_halo, so there is no halo exchange. A margin of halo can be kept
    around the result, for lagged correlations with nearby points. orcakit.filters.highpass does the same
    on the whole grid.

        Parameters:
            a, w   (float): (nz, jl+2hw, il+2hw) padded field and wet mask
            sigma  (float): Gaussian lengthscale, in grid points
            hw       (int): Halo width, at least ceil(4 sigma); the caller makes sure of it
            margin   (int): Halo cells kept around the result, at most hw; 0 for the bare
                            interior. build_noise_lscale.py keeps kmax to correlate interior points
                            with their neighbours without a second read

        Returns:
            hp (float): (nz, jl+2margin, il+2margin) field minus low-pass, exactly 0 on land

        Formula:
            f' = (f - \tilde f)\, m

        See also:
            _masked_gaussian_padded, orcakit.filters.highpass
    """
    lowpass = _masked_gaussian_padded(a, w, sigma)
    jl, il = a.shape[-2] - 2 * hw, a.shape[-1] - 2 * hw
    lo_j, hi_j = hw - margin, hw + jl + margin
    lo_i, hi_i = hw - margin, hw + il + margin
    a_m = a[:, lo_j:hi_j, lo_i:hi_i]
    lp_m = lowpass[:, lo_j:hi_j, lo_i:hi_i]
    wet_m = w[:, lo_j:hi_j, lo_i:hi_i]
    hp = (a_m - lp_m) * wet_m
    return np.ascontiguousarray(hp, dtype=acc_dtype)


def exchange_lowpass(dec, field, wet, sigma, is_meridional, point="T"):
    r"""
    exchange_lowpass: anti-alias low-pass of a finished statistic, with its halo exchanged between ranks.

    A statistic summed over many days exists only in memory, so its halo comes from the neighbouring tiles
    rather than from disk. It is applied once, before the statistic is sampled onto a coarser grid.

        Parameters:
            dec        (Decomp): This rank's decomposition
            field, wet  (float): (nz, jl, il) statistic and wet mask
            sigma       (float): Gaussian lengthscale, in grid points
            is_meridional (bool): True for the means of u and v, which flip sign at the fold;
                                  False for variances and covariance
            point         (str): 'T', 'U', 'V' or 'F', where the field lives

        Returns:
            low (float): (nz, jl, il) low-passed statistic, 0 on land

        Formula:
            \tilde s = \frac{G_\sigma * (s m)}{G_\sigma * m}\, m

        See also:
            orcakit.filters.lowpass_tile
    """
    hw = int(np.ceil(4.0 * sigma))
    a = mo.add_halo(dec, field, hw, is_meridional=is_meridional, point=point)
    w = mo.add_halo(dec, wet, hw, is_meridional=False, point=point)
    lowpass = _masked_gaussian_padded(a, w, sigma)
    jl, il = field.shape[-2], field.shape[-1]
    lp_int = lowpass[:, hw:hw + jl, hw:hw + il]
    return np.ascontiguousarray(lp_int * wet, dtype=acc_dtype)


def process_one_day(dec, uf, vf, sigma, hw):
    r"""
    process_one_day: high-pass one day's u and v on this rank's tile, reading the halo from disk.

        Parameters:
            dec  (Decomp): This rank's decomposition
            uf, vf  (str): The day's daily-baroclinic files
            sigma (float): Gaussian lengthscale, in grid points
            hw      (int): Halo width, ceil(4 sigma)

        Returns:
            hu (float): (nz, jl, il) high-pass residual of u
            hv (float): (nz, jl, il) high-pass residual of v

        Formula:
            u' = \text{highpass}(u),\qquad v' = \text{highpass}(v)

        See also:
            highpass_tile
    """
    au, wu = mo.read_slab_with_halo(dec, uf, "vozocrtx", hw, is_meridional=True, point="U")
    av, wv = mo.read_slab_with_halo(dec, vf, "vomecrty", hw, is_meridional=True, point="V")
    hu = highpass_tile(au, wu, sigma, hw)
    hv = highpass_tile(av, wv, sigma, hw)
    return hu, hv


# ============================================================================
#  Correctness check: tiled result vs a whole-domain reference, one day
# ============================================================================

def reference_highpass(path, varname, sigma, is_meridional, point="T"):
    r"""
    reference_highpass: whole-domain high-pass in one process, the reference a tiled run is checked against.

    Uses orca.smooth_2D, an independent implementation of the boundary rules the tiled reading applies.

        Parameters:
            path          (str): NetCDF file
            varname       (str): Variable
            sigma       (float): Gaussian lengthscale, in grid points
            is_meridional (bool): True if the field's sign flips across the fold
            point         (str): 'T', 'U', 'V' or 'F', where the field lives

        Returns:
            hp (float): (nz, ny, nx) high-pass residual, 0 on land

        Formula:
            f' = (f - \tilde f)\, m,\qquad \tilde f = \frac{\text{smooth\_2D}(f m)}{\text{smooth\_2D}(m)}

        See also:
            check_tiling
    """
    with netcdf.Dataset(path) as nc:
        raw = np.ma.masked_invalid(nc[varname][:, :, :])
    data = np.ma.filled(raw, 0.0).astype(acc_dtype)
    wet = (~np.ma.getmaskarray(raw) & np.isfinite(data)).astype(acc_dtype)
    lowpass = orca.smooth_2D(data * wet, sigma, isMeridional=is_meridional, point=point)
    den = orca.smooth_2D(wet, sigma, isMeridional=False, point=point)
    with np.errstate(invalid="ignore", divide="ignore"):
        lowpass = np.where(den > 0.01, lowpass / den, 0.0)
    return np.ascontiguousarray((data - lowpass) * wet, dtype=acc_dtype)


def reference_lowpass(field_whole, wet_whole, sigma, is_meridional, point="T"):
    r"""
    reference_lowpass: whole-domain low-pass in one process, the reference for exchange_lowpass.

        Parameters:
            field_whole (float): (nz, ny, nx) field
            wet_whole   (float): (nz, ny, nx) wet mask
            sigma       (float): Gaussian lengthscale, in grid points
            is_meridional (bool): True if the field's sign flips across the fold
            point         (str): 'T', 'U', 'V' or 'F', where the field lives

        Returns:
            low (float): (nz, ny, nx) low-pass, 0 on land

        Formula:
            \tilde f = \frac{\text{smooth\_2D}(f m)}{\text{smooth\_2D}(m)}\, m

        See also:
            check_tiling
    """
    lowpass = orca.smooth_2D(field_whole * wet_whole, sigma, isMeridional=is_meridional, point=point)
    den = orca.smooth_2D(wet_whole, sigma, isMeridional=False, point=point)
    with np.errstate(invalid="ignore", divide="ignore"):
        lowpass = np.where(den > 0.01, lowpass / den, 0.0)
    return np.ascontiguousarray(lowpass * wet_whole, dtype=acc_dtype)


def check_tiling(dec, uf, vf, sigma, hw, sigma_lp):
    r"""
    check_tiling: check both halo modes on one day against whole-domain references.

    The tiled high-pass, with the halo read from disk, is compared with reference_highpass; that result,
    low-passed with the halo exchanged, is compared with reference_lowpass. One day stands in for a
    statistic: this checks the code, not any physical number.

        Parameters:
            dec    (Decomp): This run's decomposition, serial or tiled
            uf, vf    (str): One day's files
            sigma   (float): High-pass lengthscale, in grid points
            hw        (int): Halo width
            sigma_lp (float): Anti-alias low-pass lengthscale, in grid points

        Returns:
            None; prints the largest and rms difference per component and mode on rank 0

        Formula:
            \max |f_{\text{tiled}} - f_{\text{ref}}|,\quad \text{rms}(f_{\text{tiled}} - f_{\text{ref}})

        See also:
            reference_highpass, reference_lowpass
    """
    t0 = time.time()
    hu, hv = process_one_day(dec, uf, vf, sigma, hw)
    _, wu = mo.read_slab(dec, uf, "vozocrtx")
    _, wv = mo.read_slab(dec, vf, "vomecrty")
    gu = mo.gather(dec, hu)
    gv = mo.gather(dec, hv)
    if dec.rank == 0:
        dt = time.time() - t0
        ref_u = reference_highpass(uf, "vozocrtx", sigma, True, "U")
        ref_v = reference_highpass(vf, "vomecrty", sigma, True, "V")
        for name, g, ref in (("u", gu, ref_u), ("v", gv, ref_v)):
            e = np.abs(g - ref)
            print(f"[check hp {name}] tiling={dec.Py}x{dec.Px} hw={hw}  "
                  f"max|delta|={e.max():.3e}  rms={np.sqrt((e ** 2).mean()):.3e}  "
                  f"(dt={dt:.1f}s)")

    # ---- mode 2: exchange_lowpass on that same high-pass field ----
    t1 = time.time()
    lu = exchange_lowpass(dec, hu, wu, sigma_lp, is_meridional=True, point="U")
    lv = exchange_lowpass(dec, hv, wv, sigma_lp, is_meridional=True, point="V")
    glu = mo.gather(dec, lu)
    glv = mo.gather(dec, lv)
    wet_whole_u = mo.gather(dec, wu)   # collective -- every rank must call gather
    wet_whole_v = mo.gather(dec, wv)
    if dec.rank == 0:
        dt = time.time() - t1
        hw_lp = int(np.ceil(4.0 * sigma_lp))
        ref_lu = reference_lowpass(gu, wet_whole_u, sigma_lp, True, "U")
        ref_lv = reference_lowpass(gv, wet_whole_v, sigma_lp, True, "V")
        for name, g, ref in (("u", glu, ref_lu), ("v", glv, ref_lv)):
            e = np.abs(g - ref)
            print(f"[check lp {name}] tiling={dec.Py}x{dec.Px} hw_lp={hw_lp}  "
                  f"max|delta|={e.max():.3e}  rms={np.sqrt((e ** 2).mean()):.3e}  "
                  f"(dt={dt:.1f}s)")


# ============================================================================
#  Downsample onto eORCA1 and write noise_sigma.nc (rank 0 only, once)
# ============================================================================

def write_noise_sigma(lat_glorys, lon_glorys, wet_u, wet_v, mean_u, mean_v, var_u, var_v, cov_uv):
    r"""
    write_noise_sigma: average the statistics onto eORCA1 and write them as NEMO's noise statistics file.

    The file has u_tavg, v_tavg, u_tvar, v_tvar and uv_tvcov, all (z, y, x), as tlu_gss.F90 reads them.
    The wet masks are read once, from the first day.

        Parameters:
            lat_glorys, lon_glorys (float): GLORYS04 coordinates
            wet_u, wet_v (float): (nz, ny_g, nx_g) wet masks
            mean_u, mean_v, var_u, var_v, cov_uv (float): (nz, ny_g, nx_g) statistics

        Returns:
            None; writes output_file

        Formula:
            \bar s_c = \frac{\sum_{p \in c} m_p s_p}{\sum_{p \in c} m_p} \text{ for each statistic}

        See also:
            main
    """
    with netcdf.Dataset(orca1_domcfg) as nc:
        lon_orca1 = nc["nav_lon"][:, :]
        lat_orca1 = nc["nav_lat"][:, :]
        mask_orca1 = nc["top_level"][0, :, :] > 0
    n_coarse = lat_orca1.size
    nz = mean_u.shape[0]

    print("assigning every GLORYS04 point to its nearest eORCA1 cell ...")
    coarse_flat_idx = orca.assign_to_coarse_grid(lat_glorys, lon_glorys, lat_orca1, lon_orca1, mask_orca1)

    def to_orca1(field3d, weight3d):
        r"""
        to_orca1: average a 3D fine-grid field over the GLORYS04 points of each eORCA1 cell.

            Parameters:
                field3d  (float): (nz, ny_g, nx_g) fine-grid field
                weight3d (float): (nz, ny_g, nx_g) weights, the wet mask

            Returns:
                field (MaskedArray): (nz, 331, 360) on eORCA1, masked on land

            Formula:
                \bar f_c = \frac{\sum_{p \in c} w_p f_p}{\sum_{p \in c} w_p}

            See also:
                orcakit.orca.bin_to_coarse
        """
        avg, _ = orca.bin_to_coarse(field3d.astype(np.float64), weight3d,
                                     coarse_flat_idx, n_coarse)
        return np.ma.array(avg.reshape(nz, *lat_orca1.shape).astype(acc_dtype),
                            mask=np.broadcast_to(~mask_orca1, (nz,) + lat_orca1.shape))

    # uv_tvcov is weighted by wet_u -- cov_uv is only meaningful (and only
    # nonzero, by construction: read_slab_with_halo's fill-then-mask leaves
    # land at exactly 0) where u itself is wet; effectively identical to
    # wet_v away from the coastline these masks differ at.
    fields = {
        "u_tavg":   to_orca1(mean_u, wet_u),
        "v_tavg":   to_orca1(mean_v, wet_v),
        "u_tvar":   to_orca1(var_u,  wet_u),
        "v_tvar":   to_orca1(var_v,  wet_v),
        "uv_tvcov": to_orca1(cov_uv, wet_u),
    }
    print(f"downsampled to eORCA1 grid: {(nz,) + lat_orca1.shape}")

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    ncdest = nc_copy(orca1_domcfg, output_file, atts=False, dims=True, vars=False)
    for name, data in fields.items():
        var_id = ncdest.createVariable(name, "f4", ("z", "y", "x"), fill_value=0.0)
        var_id[...] = np.ma.filled(data, 0.0)
    ncdest.close()
    print(f"wrote {len(fields)} fields to {output_file}")


# ============================================================================
#  Driver
# ============================================================================

def main():
    r"""
    main: command-line entry point: statistics of the GLORYS04 small-scale baroclinic residual.

    The residual is the high-passed baroclinic velocity; its mean, variance and covariance, low-passed and
    averaged onto eORCA1, give the amplitude of the stochastic noise. Runs tiled under mpirun.

        Parameters:
            None; options are read from the command line (--sigma, --sigma-lowpass, --check and others)

        Returns:
            None; writes noise_sigma.nc

        Formula:
            u' = \text{highpass}(u_{bc}),\qquad \mu,\ \sigma^2,\ c_{uv} \text{ by Welford over all days}

        See also:
            process_one_day, exchange_lowpass, write_noise_sigma
    """
    ap = argparse.ArgumentParser(
        description="MPI-tiled GLORYS04 high-pass noise statistics -> noise_sigma.nc "
                    "(tlu_gss.F90's ln_tlu_stat / cn_tlu_stat_file schema)."
    )
    ap.add_argument("--sigma", type=float, default=sigma_hp_default,
                    help=f"Gaussian high-pass sigma, gridpoints (default {sigma_hp_default})")
    ap.add_argument("--sigma-lowpass", type=float, default=sigma_lowpass_default,
                    help=f"anti-alias low-pass sigma before the eORCA1 downsample, "
                         f"gridpoints (default {sigma_lowpass_default}, matches "
                         f"build_noise_stat.py's own pre-downsample low-pass); "
                         f"exchange-halo mode, see exchange_lowpass()")
    ap.add_argument("--year", type=int, action="append", default=None,
                    help="restrict to this year (repeatable, e.g. --year 2015 --year 2016); "
                         "default: the full archive")
    ap.add_argument("--limit", type=int, default=None,
                    help="process only the first N days found (smoke test)")
    ap.add_argument("--progress-every", type=int, default=10)
    ap.add_argument("--check", action="store_true",
                    help="one-day tiled-vs-reference correctness check instead of the full run")
    args = ap.parse_args()

    serial = (not _HAVE_MPI) or MPI.COMM_WORLD.Get_size() == 1
    ufiles, vfiles = discover_files(years=args.year, limit=args.limit)
    if not ufiles:
        raise RuntimeError("no daily-baroclinic files matched --year/--limit")

    dec = mo.Decomp.serial_decomp(*_grid_size(ufiles[0])) if serial \
        else mo.Decomp.from_comm(MPI.COMM_WORLD, *_grid_size(ufiles[0]))
    hw = int(np.ceil(4.0 * args.sigma))

    if dec.rank == 0:
        tag = "serial" if serial else f"{dec.Py}x{dec.Px} ({dec.nranks} ranks)"
        print(f"grid=({dec.ny_g},{dec.nx_g})  tiling={tag}  sigma={args.sigma}  hw={hw}  "
              f"{len(ufiles)} day(s): {os.path.basename(ufiles[0])} .. {os.path.basename(ufiles[-1])}")

    if args.check:
        check_tiling(dec, ufiles[0], vfiles[0], args.sigma, hw, args.sigma_lowpass)
        return

    # -------------------------------------------------------------- day loop
    with netcdf.Dataset(ufiles[0]) as nc:
        nz = nc["vozocrtx"].shape[0]
        lat_glorys = np.asarray(nc["nav_lat"][:, :]) if dec.rank == 0 else None
        lon_glorys = np.asarray(nc["nav_lon"][:, :]) if dec.rank == 0 else None
    # static wet masks, read once (interior only, no halo needed -- these
    # are never filtered themselves): land doesn't move day to day.
    _, wet_u_tile = mo.read_slab(dec, ufiles[0], "vozocrtx")
    _, wet_v_tile = mo.read_slab(dec, vfiles[0], "vomecrty")
    agg = welford.Welford_aggregate_xy_init(np.zeros((nz, dec.jl, dec.il), dtype=acc_dtype))

    t0 = time.time()
    for i, (uf, vf) in enumerate(zip(ufiles, vfiles)):
        hu, hv = process_one_day(dec, uf, vf, args.sigma, hw)
        agg = welford.Welford_covariance(agg, hu, hv)
        if dec.rank == 0 and ((i + 1) % args.progress_every == 0 or (i + 1) == len(ufiles)):
            print(f"  [{i + 1}/{len(ufiles)}] folded ({time.time() - t0:.0f}s elapsed)")

    mean_u, mean_v, var_u, var_v, cov_uv = welford.Welford_cov_finalize(agg)
    if dec.rank == 0:
        print(f"day loop done in {time.time() - t0:.0f}s over {len(ufiles)} days")

    # ----------------------------------------- anti-alias low-pass, once
    # (exchange-halo mode -- see exchange_lowpass()'s docstring for why this
    # step, unlike the day loop, needs real MPI communication). Means flip sign
    # at the fold, variances and covariance cannot.
    t1 = time.time()
    mean_u = exchange_lowpass(dec, mean_u, wet_u_tile, args.sigma_lowpass, is_meridional=True, point="U")
    mean_v = exchange_lowpass(dec, mean_v, wet_v_tile, args.sigma_lowpass, is_meridional=True, point="V")
    var_u = exchange_lowpass(dec, var_u, wet_u_tile, args.sigma_lowpass, is_meridional=False, point="U")
    var_v = exchange_lowpass(dec, var_v, wet_v_tile, args.sigma_lowpass, is_meridional=False, point="V")
    cov_uv = exchange_lowpass(dec, cov_uv, wet_u_tile, args.sigma_lowpass, is_meridional=False, point="U")
    dec.barrier()
    if dec.rank == 0:
        print(f"anti-alias low-pass (sigma={args.sigma_lowpass}) done in {time.time() - t1:.1f}s")

    # ------------------------------------------------- gather + write, once
    g_mean_u = mo.gather(dec, mean_u)
    g_mean_v = mo.gather(dec, mean_v)
    g_var_u = mo.gather(dec, var_u)
    g_var_v = mo.gather(dec, var_v)
    g_cov_uv = mo.gather(dec, cov_uv)
    g_wet_u = mo.gather(dec, wet_u_tile)
    g_wet_v = mo.gather(dec, wet_v_tile)

    if dec.rank == 0:
        wet2d = g_wet_u[0] > 0
        print(f"var_u surface domain-avg (wet) = {g_var_u[0][wet2d].mean():.3e} m2/s2")
        write_noise_sigma(lat_glorys, lon_glorys, g_wet_u, g_wet_v,
                          g_mean_u, g_mean_v, g_var_u, g_var_v, g_cov_uv)


def _grid_size(path):
    r"""
    _grid_size: the horizontal size of a daily file; u and v share it.

        Parameters:
            path (str): File holding vozocrtx

        Returns:
            ny (int): Rows
            nx (int): Columns

        Formula:
            (n_y, n_x) = \text{shape of vozocrtx}[1:]

        See also:
            main
    """
    with netcdf.Dataset(path) as nc:
        _, ny, nx = nc["vozocrtx"].shape
    return ny, nx


if __name__ == "__main__":
    main()
