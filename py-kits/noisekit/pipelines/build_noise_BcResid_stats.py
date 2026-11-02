#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
prod/build_noise_BcResid_stats.py

Production script: the "barotropic-removed coarse-graining residual" noise
statistics for TLU. For every day of GLORYS04 u/v over the full 2-year
archive (731 days, gap-free), at the full 75-level column:

    1. read u, v in full
    2. split off the barotropic mode (e3-weighted vertical mean, via
       orca.split_barotropic + orca.e3uv_from_partial_steps) and KEEP the
       baroclinic part u_bc, v_bc
    3. coarse-grain each level of u_bc / v_bc onto the ORCA1 grid:
       Gaussian pre-smooth (sigma_cg gridpoints, representing ORCA1's own
       effective truncation) -> nearest-cell bin-average onto ORCA1 cells
       (orca.assign_to_coarse_grid / orca.bin_to_coarse) -> reconstruct on
       the fine grid.  The RESIDUAL = u_bc - that is the sub-ORCA1-scale
       baroclinic variability -- the target for the TLU sigma*dB_t term.
    4. fold the residual into a running Welford mean/variance accumulator
       (one day resident at a time -- 731 days x 75 levels will not stack)
    5. downsample the mean/variance fields onto the ORCA1 grid (bin-average
       the fine-grid statistic over each ORCA1 cell's fine points)
    6. write noise_BcResid_stats.nc (an nc_copy of the ORCA1 domain_cfg.nc,
       dimensions only, + the 4 fields ubc_resid_{mean,var}, vbc_resid_{mean,var})

Prototyped and step-by-step verified in glorys04-varnoise.ipynb (the
"Splitting GLORYS04 into an ORCA1-comparable resolved part" section and the
"Barotropic-removed residual" section); this is the same operation run over
the full 2 years at full depth.

sigma_cg = 8 GLORYS04 gridpoints: the *global* crossover between resolved
and residual variance found in the notebook. It is NOT uniform with
latitude (equatorial ~9, mid-latitude ~3-5, high-latitude ~7-13) -- a
latitude/space-dependent sigma is deliberately PARKED for now; see the
notebook and the tlu-regional-sigma-parked memory note.

DESEASONALIZE (below) is False: this computes the RAW residual
mean/variance. The notebook's year-long analysis flagged that an annual+
time-variance folds in the seasonal cycle (a predictable signal, not
noise) and should "almost certainly be deseasonalized first" -- doing that
properly needs a two-pass structure (pass 1: per-calendar-month
climatology of the residual; pass 2: Welford variance of the anomaly).
Not implemented here yet -- a flagged, obvious next step.

Memory: full-depth day ~450MB (float64) / ~225MB (float32) per component.
Welford aggregates (u/v residual, each mean+M2 at 75 x 1021 x 1442,
float32) ~1.8GB total; peak with per-day working arrays a few GB. NOT the
~150GB+ a stack-then-reduce would need.

Runtime: the batched 75-level smooth is ~9 s/component/day; ~23 s/day for
both -> order 4-5 hours over the full 2 years. A deliberate batch job, not
something to trigger by editing this file.
"""

import glob
import importlib
import os
import re
import sys
import time

import numpy as np
import netCDF4 as netcdf

from orcakit import paths
from orcakit import orca
from orcakit import filters
from noisekit import welford
from orcakit.netcdf import nc_copy
for _m in (orca, filters, welford):
    importlib.reload(_m)

from functools import partial

# -- Configuration ------------------------------------------------------------
acc_dtype = np.float32

database_dir = os.path.join(paths.get("glorys04"), "daily")
mesh_mask_file = os.path.join(paths.get("glorys04"), "static", "GLORYS2V4_mesh_mask.nc")
orca1_domcfg = paths.get("domain_orca1")
noise_stat_file = os.path.join(paths.get("outputs"), "noise_BcResid_stats.nc")

target_months = "all"       # full archive; or [(2015, 6)] etc. -- see glorys04-varnoise.ipynb
sigma_cg = 8.0              # coarse-graining pre-smooth cutoff, GLORYS04 gridpoints (see docstring)
DESEASONALIZE = False       # see docstring -- not implemented; raw residual statistics for now
progress_every = 10         # print a progress line every N days

udate_re = re.compile(r"vozocrtx_glorys04_(\d{4})-(\d{2})-(\d{2})\.nc$")
vdate_re = re.compile(r"vomecrty_glorys04_(\d{4})-(\d{2})-(\d{2})\.nc$")

if DESEASONALIZE:
    raise NotImplementedError("DESEASONALIZE is a flagged future two-pass extension; not built yet")


def discover_months(component_dir, varname):
    r"""
    discover_months: the months that have daily files of a variable.

        Parameters:
            component_dir (str): Folder of the variable's daily files
            varname       (str): Variable

        Returns:
            months (list): Sorted (year, month) pairs

        Formula:
            \{(y, m) : \exists\ v\_\text{glorys04}\_y\text{-}m\text{-}d.nc\}

        See also:
            bc_residual
    """
    pattern = re.compile(rf"{varname}_glorys04_(\d{{4}})-(\d{{2}})-\d{{2}}\.nc$")
    found = set()
    for f in glob.glob(os.path.join(component_dir, f"{varname}_glorys04_*.nc")):
        m = pattern.search(os.path.basename(f))
        if m:
            found.add((int(m.group(1)), int(m.group(2))))
    return sorted(found)


if target_months == "all":
    target_months = discover_months(os.path.join(database_dir, "vozocrtx"), "vozocrtx")
if not target_months:
    raise RuntimeError("target_months is empty -- nothing to load")

ufiles, vfiles = [], []
for yr, mo in target_months:
    ufiles += sorted(
        f for f in glob.glob(database_dir + f"/vozocrtx/vozocrtx_glorys04_{yr}-{mo:02d}-??.nc")
        if udate_re.search(os.path.basename(f))
    )
    vfiles += sorted(
        f for f in glob.glob(database_dir + f"/vomecrty/vomecrty_glorys04_{yr}-{mo:02d}-??.nc")
        if vdate_re.search(os.path.basename(f))
    )

if len(ufiles) != len(vfiles):
    raise RuntimeError(f"u/v file count mismatch: {len(ufiles)} vs {len(vfiles)}")
for uf, vf in zip(ufiles, vfiles):
    if udate_re.search(os.path.basename(uf)).groups() != vdate_re.search(os.path.basename(vf)).groups():
        raise RuntimeError(f"u/v date mismatch: {uf} vs {vf}")

n_days = len(ufiles)
print(f"{n_days} paired daily u/v snapshots, "
      f"{target_months[0][0]}-{target_months[0][1]:02d} .. {target_months[-1][0]}-{target_months[-1][1]:02d}")

# ============================================================================
# --- Static setup: e3u/e3v, the ORCA1 grid, the fine->coarse assignment ---
# ============================================================================
#
# e3u/e3v: NEMO partial-step reconstruction from mesh_mask.nc's T-only
# fields (orca.e3uv_from_partial_steps). coarse_flat_idx: for every
# GLORYS04 point, the nearest wet ORCA1 cell (orca.assign_to_coarse_grid,
# 3D-Cartesian nearest-neighbour -- fold/wrap-safe). All static.
#
with netcdf.Dataset(mesh_mask_file) as nc:
    e3t_0 = np.asarray(nc["e3t_0"][0, :])
    e3t_ps = np.asarray(nc["e3t_ps"][0, :, :])
    tmask_mm = np.asarray(nc["tmask"][0, :, :, :]).astype(bool)
    umask = np.asarray(nc["umask"][0, :, :, :]).astype(bool)
    vmask = np.asarray(nc["vmask"][0, :, :, :]).astype(bool)

e3u_3d, e3v_3d = orca.e3uv_from_partial_steps(e3t_0, e3t_ps, tmask_mm, umask, vmask)
e3u_3d = e3u_3d.astype(acc_dtype)
e3v_3d = e3v_3d.astype(acc_dtype)
umask_f = umask.astype(acc_dtype)
vmask_f = vmask.astype(acc_dtype)
nz = e3u_3d.shape[0]

with netcdf.Dataset(orca1_domcfg) as nc:
    lon_orca1 = nc["nav_lon"][:, :]
    lat_orca1 = nc["nav_lat"][:, :]
    mask_orca1 = nc["top_level"][0, :, :] > 0
n_coarse = lat_orca1.size

with netcdf.Dataset(ufiles[0]) as nc:
    lat_glorys = np.asarray(nc["nav_lat"][:, :])
    lon_glorys = np.asarray(nc["nav_lon"][:, :])

print("assigning every GLORYS04 point to its nearest ORCA1 cell ...")
coarse_flat_idx = orca.assign_to_coarse_grid(lat_glorys, lon_glorys, lat_orca1, lon_orca1, mask_orca1)
flat = coarse_flat_idx.ravel()

smooth_u = partial(orca.smooth_2D, mask=None, isMeridional=True, point="U")
smooth_v = partial(orca.smooth_2D, mask=None, isMeridional=True, point="V")
smooth_mask_u = partial(orca.smooth_2D, mask=None, isMeridional=False, point="U")
smooth_mask_v = partial(orca.smooth_2D, mask=None, isMeridional=False, point="V")
print(f"e3u/e3v, ORCA1 grid, fine->coarse map ready ({nz} levels)")

# ============================================================================
# --- MAIN STREAMING PASS -------------------------------------------------
# ============================================================================
#
# Per day: read full column -> strip barotropic (keep baroclinic) ->
# per-level coarse-grain -> residual -> fold into Welford. One day resident
# at a time.
#
t0 = time.time()
u_agg = welford.Welford_aggregate_init(np.zeros((nz,) + lat_glorys.shape, dtype=acc_dtype))
v_agg = welford.Welford_aggregate_init(np.zeros((nz,) + lat_glorys.shape, dtype=acc_dtype))


def bc_residual(path, varname, e3, wetmask_f, smoother, smoother_mask):
    r"""
    bc_residual: one day's baroclinic velocity minus its coarse-grained version, on the fine grid.

    The part of the baroclinic flow that a 1-degree model cannot represent: each level is smoothed,
    averaged over the fine points of each eORCA1 cell, put back on the fine grid and subtracted.

        Parameters:
            path        (str): Daily file
            varname     (str): 'vozocrtx' or 'vomecrty'
            e3        (float): (nz, ny, nx) thicknesses at the field's points
            wetmask_f (float): (nz, ny, nx) fine-grid wet mask
            smoother (callable): Smoothing of the field
            smoother_mask (callable): Smoothing of the mask, without sign flip

        Returns:
            residual (float): (nz, ny, nx) sub-eORCA1 residual, 0 on land

        Formula:
            r = \big(u' - \mathcal{P}\,\overline{G * u'}^{\,\text{cell}}\big)\, m,\qquad u' = u - \bar u^{z}

        See also:
            to_orca1, orcakit.orca.coarse_grain
    """
    with netcdf.Dataset(path) as nc:
        full = nc[varname][:, :, :].filled(0.0).astype(acc_dtype)
    _, bc = orca.split_barotropic(full, wetmask_f, e3=e3)          # (nz, y, x)
    del full
    smoothed = filters.masked_smooth(bc, wetmask_f, sigma_cg, smoother, smoother_mask)
    avg, _ = orca.bin_to_coarse(smoothed, wetmask_f, coarse_flat_idx, n_coarse)   # (nz, n_coarse)
    coarse_on_fine = avg[np.arange(nz)[:, None], flat[None, :]].reshape(bc.shape)
    return ((bc - coarse_on_fine) * wetmask_f).astype(acc_dtype)


for i, (uf, vf) in enumerate(zip(ufiles, vfiles)):
    u_agg = welford.Welford_variance(u_agg, bc_residual(uf, "vozocrtx", e3u_3d, umask_f, smooth_u, smooth_mask_u))
    v_agg = welford.Welford_variance(v_agg, bc_residual(vf, "vomecrty", e3v_3d, vmask_f, smooth_v, smooth_mask_v))
    if (i + 1) % progress_every == 0 or (i + 1) == n_days:
        print(f"  [{i + 1}/{n_days}] folded ({time.time() - t0:.0f}s elapsed)")

ubc_resid_mean, ubc_resid_var, _ = welford.Welford_finalize(u_agg)
vbc_resid_mean, vbc_resid_var, _ = welford.Welford_finalize(v_agg)
print(f"streaming pass done in {time.time() - t0:.0f}s over {n_days} days")

wet2d_u = umask.any(axis=0)
wet2d_v = vmask.any(axis=0)
print(f"ubc_resid_var  surface domain-avg  {ubc_resid_var[0][wet2d_u].mean():.3e} m2/s2")
print(f"vbc_resid_var  surface domain-avg  {vbc_resid_var[0][wet2d_v].mean():.3e} m2/s2")

# ============================================================================
# --- Downsample the statistic fields onto the ORCA1 grid ----------------
# ============================================================================
#
# Each ORCA1 cell gets the (wet-weighted) average of the fine-grid
# statistic over its own GLORYS04 points -- the natural coarse-grid
# representative of a sub-grid variance, using the SAME fine->coarse
# assignment as the split itself.
#
def to_orca1(field3d_fine, wetmask_f):
    r"""
    to_orca1: average a fine-grid statistic over the GLORYS04 points of each eORCA1 cell.

        Parameters:
            field3d_fine (float): (nz, ny, nx) fine-grid statistic
            wetmask_f    (float): (nz, ny, nx) fine-grid wet mask, the weights

        Returns:
            field (MaskedArray): (nz, 331, 360) on eORCA1, masked on land

        Formula:
            \bar f_c = \frac{\sum_{p \in c} m_p f_p}{\sum_{p \in c} m_p}

        See also:
            orcakit.orca.bin_to_coarse
    """
    avg, _ = orca.bin_to_coarse(field3d_fine.astype(np.float64), wetmask_f,
                                 coarse_flat_idx, n_coarse)          # (nz, n_coarse)
    ds = avg.reshape(nz, *lat_orca1.shape).astype(acc_dtype)
    return np.ma.array(ds, mask=np.broadcast_to(~mask_orca1, ds.shape))


fields = {
    "ubc_resid_mean": to_orca1(ubc_resid_mean, umask_f),
    "ubc_resid_var":  to_orca1(ubc_resid_var,  umask_f),
    "vbc_resid_mean": to_orca1(vbc_resid_mean, vmask_f),
    "vbc_resid_var":  to_orca1(vbc_resid_var,  vmask_f),
}
print(f"downsampled to ORCA1 grid: {(nz,) + lat_orca1.shape}")

# ============================================================================
# --- Write noise_BcResid_stats.nc --------------------------------------
# ============================================================================
#
# ORCA1's own "z" dimension is 75 levels and shares GLORYS04's deptht grid
# bit-for-bit (checked while building build_noise_BcBtS_stats.py), so the
# fields go straight onto ("z", "y", "x") -- horizontal downsample only,
# no vertical remapping.
#
os.makedirs(os.path.dirname(noise_stat_file), exist_ok=True)
ncdest = nc_copy(orca1_domcfg, noise_stat_file, atts=False, dims=True, vars=False)
for name, data in fields.items():
    var_id = ncdest.createVariable(name, "f4", ("z", "y", "x"), fill_value=0.0)
    var_id[...] = np.ma.filled(data, 0.0)
ncdest.close()
print(f"wrote {len(fields)} fields to {noise_stat_file}")
