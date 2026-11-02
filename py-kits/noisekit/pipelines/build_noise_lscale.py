#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_noise_lscale.py

MPI-tiled production script: measure the REAL spatial correlation length of
GLORYS04's barotropic-removed, high-passed velocity residual (the same
residual build_noise_sigma.py computes the amplitude of) and write it as
noise_lscale.nc -- the length-scale/anisotropy counterpart of noise_sigma.nc,
per the noise-generation design notes, sec 4.2.7 ("OFFLINE -- build_noise_lscale.py").

Method: time as the ensemble dimension, at each FIXED spatial lag
-----------------------------------------------------------------
The design notes, sec 4.1, measured Lx/Ly by pooling a 2-point spatial
correlation over a region box -- a few large boxes, not a per-gridpoint map.
A naive per-gridpoint version of that (a small spatial sliding window
centred on each point) would be poorly sampled (~2*halfwidth independent
values per estimate) and would mix together whatever physical regimes fall
inside that one window.

Instead: for a fixed spatial lag k, the quantity
    rho_x(k; i,j) = < u'(i,j,t) u'(i+k,j,t) >_t  /  < u'(i,j,t)^2 >_t
is a genuine per-gridpoint statistic, averaged over the WHOLE archive (365+
independent days) at that one fixed pair of points -- far better sampled,
and it reuses build_noise_sigma.py's exact day-loop/halo-read machinery
almost unchanged. The halo already read for the Gaussian high-pass
(hw=ceil(4*sigma_hp)=32) is *more* than the lag range this script needs
(kmax~25) -- so lagged correlation costs no extra disk I/O or MPI traffic
beyond what build_noise_sigma.py already pays per day; see
highpass_tile(..., margin=kmax) (a small additive extension to that
function, backward compatible -- margin=0 there still means exactly what it
always meant).

u' and v' are pooled together (not treated as two separate length-scale
maps): the actual thing being measured is the spatial scale of the
underlying eddy STREAMFUNCTION (the design notes' sec 4.2 curl-of-
streamfunction construction) -- both u' and v' are projections of the SAME
scalar field, so pooling their lagged products before dividing (rather than
computing Lx from u' and v' separately and averaging two already-independent
ratios) is the statistically coherent thing to do, not just a convenience.

2-D, not 3-D: this is the length scale of the CONSTRUCTION operator's
smoothing, applied to a 2-D streamfunction -- unlike noise_sigma.nc's
u_tvar/v_tvar (a 3-D amplitude field), the design notes' own pseudocode
for w_x/w_y (sec 4.2.7 step 2) has no z index. So this script only ever
touches the SURFACE level (mpi_orca.read_slab_with_halo(..., level=0)) --
75x cheaper than build_noise_sigma.py's full-depth day loop, which is why
this can comfortably run over the full 2-year archive rather than one year.

Known limitation (documented, not fixed): a coastal wet cell whose k-lag
neighbour is on land gets rho_x/rho_y artificially pulled toward 0 at that
lag (the neighbour's high-pass residual is exactly 0 there, not "unknown"),
biasing the measured length scale short in a narrow coastal band. A properly
wet-pair-weighted estimator (numerator/denominator accumulated separately,
masked_smooth-style) would fix this; not implemented here -- the bias is
narrow-band and this project has accepted comparable edge simplifications
elsewhere (e.g. Step G/Step H's Neumann-not-fold boundary in
the Chebyshev-Jacobi prototype notebook).
"""

import argparse
import os
import sys
import time

import numpy as np
import netCDF4 as netcdf

HERE = os.path.dirname(os.path.abspath(__file__))

from orcakit import paths
from orcakit import orca
from orcakit.netcdf import nc_copy

from orcakit import mpi_orca as mo
# reuse the already-validated Gaussian high-pass + file discovery, rather
# than re-implementing them -- see build_noise_sigma.py's own module
# docstring for why the halo-read (not exchange) mode is used for the day
# loop, and highpass_tile's `margin` parameter (added for this script,
# backward compatible: margin=0 there is unchanged from before).
from noisekit.pipelines import build_noise_sigma as bns

try:
    from mpi4py import MPI
    _HAVE_MPI = True
except Exception:
    _HAVE_MPI = False

# -- Configuration ----------------------------------------------------------
acc_dtype = np.float32
mesh_mask_file = os.path.join(paths.get("glorys04"), "static", "GLORYS2V4_mesh_mask.nc")
orca1_domcfg = paths.get("domain_orca1")
output_file = os.path.join(paths.get("outputs"), "noise_lscale.nc")

sigma_hp_default = 8.0   # SAME scale-cut as build_noise_sigma.py -- "residual" must mean the same thing in both files
kmax_default = 25        # max lag, GLORYS04 gridpoints; must be < ceil(4*sigma_hp) (32 at the default sigma) --
                          # comfortably above the ~16 gp equatorial Lx the design notes' sec 4.1 already measured
r_clamp = 5.0             # sanity clamp on the anisotropy ratio (either direction) -- pathological-outlier guard,
                          # not a physical ceiling (cf. Step H's r_max=3, a DIFFERENT, tighter, eORCA1-grid-specific
                          # cap on a downstream-consumed quantity; this is just noise-floor protection on the raw
                          # measurement, values above this reflect near-zero-variance estimation noise, not eddies)


# ============================================================================
#  Per-day, per-tile: lagged products, pooling u' and v'
# ============================================================================

def accumulate_lag_products(hu, hv, kmax, Px, Py, P0_pad):
    r"""
    accumulate_lag_products: add one day's high-pass residual to the running sums of lagged products, in place.

    u' and v' are pooled into the same sums. The sum of squares covers the whole padded tile, not only
    its interior, so each lag can be normalised by the variances at both its ends: normalising by the
    reference point alone gave correlations up to 27, since a correlation is bounded by 1 only with the
    symmetric form.

        Parameters:
            hu, hv   (float): (1, jl+2kmax, il+2kmax) residual, with a kmax-wide margin
            kmax       (int): Largest lag, in grid points
            Px, Py   (float): (kmax, jl, il) running sums of lags 1 to kmax in x and y
            P0_pad   (float): (jl+2kmax, il+2kmax) running sum of u'^2 + v'^2

        Returns:
            None; Px, Py and P0_pad are updated in place

        Formula:
            P_x(k) \mathrel{+}= u'_{i}\,u'_{i+k} + v'_{i}\,v'_{i+k},\qquad P_y(k) \text{ likewise in } j,\qquad P_0 \mathrel{+}= u'^2 + v'^2

        See also:
            main, noisekit.correlation.spatial_correlation
    """
    hu0 = hu[0]
    hv0 = hv[0]
    jl = hu0.shape[0] - 2 * kmax
    il = hu0.shape[1] - 2 * kmax

    P0_pad += hu0 * hu0 + hv0 * hv0    # self product over the WHOLE padded tile, not just the interior

    hu_c = hu0[kmax:kmax + jl, kmax:kmax + il]
    hv_c = hv0[kmax:kmax + jl, kmax:kmax + il]

    for k in range(1, kmax + 1):
        hu_xk = hu0[kmax:kmax + jl, kmax + k:kmax + k + il]
        hv_xk = hv0[kmax:kmax + jl, kmax + k:kmax + k + il]
        Px[k - 1] += hu_c * hu_xk + hv_c * hv_xk

        hu_yk = hu0[kmax + k:kmax + k + jl, kmax:kmax + il]
        hv_yk = hv0[kmax + k:kmax + k + jl, kmax:kmax + il]
        Py[k - 1] += hu_c * hu_yk + hv_c * hv_yk


# ============================================================================
#  Vectorised e-folding lag, over the whole (jl, il) map at once
# ============================================================================

def efold_map(rho, kmax):
    r"""
    efold_map: the lag at which each point's correlation first drops below 1/e.

    Linear interpolation between the two lags around the crossing; lag 0 is 1. Where it never drops below
    1/e within kmax, the length is at least kmax: it is returned as kmax and flagged as censored.

        Parameters:
            rho  (float): (kmax, jl, il) correlation at lags 1 to kmax
            kmax   (int): Largest lag

        Returns:
            L        (float): (jl, il) e-folding lag, in grid points; exactly kmax where censored
            censored  (bool): (jl, il) True where rho never dropped below 1/e: the true length is
                              at least kmax

        Formula:
            L = (k-1) + \frac{\rho_{k-1} - e^{-1}}{\rho_{k-1} - \rho_k},\qquad k = \min\{k : \rho_k < e^{-1}\}

        See also:
            noisekit.correlation.efolding_lag
    """
    below = rho < (1.0 / np.e)                    # (kmax, jl, il)
    censored = ~below.any(axis=0)
    first_k = np.argmax(below, axis=0) + 1         # +1: rho index 0 here IS lag 1 (lag 0 excluded, see docstring)
    # bracket [k-1, k]; rho at k-1 is 1.0 (the implicit lag-0 value) when first_k==1
    rho_k = np.take_along_axis(rho, (first_k - 1)[None, :, :], axis=0)[0]
    rho_km1 = np.where(first_k == 1, 1.0,
                        np.take_along_axis(rho, (first_k - 2).clip(min=0)[None, :, :], axis=0)[0])
    denom = rho_km1 - rho_k
    frac = np.where(denom > 1e-12, (rho_km1 - 1.0 / np.e) / denom, 0.0)
    L = np.where(censored, float(kmax), (first_k - 1) + frac)
    return L.astype(acc_dtype), censored


# ============================================================================
#  Driver
# ============================================================================

def main():
    r"""
    main: command-line entry point: measure the spatial correlation length of the GLORYS04 small-scale residual.

    The residual is the high-passed baroclinic velocity. Its lagged correlations in x and y over time give
    the e-folding lengths and their ratio, written as the length-scale input of the stochastic noise.
    Runs tiled under mpirun; its --sigma must match the one used for the noise amplitude.

        Parameters:
            None; options are read from the command line (--sigma, which must match
            build_noise_sigma.py's, --kmax, --check and others)

        Returns:
            None; writes noise_lscale.nc

        Formula:
            \rho_x(k) = \frac{P_x(k)}{\sqrt{P_0(i)\,P_0(i+k)}},\qquad L_x = \text{efold}(\rho_x)\, e_{1t}

        See also:
            accumulate_lag_products, efold_map, write_noise_lscale
    """
    ap = argparse.ArgumentParser(
        description="MPI-tiled GLORYS04 spatial-correlation-length measurement -> noise_lscale.nc "
                    "(the length-scale/anisotropy counterpart of noise_sigma.nc)."
    )
    ap.add_argument("--sigma", type=float, default=sigma_hp_default,
                    help=f"Gaussian high-pass sigma, gridpoints (default {sigma_hp_default}, "
                         f"MUST match build_noise_sigma.py's --sigma for 'the residual' to mean the same thing)")
    ap.add_argument("--kmax", type=int, default=kmax_default,
                    help=f"max correlation lag, GLORYS04 gridpoints (default {kmax_default})")
    ap.add_argument("--year", type=int, action="append", default=None,
                    help="restrict to this year (repeatable); default: the full archive")
    ap.add_argument("--limit", type=int, default=None,
                    help="process only the first N days found (smoke test)")
    ap.add_argument("--progress-every", type=int, default=20)
    ap.add_argument("--check", action="store_true",
                    help="tiled-vs-serial tiling-invariance check on a short window instead of the full run")
    args = ap.parse_args()

    serial = (not _HAVE_MPI) or MPI.COMM_WORLD.Get_size() == 1
    ufiles, vfiles = bns.discover_files(years=args.year, limit=args.limit)
    if not ufiles:
        raise RuntimeError("no daily-baroclinic files matched --year/--limit")

    hw = int(np.ceil(4.0 * args.sigma))
    if args.kmax >= hw:
        raise ValueError(f"--kmax ({args.kmax}) must be < the high-pass halo width hw=ceil(4*sigma)={hw}, "
                          f"or the lag correlation would reach past the halo read_slab_with_halo provides")

    dec = mo.Decomp.serial_decomp(*bns._grid_size(ufiles[0])) if serial \
        else mo.Decomp.from_comm(MPI.COMM_WORLD, *bns._grid_size(ufiles[0]))

    if dec.rank == 0:
        tag = "serial" if serial else f"{dec.Py}x{dec.Px} ({dec.nranks} ranks)"
        print(f"grid=({dec.ny_g},{dec.nx_g})  tiling={tag}  sigma={args.sigma}  hw={hw}  kmax={args.kmax}  "
              f"{len(ufiles)} day(s): {os.path.basename(ufiles[0])} .. {os.path.basename(ufiles[-1])}")

    # static wet mask (surface only), read once -- land doesn't move day to day
    _, wet_u_tile = mo.read_slab(dec, ufiles[0], "vozocrtx", level=0)

    Px = np.zeros((args.kmax, dec.jl, dec.il), dtype=np.float64)   # float64 accumulator: many days summed
    Py = np.zeros((args.kmax, dec.jl, dec.il), dtype=np.float64)
    P0_pad = np.zeros((dec.jl + 2 * args.kmax, dec.il + 2 * args.kmax), dtype=np.float64)

    t0 = time.time()
    for i, (uf, vf) in enumerate(zip(ufiles, vfiles)):
        au, wu = mo.read_slab_with_halo(dec, uf, "vozocrtx", hw, is_meridional=True, level=0, point="U")
        av, wv = mo.read_slab_with_halo(dec, vf, "vomecrty", hw, is_meridional=True, level=0, point="V")
        hu = bns.highpass_tile(au, wu, args.sigma, hw, margin=args.kmax)
        hv = bns.highpass_tile(av, wv, args.sigma, hw, margin=args.kmax)
        accumulate_lag_products(hu, hv, args.kmax, Px, Py, P0_pad)
        if dec.rank == 0 and ((i + 1) % args.progress_every == 0 or (i + 1) == len(ufiles)):
            print(f"  [{i + 1}/{len(ufiles)}] folded ({time.time() - t0:.0f}s elapsed)")

    n_days = len(ufiles)
    Px /= n_days
    Py /= n_days
    P0_pad /= n_days
    if dec.rank == 0:
        print(f"day loop done in {time.time() - t0:.0f}s over {n_days} days")

    # Symmetric (Pearson) normalisation: sqrt(var at the reference point *
    # var at the point k away), NOT the reference point's variance alone --
    # see accumulate_lag_products' docstring for why the asymmetric version
    # is wrong (violates |rho|<=1 near sharp variance gradients, e.g. a
    # boundary current). P0_int is P0_pad's own interior (= lag-0 "self"
    # variance at the reference point); P0_shift_x/_y are P0_pad shifted by
    # k, read directly out of the SAME padded accumulator (already covers
    # the margin, no extra halo or I/O needed).
    kmax = args.kmax
    jl, il = dec.jl, dec.il
    P0_int = P0_pad[kmax:kmax + jl, kmax:kmax + il]
    rho_x = np.empty((kmax, jl, il), dtype=acc_dtype)
    rho_y = np.empty((kmax, jl, il), dtype=acc_dtype)
    with np.errstate(invalid="ignore", divide="ignore"):
        for k in range(1, kmax + 1):
            P0_shift_x = P0_pad[kmax:kmax + jl, kmax + k:kmax + k + il]
            denom_x = np.sqrt(P0_int * P0_shift_x)
            rho_x[k - 1] = np.where(denom_x > 1e-20, Px[k - 1] / denom_x, 0.0)

            P0_shift_y = P0_pad[kmax + k:kmax + k + jl, kmax:kmax + il]
            denom_y = np.sqrt(P0_int * P0_shift_y)
            rho_y[k - 1] = np.where(denom_y > 1e-20, Py[k - 1] / denom_y, 0.0)

    if args.check:
        check_tiling(dec, rho_x, rho_y, args.kmax)
        return

    Lx_gp, censored_x = efold_map(rho_x, args.kmax)
    Ly_gp, censored_y = efold_map(rho_y, args.kmax)

    # every gather below is collective -- call unconditionally on every
    # rank (gather itself returns None on non-root, per its own docstring)
    g_censored_x = mo.gather(dec, censored_x.astype(acc_dtype))
    g_censored_y = mo.gather(dec, censored_y.astype(acc_dtype))
    g_Lx_gp = mo.gather(dec, Lx_gp)
    g_Ly_gp = mo.gather(dec, Ly_gp)
    g_wet = mo.gather(dec, wet_u_tile)

    if dec.rank == 0:
        wet2d = g_wet[0] > 0
        pct_cx = 100 * g_censored_x[wet2d].mean()
        pct_cy = 100 * g_censored_y[wet2d].mean()
        print(f"censored (never dropped below 1/e within kmax={args.kmax}, wet cells only): "
              f"x {pct_cx:.1f}%, y {pct_cy:.1f}%")
        write_noise_lscale(g_Lx_gp, g_Ly_gp, g_wet)


# ============================================================================
#  Correctness check: tiling-invariance (serial vs tiled must agree exactly)
# ============================================================================

def check_tiling(dec, rho_x, rho_y, kmax):
    r"""
    check_tiling: print a summary of the lag correlations on rank 0, to compare a serial and a tiled run.

    There is no independent implementation of this statistic, so this checks that the result does not
    depend on the tiling, not that it is right.

        Parameters:
            dec      (Decomp): This run's decomposition
            rho_x, rho_y (float): (kmax, jl, il) this rank's lag correlations
            kmax        (int): Largest lag

        Returns:
            None; prints the minimum, mean and maximum of rho_x and rho_y at lags 1 and kmax

        Formula:
            \min,\ \text{mean},\ \max \text{ of } \rho_x, \rho_y \text{ at lags } 1 \text{ and } k_{max}

        See also:
            main
    """
    g_rho_x = mo.gather(dec, rho_x)
    g_rho_y = mo.gather(dec, rho_y)
    if dec.rank == 0:
        for name, g in (("rho_x", g_rho_x), ("rho_y", g_rho_y)):
            print(f"[check {name}] tiling={dec.Py}x{dec.Px}  lag1: min={g[0].min():.4f} mean={g[0].mean():.4f} "
                  f"max={g[0].max():.4f}   lag{kmax}: min={g[-1].min():.4f} mean={g[-1].mean():.4f} "
                  f"max={g[-1].max():.4f}")


# ============================================================================
#  Physical units, downsample onto eORCA1, write noise_lscale.nc
# ============================================================================

def write_noise_lscale(Lx_gp, Ly_gp, wet):
    r"""
    write_noise_lscale: write the correlation lengths on eORCA1.

    The lags become distances with GLORYS04's cell widths, are clamped, averaged onto eORCA1, and written
    as Lx, Ly in metres and their ratio, all 2D.

        Parameters:
            Lx_gp, Ly_gp (float): (ny_g, nx_g) e-folding lags, in grid points
            wet          (float): (1, ny_g, nx_g) wet mask, from read_slab

        Returns:
            None; writes output_file

        Formula:
            L_x = \bar{(k_x\, e_{1t})}^{\,\text{cell}},\quad L_y = \bar{(k_y\, e_{2t})}^{\,\text{cell}},\quad r = L_x / L_y

        See also:
            efold_map
    """
    with netcdf.Dataset(mesh_mask_file) as nc:
        e1t = np.asarray(nc["e1t"][0])
        e2t = np.asarray(nc["e2t"][0])
        lat_glorys = np.asarray(nc["nav_lat"][:, :])
        lon_glorys = np.asarray(nc["nav_lon"][:, :])

    Lx_phys = Lx_gp * e1t
    Ly_phys = Ly_gp * e2t
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.where(Ly_phys > 1.0, Lx_phys / Ly_phys, 1.0)
    r = np.clip(r, 1.0 / r_clamp, r_clamp)
    print(f"physical length scale: Lx {Lx_phys[wet[0] > 0].mean() / 1e3:.0f} km, "
          f"Ly {Ly_phys[wet[0] > 0].mean() / 1e3:.0f} km (domain-mean, wet cells), "
          f"r domain-mean {r[wet[0] > 0].mean():.2f} (clamped to [{1 / r_clamp:.2f}, {r_clamp:.2f}])")

    with netcdf.Dataset(orca1_domcfg) as nc:
        lon_orca1 = nc["nav_lon"][:, :]
        lat_orca1 = nc["nav_lat"][:, :]
        mask_orca1 = nc["top_level"][0, :, :] > 0
    n_coarse = lat_orca1.size

    print("assigning every GLORYS04 point to its nearest eORCA1 cell ...")
    coarse_flat_idx = orca.assign_to_coarse_grid(lat_glorys, lon_glorys, lat_orca1, lon_orca1, mask_orca1)
    wet2d = wet[0]

    def to_orca1(field2d):
        r"""
        to_orca1: average a 2D fine-grid field over the GLORYS04 points of each eORCA1 cell.

            Parameters:
                field2d (float): (ny_g, nx_g) fine-grid field

            Returns:
                field (MaskedArray): (331, 360) on eORCA1, masked on land

            Formula:
                \bar f_c = \frac{\sum_{p \in c} m_p f_p}{\sum_{p \in c} m_p}

            See also:
                orcakit.orca.bin_to_coarse
        """
        avg, _ = orca.bin_to_coarse(field2d.astype(np.float64), wet2d, coarse_flat_idx, n_coarse)
        return np.ma.array(avg.reshape(lat_orca1.shape).astype(acc_dtype), mask=~mask_orca1)

    Lx_ds = to_orca1(Lx_phys)
    Ly_ds = to_orca1(Ly_phys)
    # r recomputed on the COARSE fields, not itself downsampled -- averaging a
    # ratio directly would weight it incorrectly relative to averaging its
    # numerator and denominator separately (the same reasoning as u'/v' pooling
    # in accumulate_lag_products: combine before dividing, not after)
    with np.errstate(invalid="ignore", divide="ignore"):
        r_ds = np.ma.where(Ly_ds > 1.0, Lx_ds / Ly_ds, 1.0)
    r_ds = np.clip(r_ds, 1.0 / r_clamp, r_clamp)
    print(f"downsampled to eORCA1 grid: {lat_orca1.shape}")

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    ncdest = nc_copy(orca1_domcfg, output_file, atts=False, dims=True, vars=False)
    fields = {"Lx": Lx_ds,
              "Ly": Ly_ds,
              "r": r_ds}
    for name, data in fields.items():
        var_id = ncdest.createVariable(name, "f4", ("y", "x"), fill_value=0.0)
        var_id[...] = np.ma.filled(data, 0.0)
    ncdest.close()
    print(f"wrote {len(fields)} fields to {output_file}")


if __name__ == "__main__":
    main()
