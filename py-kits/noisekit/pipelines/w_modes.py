#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
W1: vertical-velocity modes by DMD, for tlu_gss.F90 (CMECS_wmodes.nc).

Computes w by integrating the horizontal divergence upward from the bottom,
decomposes it with exact DMD, keeps the modes near the unit circle and samples
them onto eORCA1.

Ported from the earlier CMEMS w-noise notebook. This deliberately differs
from that notebook, so CMECS_wmodes.nc cannot be reproduced with it:
  - the notebook differenced v[1:] - v[:-1] on a flipped grid, which gives the
    wrong sign of dv/dy, and did not use metric-weighted fluxes;
  - the notebook's wavg is 1.0 everywhere; here it is the mean of w.

    python -m noisekit.pipelines.w_modes --out w_modes.nc --start 2015-01 --end 2019-12
"""
import argparse
import importlib
import os

import netCDF4
import numpy as np

from datakit import clean
from modredkit import dmd, modes
from orcakit import paths
from orcakit import grid_c, netcdf, orca

DOMAIN_ORCA1 = paths.get("domain_orca1")
DT_MONTH = 31 * 86400


def snapshot_w(u, v, m):
    r"""
    snapshot_w: w from one pair of 3D velocities, without partial steps.

    Uses the 1D reference thicknesses and a divergence without thicknesses, so near the sea floor it
    differs from orcakit.grid_c.vertical_velocity.

        Parameters:
            u (float): (nz, ny, nx) U-point velocity, land as 0
            v (float): (nz, ny, nx) V-point velocity, land as 0
            m  (dict): e1t, e2t, e1v, e2u, e3t and tmask

        Returns:
            w (float): (nz, ny, nx) vertical velocity at the top of T-cells

        Formula:
            w_k = w_{k+1} - e_{3t,k}^0\, (\nabla_h \cdot \mathbf{u})_k\, t_k

        See also:
            orcakit.grid_c.vertical_velocity, orcakit.grid_c.w_from_divergence
    """
    div = grid_c.hDiv_Vgrid(u, v, m["e1t"], m["e2t"], m["e1v"], m["e2u"], orca.halo)
    return grid_c.w_from_divergence(div, m["e3t"], m["tmask"])


def run(out, start, end, dataset="cmems", rank=100.0, tol=0.01, map_file=None, domain=DOMAIN_ORCA1):
    r"""
    run: W1: DMD modes of the vertical velocity, as NEMO's noise input.

    w is computed from u and v, decomposed by exact DMD, and the modes near the unit circle are sampled
    onto eORCA1 and written in the format of tlu_gss.F90.

        Parameters:
            out      (str): File to write
            start    (str): First date
            end      (str): Last date
            dataset  (str): datakit module name
            rank   (float): Rank, or share of information to keep
            tol    (float): Largest distance of |eig| from 1 kept
            map_file (str): Index map cache, default <out dir>/<dataset>_to_eorca1.nc
            domain   (str): eORCA1 domain_cfg

        Returns:
            None; writes out

        Formula:
            w(t) \approx \bar w + \mathrm{Re}\sum_{|\,|\lambda_k| - 1| \le \text{tol}} b_k \phi_k\, \lambda_k^{\,t}

        See also:
            snapshot_w, modredkit.dmd.exactDMD
    """
    data = importlib.import_module(f"datakit.{dataset}")
    e1t, e2t, e1v, e2u, e3t_0, tmask = data.mesh("e1t", "e2t", "e1v", "e2u", "e3t_0", "tmask")
    m = {"e1t": e1t,
         "e2t": e2t,
         "e1v": e1v,
         "e2u": e2u,
         "tmask": tmask,
         "e3t": np.ones_like(tmask)*e3t_0[:, None, None]}
    ws = [snapshot_w(clean(u, fill=0.0).astype(np.float64), clean(v, fill=0.0).astype(np.float64), m)
          for _, u, v in data.stream("vozocrtx", "vomecrty", start=start, end=end, level=None)]
    w = np.stack(ws)
    del ws
    wet = np.ones(w.shape[1:], dtype=bool)
    X = modes.pack([w], [wet])
    mean = X.mean(axis=1, keepdims=True)
    A = X - mean
    mds, eig, _ = dmd.exactDMD(A[:, :-1], A[:, 1:], rank)
    ct = dmd.DMD_continuous_eig(eig, DT_MONTH)
    amp = dmd.DMD_amplitudes_pinv(mds, A[:, 0])
    keep = modes.mask_unit_circle(eig, tol)
    mds, ct, amp = modes.apply_mask(keep, mds, ct, amp)
    print(f"{keep.sum()} of {eig.size} modes kept", flush=True)

    with netCDF4.Dataset(domain) as nc:
        lat, lon = nc["nav_lat"][:, :], nc["nav_lon"][:, :]
        mask = nc["top_level"][0, :, :] > 0
    lat_f, lon_f = data.mesh("nav_lat", "nav_lon")
    map_file = map_file or os.path.join(os.path.dirname(os.path.abspath(out)), f"{dataset}_to_eorca1.nc")
    ilat, ilon = orca.load_or_make_map(map_file, lat, lon, lat_f, lon_f, mask=mask)
    wm = modes.unpack(mean, [wet])[0][0]
    wr = modes.unpack(mds*amp, [wet])[0]
    netcdf.write_modes(out, {"wavg": orca.downsample_field(wm, ilat, ilon, mask)},
                       {"wmode": orca.downsample_field(wr, ilat, ilon, mask)}, ct.imag)
    print("wrote", out, flush=True)


def main():
    r"""
    main: command-line entry point of W1.

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
    ap.add_argument("--start", required=True)
    ap.add_argument("--end", required=True)
    ap.add_argument("--dataset", default="cmems")
    ap.add_argument("--rank", type=float, default=100.0)
    ap.add_argument("--tol", type=float, default=0.01, help="distance from the unit circle")
    a = ap.parse_args()
    run(a.out, a.start, a.end, a.dataset, a.rank, a.tol)


if __name__ == "__main__":
    main()
