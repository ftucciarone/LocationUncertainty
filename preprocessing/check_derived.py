#!/usr/bin/env python3
"""
Check the output of vertical_velocity.py (--step w) or geostrophic.py (--step geostrophic), for any source.

Both: land is NaN where the mesh says so, and the seam and north fold hold.
w: compared with the source's own w when it has one, level by level.
geostrophic: compared with the surface velocity of the source, by latitude band, when the dates overlap.

    python check_derived.py --step w --source glorys12 --dst /path/to/out --start 1993-01 --end 1993-01
"""
import argparse
import os
import sys

import numpy as np

import datakit
from orcakit import orca


def duplicates(field, point, sign):
    """Largest difference across the east-west seam and the north fold. Fold rows follow orca.FOLD."""
    nx = field.shape[-1]
    row, _ = orca.FOLD[point]
    seam = np.nanmax(np.abs(field[..., 0] - field[..., nx - 2]))
    ny = field.shape[-2]
    mirrored = field[..., ny - row + 1, :][..., orca.fold_columns(np.arange(nx), nx, point)]
    return max(seam, np.nanmax(np.abs(field[..., ny - 1, :] - sign*mirrored)))


def correlation(a, b):
    good = np.isfinite(a) & np.isfinite(b)
    return np.corrcoef(a[good], b[good])[0, 1], np.sqrt(np.mean((a[good] - b[good])**2)), np.sqrt(np.mean(b[good]**2))


def check_w(src, out, date, tol):
    written = out.groups("w", start=date, end=date)[0][0][1]
    w = out.snapshot(written, "w", np.float64)
    tmask = src.mesh("tmask").astype(bool)
    land = bool((np.isnan(w) == ~tmask).all())
    dup = duplicates(w[0], "T", 1)
    print(f"  {date}: land ok {land}, seam and fold {dup:.2e}, surface rms {np.nanstd(w[0]):.2e} max {np.nanmax(np.abs(w[0])):.2e} m/s", flush=True)
    try:
        raw = src.snapshot(src.groups("w", start=date, end=date)[0][0][1], "w", np.float64)
    except (FileNotFoundError, IndexError):
        print("  source has no w", flush=True)
        return land and dup < tol
    lat = src.grid()["lat"]
    for k in (0, 5, 15, 25, 35, 45):
        if k < w.shape[0]:
            r, rmse, rms = correlation(w[k], raw[k])
            band = np.abs(lat) < 60
            rb, *_ = correlation(w[k][band], raw[k][band])
            print(f"    level {k}: corr {r:.3f} (|lat|<60: {rb:.3f}), rms difference {rmse:.2e}, rms source {rms:.2e}", flush=True)
    return land and dup < tol


def check_geostrophic(src, out, date, min_lat, tol):
    grid_lat, gphiu, gphiv, umask, vmask = src.grid()["lat"], *src.mesh("gphiu", "gphiv", "umask", "vmask")
    good = True
    for var, name, lat, mask, point in (("vozocrtx", "u", gphiu, umask[0], "U"), ("vomecrty", "v", gphiv, vmask[0], "V")):
        g = out.snapshot(out.groups(var, start=date, end=date)[0][0][1], var, np.float64)
        wet = mask.astype(bool) & (np.abs(lat) >= min_lat)
        wet[0] = False
        extra, absent = int((np.isfinite(g) & ~mask.astype(bool)).sum()), int((np.isnan(g) & wet).sum())
        band = int(np.isfinite(g[np.abs(lat) < min_lat]).sum())
        dup = duplicates(g, point, -1)
        print(f"  {date} {var}: on land {extra}, wet without value {absent}, in equatorial band {band}, seam and fold {dup:.2e}, max {np.nanmax(np.abs(g)):.2f} m/s", flush=True)
        good &= extra == 0 and absent == 0 and band == 0 and dup < tol
        found = src.groups(var, start=date, end=date)[0]
        if not found:
            print(f"    no surface {name} for this date", flush=True)
            continue
        raw = src.snapshot(found[0][1], var, np.float64)[0]
        for lo, hi in ((min_lat, 15), (15, 40), (40, 90)):
            sel = (np.abs(lat) >= lo) & (np.abs(lat) < hi)
            r, rmse, rms = correlation(g[sel], raw[sel])
            print(f"    {lo:g}-{hi:g} deg: corr {r:.3f}, rms difference {rmse:.3f}, rms surface {rms:.3f} m/s", flush=True)
    return good


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--step", required=True, choices=("w", "geostrophic"))
    ap.add_argument("--source", required=True)
    ap.add_argument("--dst", required=True)
    ap.add_argument("--start")
    ap.add_argument("--end")
    ap.add_argument("--min-lat", type=float, default=5.0)
    ap.add_argument("--tol", type=float, default=1e-3, help="largest difference across seam and fold, m/s")
    a = ap.parse_args()
    src = datakit.source(a.source)
    out = datakit.derived(os.path.join(a.dst, a.step), src)
    good = True
    variable = "w" if a.step == "w" else "vozocrtx"
    for date, _ in out.groups(variable, start=a.start, end=a.end)[0]:
        good &= check_w(src, out, date, a.tol) if a.step == "w" else check_geostrophic(src, out, date, a.min_lat, a.tol)
    print("ALL OK" if good else "FAILED")
    sys.exit(0 if good else 1)


if __name__ == "__main__":
    main()
