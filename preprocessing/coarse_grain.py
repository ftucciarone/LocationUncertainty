#!/usr/bin/env python3
"""
PP5: average any fields of any source onto a target grid, snapshot by snapshot, or sample the nearest point.

Each point of the source goes to its nearest wet target cell of the same kind (T, U or V), weighted by
its cell area. Repeated seam and fold points are left out. Levels stay those of the source.
Writes <dst>/coarse/<var>/, under the same names, with the coordinates of the target. No data is NaN.

    mpirun -np 16 python coarse_grain.py --source /path/to/out/lowpass --like glorys04 --dst /path/to/out
"""
import os

import netCDF4
import numpy as np

import _common
import datakit
from datakit import NAMES, clean
from orcakit import orca, paths


def target_grid(path):
    """lat, lon and surface wet mask of the target at T, U and V points, and the shape of its grid as "lat"."""
    with netCDF4.Dataset(path) as nc:
        tmask = (nc["top_level"][0, :, :] > 0).data
        grid = {"lat": nc["nav_lat"][:].data,
                "lon": nc["nav_lon"][:].data}
        for p, mask in (("T", tmask), ("U", tmask & np.roll(tmask, -1, axis=1)), ("V", tmask & np.roll(tmask, -1, axis=0))):
            grid[p] = (np.squeeze(nc[f"gphi{p.lower()}"][:].data), np.squeeze(nc[f"glam{p.lower()}"][:].data), mask)
    return grid


def setup(src, a, names):
    """For each field, the map onto the target and the weights."""
    target = target_grid(a.target)
    grid = src.grid()
    maps = {}
    for point in {datakit.FIELDS[v][0] for v in names}:
        p = point.lower()
        lat, lon = src.mesh(f"gphi{p}", f"glam{p}")
        tlat, tlon, tmask = target[point]
        unique = np.ones(lat.shape, bool)
        unique[-2:], unique[:, -2:] = False, False
        if a.method == "bin":
            maps[point] = dict(index=orca.assign_to_coarse_grid(lat, lon, tlat, tlon, tmask), weight=grid[f"e1{p}"]*grid[f"e2{p}"]*unique)
        else:
            ilat, ilon = orca.find_nearest_grid(tlat, tlon, lat, lon, mask=tmask)
            maps[point] = dict(ilat=ilat, ilon=ilon, mask=tmask)
    return target, maps


def regrid(field, m, shape):
    if "index" in m:
        return orca.coarse_grain(field, m["weight"], m["index"], shape)
    out = orca.downsample_field(field, m["ilat"], m["ilon"], m["mask"])
    return np.where(m["mask"], out, np.nan)


def process(src, date, paths_, names, context, a, dtype):
    """Regrid the fields of one date and write them."""
    target, maps = context
    for var, path in zip(names, paths_):
        field = clean(src.snapshot(path, var, dtype))
        out = regrid(field, maps[datakit.FIELDS[var][0]], target["lat"].shape)
        attrs = datakit.provenance(f"PP5 coarse_grain: {a.method} onto {os.path.basename(a.target)}", path, source=src.name, target=os.path.basename(a.target))
        tlat, tlon, _ = target[datakit.FIELDS[var][0]]
        datakit.write_field(os.path.join(a.dst, "coarse"), var, date, out, tlat, tlon, src.time_hours(path), src.tag, src.depth(path), dtype, attrs,
                            depth_name="depthw" if var == "vovecrtz" else "deptht")


def main():
    ap = _common.arguments(__doc__)
    ap.add_argument("--target", default=paths.get("domain_orca1"), help="domain_cfg of the target grid")
    ap.add_argument("--method", default="bin", choices=("bin", "nearest"))
    ap.add_argument("--vars", nargs="+", default=["u", "v"], help="common or source names")
    a = ap.parse_args()
    names = [NAMES.get(v, v) for v in a.vars]
    _common.run(a, names, lambda tag, date: {os.path.join(a.dst, "coarse", v, f"{v}_{tag}_{date}.nc"): (v,) for v in names}, lambda src, dtype: setup(src, a, names),
                lambda src, date, paths_, context, dtype: process(src, date, paths_, names, context, a, dtype))


if __name__ == "__main__":
    main()
