#!/usr/bin/env python3
"""
PP3: surface geostrophic velocity from the SSH of any source, in parallel in time.

    Ug = -(g/f) d(ssh)/dy on U-points,    Vg = +(g/f) d(ssh)/dx on V-points

with SSH averaged to F-points over wet T-points only. NaN on land and for |lat| < --min-lat, where
geostrophy does not hold. No ice-cover mask, so velocities under ice are kept.
Writes <dst>/geostrophic/<var>/<var>_<tag>_<date>.nc, readable with datakit.derived.

    mpirun -np 16 python geostrophic.py --source cmems --start 2015-01 --end 2015-12 --dst /path/to/out
"""
import os

import _common
import datakit
from orcakit import grid_c, orca

VARS = ("vozocrtx", "vomecrty")


def output_path(dst, var, tag, date):
    return os.path.join(dst, "geostrophic", var, f"{var}_{tag}_{date}.nc")


def metrics(src, dtype):
    grid = src.grid()
    tmask, umask, vmask, gphiu, gphiv = src.mesh("tmask", "umask", "vmask", "gphiu", "gphiv")
    return dict(tmask=tmask[0].astype(bool), umask=umask[0].astype(bool), vmask=vmask[0].astype(bool), gphiu=gphiu, gphiv=gphiv,
                **{k: grid[k] for k in ("e1t", "e2t", "e1v", "e2u")})


def process(src, date, paths, m, dst, dtype, min_lat):
    """Compute and write Ug and Vg of one date."""
    path = paths[0]
    fields = grid_c.geostrophic_uv(src.snapshot(path, "ssh", dtype), exchange=orca.halo, min_lat=min_lat, **m)
    attrs = datakit.provenance(f"PP3 geostrophic: surface geostrophic velocity from SSH, NaN below {min_lat} deg", path, source=src.name, min_lat=min_lat)
    for var, field in zip(VARS, fields):
        lat, lon = src.point_coords(datakit.FIELDS[var][0])
        datakit.write_field(os.path.join(dst, "geostrophic"), var, date, field, lat, lon, src.time_hours(path), src.tag, dtype=dtype, attrs=attrs)


def main():
    ap = _common.arguments(__doc__)
    ap.add_argument("--min-lat", type=float, default=5.0, help="NaN for |lat| below this, in degrees")
    a = ap.parse_args()
    _common.run(a, ["sossheig"], lambda tag, date: {output_path(a.dst, v, tag, date): (v,) for v in VARS}, metrics,
                lambda src, date, paths, m, dtype: process(src, date, paths, m, a.dst, dtype, a.min_lat))


if __name__ == "__main__":
    main()
