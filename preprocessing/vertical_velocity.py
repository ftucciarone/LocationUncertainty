#!/usr/bin/env python3
"""
PP2: vertical velocity from the horizontal velocity of any source, by continuity, in parallel in time.

    w(k) = w(k+1) - [ di(e2u e3u u) + dj(e1v e3v v) ](k) / (e1t e2t),    zero below the deepest wet level

with the thicknesses of the mesh, fixed in time, as in baroclinic.py. w(k) sits at the top of T-cell k,
so the files use the W levels of the mesh. Uses the full velocity.
Writes <dst>/w/vovecrtz/vovecrtz_<tag>_<date>.nc, readable with datakit.derived. Land is NaN.

    mpirun -np 16 python vertical_velocity.py --source cmems --start 2015-01 --end 2015-12 --dst /path/to/out
"""
import os

import numpy as np

import _common
import datakit
from datakit import clean
from orcakit import grid_c, orca


def output_path(dst, tag, date):
    return os.path.join(dst, "w", "vovecrtz", f"vovecrtz_{tag}_{date}.nc")


def metrics(src, dtype):
    e3u, e3v, umask, vmask = src.vertical(np.float64)
    grid = src.grid()
    return dict(e3u=e3u, e3v=e3v, umask=umask, vmask=vmask, tmask=src.point_mask("T"), **{k: grid[k] for k in ("e1t", "e2t", "e1v", "e2u")})


def process(src, date, paths, m, dst, dtype):
    """Compute and write w of one date."""
    upath, vpath = paths
    w = grid_c.vertical_velocity(clean(src.snapshot(upath, "u", np.float64)), clean(src.snapshot(vpath, "v", np.float64)), exchange=orca.halo, **m)
    lat, lon = src.point_coords("T")
    attrs = datakit.provenance("PP2 vertical_velocity: w from continuity, zero at the bottom, e3u and e3v from the mesh", {"u": upath,
                                                                                                                           "v": vpath}, source=src.name)
    datakit.write_field(os.path.join(dst, "w"), "w", date, w, lat, lon, src.time_hours(upath), src.tag, src.mesh("gdepw_0"), dtype, attrs, depth_name="depthw")


def main():
    a = _common.arguments(__doc__).parse_args()
    _common.run(a, ["vozocrtx", "vomecrty"], lambda tag, date: {output_path(a.dst, tag, date): ("vovecrtz",)}, metrics,
                lambda src, date, paths, m, dtype: process(src, date, paths, m, a.dst, dtype))


if __name__ == "__main__":
    main()
