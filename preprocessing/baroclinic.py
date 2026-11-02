#!/usr/bin/env python3
"""
PP1: split u and v from any source into barotropic and baroclinic parts, in parallel in time.

    ubt = sum_k(u e3u) / sum_k(e3u),    ubc = u - ubt

and likewise for v with e3v. The thicknesses e3u and e3v are those of the mesh, fixed in time.
Writes <dst>/split/<var>/<var>_<tag>_<date>.nc, readable with datakit.derived: each file holds the
baroclinic part under <var> and the barotropic part under NEMO's name for it (vobtcrtx, vobtcrty).
Dates are split round-robin across ranks. Existing files are skipped, so a stopped run can resume,
unless --overwrite is given. A resume stops if an existing file has another precision or layout.

    mpirun -np 16 python baroclinic.py --source glorys12 --start 1993-01 --end 1993-12 --dst /path/to/out
"""
import os

import numpy as np

import _common
import datakit
from datakit import clean
from orcakit import orca

VARS = datakit.BAROTROPIC


def split(field, e3, mask):
    """Barotropic (ny, nx) and baroclinic (nz, ny, nx) parts of a field that may hold NaN on land. Land is NaN."""
    bt, bc = orca.split_barotropic(clean(field, fill=0.0), mask, e3=e3)
    return np.where(mask.any(axis=0), bt, np.nan), np.where(mask, bc, np.nan)


def output_path(dst, var, tag, date):
    return os.path.join(dst, "split", var, f"{var}_{tag}_{date}.nc")


def metrics(src, dtype):
    e3u, e3v, umask, vmask = src.vertical(dtype)
    return {"vozocrtx": (e3u, umask),
            "vomecrty": (e3v, vmask)}


def process(src, date, paths, metric, dst, dtype):
    """Split one date and write its two files."""
    for var, path in zip(VARS, paths):
        e3, mask = metric[var]
        bt, bc = split(src.snapshot(path, var, dtype), e3, mask)
        lat, lon = src.point_coords(datakit.FIELDS[var][0])
        attrs = datakit.provenance("PP1 baroclinic: u - sum_k(u e3)/sum_k(e3), e3 from the mesh", path, source=src.name,
                                   content=f"{var}: baroclinic part; {VARS[var]}: barotropic part")
        datakit.write_field(os.path.join(dst, "split"), var, date, bc, lat, lon, src.time_hours(path), src.tag, src.depth(path), dtype, attrs, {VARS[var]: bt})


def main():
    a = _common.arguments(__doc__).parse_args()
    _common.run(a, list(VARS), lambda tag, date: {output_path(a.dst, v, tag, date): (v, VARS[v]) for v in VARS}, metrics,
                lambda src, date, paths, metric, dtype: process(src, date, paths, metric, a.dst, dtype))


if __name__ == "__main__":
    main()
