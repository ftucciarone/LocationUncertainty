#!/usr/bin/env python3
"""
PP4: Gaussian low-pass or high-pass in space of any fields of any source, snapshot by snapshot.

Writes <dst>/lowpass/<var>/ or <dst>/highpass/<var>/, one file per date, under the same names.
Land is NaN. Chain steps by giving a derived folder as --source, with --like for its grid.

    mpirun -np 16 python spatial_filter.py --source glorys04 --keep low --sigma 8 --dst /path/to/out
    mpirun -np 16 python spatial_filter.py --source /path/to/split/split --like glorys04 --keep high --sigma 8 --dst /path/to/out
"""
import os

import _common
import datakit
from datakit import NAMES, clean
from orcakit import filters


def process(src, date, paths, names, a, dtype):
    """Filter the fields of one date and write them."""
    apply = filters.lowpass if a.keep == "low" else filters.highpass
    for var, path in zip(names, paths):
        point, vector = datakit.FIELDS[var]
        field = clean(src.snapshot(path, var, dtype))
        mask = src.point_mask(point)
        out = apply(field, mask if field.ndim == 3 else mask[0], a.sigma, point, vector)
        lat, lon = src.point_coords(point)
        attrs = datakit.provenance(f"PP4 spatial_filter: Gaussian {a.keep}-pass, sigma {a.sigma} grid points", path, source=src.name, sigma=a.sigma)
        datakit.write_field(os.path.join(a.dst, f"{a.keep}pass"), var, date, out, lat, lon, src.time_hours(path), src.tag, src.depth(path), dtype, attrs,
                            depth_name="depthw" if var == "vovecrtz" else "deptht")


def main():
    ap = _common.arguments(__doc__)
    ap.add_argument("--keep", required=True, choices=("low", "high"))
    ap.add_argument("--sigma", required=True, type=float, help="in grid points of the source")
    ap.add_argument("--vars", nargs="+", default=["u", "v"], help="common or source names")
    a = ap.parse_args()
    names = [NAMES.get(v, v) for v in a.vars]
    _common.run(a, names, lambda tag, date: {os.path.join(a.dst, f"{a.keep}pass", v, f"{v}_{tag}_{date}.nc"): (v,) for v in names}, lambda src, dtype: None,
                lambda src, date, paths, context, dtype: process(src, date, paths, names, a, dtype))


if __name__ == "__main__":
    main()
