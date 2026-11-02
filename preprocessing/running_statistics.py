#!/usr/bin/env python3
"""
PP7: variance, standard deviation or mean of any fields of any source over a sliding window of snapshots,
or the covariance of two velocity components (--stat cov), written as uv_tvcov by default. For the covariance,
U and V are first averaged onto T-points, U in x and V in y, with land faces as 0; land T-points are NaN.

Each output is dated by the newest snapshot of its window. Every window is read on its own, so dates
can be split across ranks; each rank holds one window. Writes <dst>/running_<stat>/<var>/, under the same names.

    mpirun -np 16 python running_statistics.py --source glorys04 --window 21 --stride 7 --stat std --dst /path/to/out
"""
import os

import numpy as np

import _common
import datakit
from datakit import NAMES, clean
from noisekit.running import running_cov, running_stat
from orcakit import grid_c, orca


def write(src, path, var, date, out, a, dtype, inputs=None):
    lat, lon = src.point_coords(datakit.FIELDS[var][0])
    what = "cov of u, v averaged onto T-points" if a.stat == "cov" else a.stat
    step = f"PP7 running_statistics: {what} over {a.window} snapshots, stride {a.stride}"
    inputs = inputs or path
    attrs = datakit.provenance(step, inputs, source=src.name)
    datakit.write_field(os.path.join(a.dst, f"running_{a.stat}"), var, date, out, lat, lon, src.time_hours(path), src.tag, src.depth(path), dtype, attrs,
                        depth_name="depthw" if var == "vovecrtz" else "deptht")


def process(src, date, position, names, a, dtype):
    """Compute and write the statistic of the window ending at one date."""
    group = position[date]
    if a.stat == "cov":
        tmask = src.point_mask("T")

        def on_t(g):
            """
            on_t: read one date's u and v and average them onto T-points, so their covariance pairs co-located values.

                Parameters:
                    g (tuple): (date, u path, v path)

                Returns:
                    ut (float): u averaged onto T-points, land faces as 0; NaN on land T-points
                    vt (float): v averaged onto T-points, same conventions
            """
            u, v = (np.nan_to_num(clean(src.snapshot(path, var, dtype))) for path, var in zip(g[1:], names))
            wet = tmask if u.ndim == 3 else tmask[0]
            return (np.where(wet, grid_c.interp_UV2T(u, orca.halo, dir="lon"), np.nan),
                    np.where(wet, grid_c.interp_UV2T(v, orca.halo, dir="lat"), np.nan))

        pairs = (on_t(g) for g in group)
        _, out = next(running_cov(pairs, a.window))
        write(src, group[-1][1], a.cov_name, date, out, a, dtype, {names[0]: group[-1][1],
                                                                   names[1]: group[-1][2]})
        return
    for n, var in enumerate(names):
        snapshots = (clean(src.snapshot(g[n + 1], var, dtype)) for g in group)
        _, out = next(running_stat(snapshots, a.window, stat=a.stat))
        write(src, group[-1][n + 1], var, date, out, a, dtype)


def main():
    ap = _common.arguments(__doc__)
    ap.add_argument("--window", required=True, type=int, help="in snapshots")
    ap.add_argument("--stride", type=int, default=1, help="in snapshots")
    ap.add_argument("--stat", default="var", choices=("var", "std", "mean", "cov"))
    ap.add_argument("--cov-name", default="uv_tvcov", help="variable name of the covariance")
    ap.add_argument("--vars", nargs="+", default=["u", "v"], help="common or source names")
    a = ap.parse_args()
    names = [NAMES.get(v, v) for v in a.vars]
    if a.stat == "cov" and len(names) != 2:
        ap.error("--stat cov needs two --vars")
    written = [a.cov_name] if a.stat == "cov" else names
    position = {}

    def keep(groups):
        for i in range(a.window - 1, len(groups), a.stride):
            position[groups[i][0]] = groups[i - a.window + 1:i + 1]
        return groups[a.window - 1::a.stride]

    _common.run(a, names, lambda tag, date: {os.path.join(a.dst, f"running_{a.stat}", v, f"{v}_{tag}_{date}.nc"): (v,) for v in written}, lambda src, dtype: None,
                lambda src, date, paths, context, dtype: process(src, date, position, names, a, dtype), keep=keep)


if __name__ == "__main__":
    main()
