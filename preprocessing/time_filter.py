#!/usr/bin/env python3
"""
PP6: filter any fields of any source in time, optionally with a Gaussian in space, in parallel in space.

UNTESTED: written and passing its synthetic tests, but never checked on real data. Check its output
before using it.

Each rank holds the whole record of its tile, with the steps in the order given by --order. The space
filter exchanges a halo between tiles. Each date is then gathered on one rank, so files are written in parallel.
Run on a square number of ranks. Writes <dst>/<name>/<var>/, one file per date, under the same names.

    mpirun -np 16 python time_filter.py --source glorys04 --method deseason --dst /path/to/out
    mpirun -np 16 python time_filter.py --source cmems --method highpass --cutoff 365 --space low --sigma 8 --order time-space --dst /path/to/out
"""
import os
import sys
import time

import numpy as np

import _common
import datakit
from datakit import NAMES, clean
from orcakit import filters, mpi_orca


def step_name(a):
    name = "deseason" if a.method == "deseason" else f"time_{a.method}"
    return f"{name}_{a.space}pass" if a.space else name


def describe(a):
    what = "deseasonalised, trend and %d harmonics of the year" % a.harmonics if a.method == "deseason" else f"Butterworth {a.method}, cutoff {a.cutoff} days, order {a.order_butter}"
    if a.space:
        what += f"; Gaussian {a.space}-pass in space, sigma {a.sigma}, order {a.order}"
    return what


def time_step(series, times, a):
    """Filter (nt, ..., jl, il) in time, level by level to limit the memory. Same precision."""
    dt = float(np.median(np.diff(times)))

    def one(f):
        if a.method == "deseason":
            out = filters.deseason(f, times, a.harmonics)[0]
            return out + f.mean(axis=0) if a.keep_mean else out
        low = filters.time_lowpass(f, dt, a.cutoff, a.order_butter)
        return low if a.method == "lowpass" else f - low

    if series.ndim == 3:
        return one(series).astype(series.dtype)
    for k in range(series.shape[1]):
        series[:, k] = one(series[:, k])
    return series


def space_step(dec, field, mask, a, point, vector):
    low = filters.lowpass_tile(dec, field, mask, a.sigma, point, vector)
    return low if a.space == "low" else field - low


def write_dates(dec, comm, series, src, var, group, n, root, a, dtype):
    """Gather each date on one rank in turn and write it there."""
    batch = max(dec.Px*dec.Py, 1)
    for first in range(0, len(group), batch):
        mine = None
        for r, t in enumerate(range(first, min(first + batch, len(group)))):
            full = mpi_orca.gather(dec, series[t], root=r)
            if dec.serial or dec.rank == r:
                mine = (t, full)
        if mine is None:
            continue
        t, full = mine
        date, path = group[t][0], group[t][n + 1]
        lat, lon = src.point_coords(datakit.FIELDS[var][0])
        attrs = datakit.provenance(f"PP6 time_filter: {describe(a)}", path, source=src.name, method=a.method, order=a.order, sigma=a.sigma or 0.0, cutoff_days=a.cutoff or 0.0)
        datakit.write_field(root, var, date, full, lat, lon, src.time_hours(path), src.tag, src.depth(path), dtype, attrs, depth_name="depthw" if var == "vovecrtz" else "deptht")
    if comm is not None:
        comm.Barrier()


def main():
    ap = _common.arguments(__doc__)
    ap.add_argument("--method", required=True, choices=("deseason", "lowpass", "highpass"))
    ap.add_argument("--harmonics", type=int, default=2)
    ap.add_argument("--keep-mean", action="store_true", help="deseason: keep the time mean")
    ap.add_argument("--cutoff", type=float, help="Butterworth: cutoff period, days")
    ap.add_argument("--order-butter", type=int, default=5)
    ap.add_argument("--space", choices=("low", "high"), help="also filter in space")
    ap.add_argument("--sigma", type=float, help="in grid points of the source")
    ap.add_argument("--order", default="space-time", choices=("space-time", "time-space"))
    ap.add_argument("--vars", nargs="+", default=["u", "v"], help="common or source names")
    a = ap.parse_args()
    if a.method != "deseason" and not a.cutoff:
        ap.error("--cutoff is needed by a Butterworth filter")
    if bool(a.space) != bool(a.sigma):
        ap.error("--space and --sigma go together")
    names = [NAMES.get(v, v) for v in a.vars]
    src = _common.open_source(a)
    ny, nx = src.grid()["lat"].shape
    dec = mpi_orca.Decomp.auto(ny, nx)
    comm, rank = dec.cart_comm, dec.rank
    dtype = np.dtype(a.dtype)
    root = os.path.join(a.dst, step_name(a))
    group, missing = src.groups(*names, start=a.start, end=a.end)
    outputs = [(os.path.join(root, v, f"{v}_{src.tag}_{g[0]}.nc"), v) for g in group for v in names]
    if rank == 0:
        print(f"{src.name}: {len(group)} dates on {dec.Px} x {dec.Py} tiles, to {root}", flush=True)
        if missing:
            print(f"skipped, missing an input: {missing}", flush=True)
    done = all(os.path.exists(p) for p, _ in outputs)
    if not a.overwrite:
        bad = [p for p, v in outputs if os.path.exists(p) and _common.incompatible(p, (v,), dtype)]
        if bad:
            if rank == 0:
                print(f"{len(bad)} existing files are not {a.dtype}, for example {bad[0]}. Use --overwrite or another --dst.", flush=True)
            sys.exit(1)
        if done:
            if rank == 0:
                print("all outputs exist; nothing to do", flush=True)
            return
    if a.dry_run:
        if rank == 0:
            print(f"would write {len(outputs)} files", flush=True)
        return
    jsl, isl = dec.slices()
    box = (jsl.start, jsl.stop, isl.start, isl.stop)
    times = np.array([src.time_hours(g[1])/24 for g in group])
    t0 = time.time()
    for n, var in enumerate(names):
        point, vector = datakit.FIELDS[var]
        mask = src.point_mask(point)[:, jsl, isl]
        series = []
        for g in group:
            field = clean(datakit.read_one(g[n + 1], var, None, box, dtype))
            m = mask if field.ndim == 3 else mask[0]
            series.append(space_step(dec, field, m, a, point, vector) if a.space and a.order == "space-time" else field)
        series = time_step(np.stack(series), times, a)
        if rank == 0:
            print(f"{var}: filtered, {time.time() - t0:.0f} s", flush=True)
        if a.space and a.order == "time-space":
            m = mask if series.ndim == 4 else mask[0]
            for t in range(series.shape[0]):
                series[t] = space_step(dec, series[t], m, a, point, vector)
        write_dates(dec, comm, series, src, var, group, n, root, a, dtype)
        del series
        if rank == 0:
            print(f"{var}: written, {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
