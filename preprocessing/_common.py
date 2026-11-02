"""What the preprocessing steps share: the options, resuming, and the dates split across ranks."""
import argparse
import os
import sys
import time

import netCDF4
import numpy as np

import datakit


def arguments(description):
    ap = argparse.ArgumentParser(description=description, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", required=True, help="glorys04, glorys12, cmems, or the folder of a derived source")
    ap.add_argument("--like", help="for a derived --source: the raw source whose grid it is on")
    ap.add_argument("--dst", required=True)
    ap.add_argument("--start")
    ap.add_argument("--end")
    ap.add_argument("--dtype", default="float32", choices=("float32", "float64"), help="saved precision")
    ap.add_argument("--overwrite", action="store_true", help="recompute existing files, which are otherwise skipped")
    ap.add_argument("--dry-run", action="store_true")
    return ap


def open_source(a):
    """The raw source named by --source, or the derived one in that folder, on the grid of --like."""
    return datakit.derived(a.source, a.like) if a.like else datakit.source(a.source)


def incompatible(path, names, dtype):
    """Whether an existing file lacks one of the variables or has another precision."""
    with netCDF4.Dataset(path) as nc:
        return any(name not in nc.variables or nc[name].dtype != dtype for name in names)


def run(a, inputs, outputs, setup, work, keep=None):
    """
    Process the dates of source a.source, round-robin across ranks.

        inputs  : variables a date needs
        outputs : outputs(tag, date) gives {path: variables in that file}
        setup   : setup(src, dtype) gives what work needs; called only if a rank has work
        work    : work(src, date, paths, context, dtype)
        keep    : keep(groups) gives the dates to write, default all
    """
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank, size = comm.Get_rank(), comm.Get_size()
    src = open_source(a)
    groups, missing = src.groups(*inputs, start=a.start, end=a.end)
    groups = keep(groups) if keep else groups
    if rank == 0:
        print(f"{src.name}: {len(groups)} dates on {size} ranks, to {a.dst}", flush=True)
        if missing:
            print(f"skipped, missing an input: {missing}", flush=True)
    dtype = np.dtype(a.dtype)
    planned = [(g, outputs(src.tag, g[0])) for g in groups[rank::size]]
    bad = [] if a.overwrite else [p for _, out in planned for p, names in out.items() if os.path.exists(p) and incompatible(p, names, dtype)]
    bad = sum(comm.allgather(bad), [])
    if bad:
        if rank == 0:
            print(f"{len(bad)} existing files are not {a.dtype} with all their variables, for example {bad[0]}. Use --overwrite or another --dst.", flush=True)
        sys.exit(1)
    mine = [g for g, out in planned if a.overwrite or not all(os.path.exists(p) for p in out)]
    if a.dry_run:
        print(f"rank {rank} would write {[g[0] for g in mine]}", flush=True)
        return
    t0 = time.time()
    if mine:
        context = setup(src, dtype)
        for n, (date, *paths) in enumerate(mine, 1):
            work(src, date, paths, context, dtype)
            print(f"rank {rank}: {date} ({n}/{len(mine)}, {time.time() - t0:.0f} s)", flush=True)
    comm.Barrier()
    if rank == 0:
        print(f"done in {time.time() - t0:.0f} s", flush=True)
