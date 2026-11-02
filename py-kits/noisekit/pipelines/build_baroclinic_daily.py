#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_baroclinic_daily.py

Standalone MPI batch job: write a per-day, full-depth copy of the GLORYS04
u / v daily snapshots with the BAROTROPIC (transport-weighted depth-mean)
component removed.

    u_bc(x, z, t) = u(x, z, t) - ubar(x, t),
    ubar          = sum_k u * e3u / sum_k e3u          (per water column)

and likewise for v with e3v.  ubar is the vertical-mean flow; u_bc carries
zero net transport by construction.

Output tree (a sibling of the source daily/ tree, identical layout,
dimensions, variable names, coordinates, dtype and _FillValue -- a drop-in
replacement any GLORYS04-reading pipeline can point database_dir at):

    <storage>/data-OPeNDAP-Mercator/glorys04/daily-baroclinic/
        vozocrtx/vozocrtx_glorys04_YYYY-MM-DD.nc
        vomecrty/vomecrty_glorys04_YYYY-MM-DD.nc

Land and below-seafloor cells are NaN (matching the source _FillValue).

Vertical metric: e3u, e3v are reconstructed from the old-format GLORYS2V4
mesh_mask (T-point e3t_0 + e3t_ps + tmask only, no 3-D e3) via
orca.e3uv_from_partial_steps -- the same helper the build_noise_*_stats.py scripts
use.  The barotropic split itself is orca.split_barotropic.

Parallelism (mpi4py): the day list is split round-robin across ranks; each
rank reads, computes and writes its own files independently -- no two ranks
touch the same output file, so plain (serial-per-file) netCDF writes are
enough, no parallel HDF5.  Each output is written to "<name>.tmp<rank>" then
os.replace()d into place, so an interrupted run is safely resumable
(existing outputs are skipped unless --overwrite).

Every rank independently reads the mesh_mask and rebuilds e3u/e3v (~1-2 s,
one 660 MB read served from page cache after the first).  Peak memory per
rank is ~3-4 GB (e3u_3d + e3v_3d + a full-depth working column for u and
v), so keep -np below node_RAM_GB / 4.

Run:
    # one-time, into the project venv (openmpi 4.1.6 mpicc is on PATH):
    MPICC=$(command -v mpicc) \
      pip install mpi4py

    mpirun -np 16 \
      python -m noisekit.pipelines.build_baroclinic_daily

    # options: --overwrite  --years 2015,2016  --limit N  --dst PATH  --dry-run

Full 2 years (731 days x u+v, full 75 levels): ~15 s/day serial
(read 2x ~2.3 s + split + write 2x zlib), so ~3 h on 1 rank, ~12 min on
16 ranks.  ~250-300 GB of output (source dtype is float64).

Prototyped in the GLORYS04 variance notebook ("Remove the barotropic mode"
section); this is that operation over the whole archive, in parallel, with
the result persisted as a dataset instead of held in RAM.
"""

import argparse
import glob
import importlib
import os
import re
import sys
import time
import warnings

import numpy as np
import netCDF4 as netcdf

# netCDF4 + numpy >= 2.5: assigning to a variable triggers an in-place
# ndarray.shape set inside netCDF4; harmless, but noisy over thousands of writes.
warnings.filterwarnings("ignore", category=DeprecationWarning, module="netCDF4")

try:
    from mpi4py import MPI
except ImportError:
    sys.stderr.write(
        "mpi4py is not installed in this interpreter.\n"
        "  MPICC=$(command -v mpicc) "
        "pip install mpi4py\n"
    )
    sys.exit(1)

from orcakit import paths
from orcakit import orca
importlib.reload(orca)

# -- Configuration ----------------------------------------------------------
acc_dtype = np.float32                                          # working precision (output keeps source dtype)

database_dir = os.path.join(paths.get("glorys04"), "daily")
default_dst = os.path.join(paths.get("glorys04"), "daily-baroclinic")
mesh_mask_file = os.path.join(paths.get("glorys04"), "static", "GLORYS2V4_mesh_mask.nc")

VARS = (("vozocrtx", "vozocrtx"), ("vomecrty", "vomecrty"))     # (subdir, varname)
udate_re = re.compile(r"vozocrtx_glorys04_(\d{4})-(\d{2})-(\d{2})\.nc$")
vdate_re = re.compile(r"vomecrty_glorys04_(\d{4})-(\d{2})-(\d{2})\.nc$")

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


def log(msg, all_ranks=False):
    r"""
    log: print a message tagged with the rank.

        Parameters:
            msg        (str): Message
            all_ranks (bool): Print from every rank, not only rank 0

        Returns:
            None

        Formula:
            [\text{rank}/\text{size}]\ \text{msg}

        See also:
            main
    """
    if all_ranks or rank == 0:
        sys.stdout.write(f"[{rank:03d}/{size}] {msg}\n")
        sys.stdout.flush()


# -- CLI ------------------------------------------------------------------
def parse_args():
    r"""
    parse_args: the command-line options of this script.

        Parameters:
            None

        Returns:
            args (Namespace): dst, years, overwrite, limit and dry_run

        Formula:
            \text{dst},\ \text{years},\ \text{overwrite},\ \text{limit},\ \text{dry\_run}

        See also:
            main
    """
    p = argparse.ArgumentParser(description=__doc__.splitlines()[3])
    p.add_argument("--dst", default=default_dst, help="output root (default: %(default)s)")
    p.add_argument("--years", default="all",
                   help="comma list of years to process, e.g. 2015,2016 (default: all on disk)")
    p.add_argument("--overwrite", action="store_true",
                   help="rewrite output files that already exist")
    p.add_argument("--limit", type=int, default=0,
                   help="process only the first N days (test runs)")
    p.add_argument("--dry-run", action="store_true",
                   help="list what each rank would do, write nothing")
    return p.parse_args()


# -- file discovery -----------------------------------------------------
def discover_pairs(years):
    r"""
    discover_pairs: the daily GLORYS04 u and v files with matching dates.

    A day whose v file is missing is skipped with a warning.

        Parameters:
            years (str): 'all', or a comma list of years such as '2015,2016'

        Returns:
            pairs (list): (u path, v path) per day, sorted by date

        Formula:
            \{(u_d, v_d) : d \in \text{dates}(u) \cap \text{dates}(v),\ \text{year}(d) \in \text{years}\}

        See also:
            main
    """
    ufiles = sorted(f for f in glob.glob(os.path.join(database_dir, "vozocrtx",
                                                      "vozocrtx_glorys04_*.nc"))
                    if udate_re.search(os.path.basename(f)))
    if years != "all":
        keep = {int(y) for y in years.split(",")}
        ufiles = [f for f in ufiles
                  if int(udate_re.search(os.path.basename(f)).group(1)) in keep]
    pairs = []
    for uf in ufiles:
        y, m, d = udate_re.search(os.path.basename(uf)).groups()
        vf = os.path.join(database_dir, "vomecrty", f"vomecrty_glorys04_{y}-{m}-{d}.nc")
        if os.path.exists(vf):
            pairs.append((uf, vf))
        else:
            log(f"WARNING: no vomecrty for {os.path.basename(uf)} -- day skipped", all_ranks=True)
    return pairs


# -- vertical metric (per rank) ---------------------------------------
def build_metric():
    r"""
    build_metric: the U and V cell thicknesses and wet masks of the GLORYS04 mesh.

        Parameters:
            None

        Returns:
            e3u   (float): (nz, ny, nx) U-cell thicknesses, from partial steps
            e3v   (float): (nz, ny, nx) V-cell thicknesses, from partial steps
            umask  (bool): (nz, ny, nx) U-point wet mask
            vmask  (bool): (nz, ny, nx) V-point wet mask

        Formula:
            e_{3u},\ e_{3v} = \text{orca.e3uv\_from\_partial\_steps}(e_{3t}^0, e_{3t}^{ps}, t, u, v)

        See also:
            orcakit.orca.e3uv_from_partial_steps
    """
    with netcdf.Dataset(mesh_mask_file) as nc:
        e3t_0 = np.asarray(nc["e3t_0"][0, :])
        e3t_ps = np.asarray(nc["e3t_ps"][0, :, :])
        tmask = np.asarray(nc["tmask"][0, :, :, :]).astype(bool)
        umask = np.asarray(nc["umask"][0, :, :, :]).astype(bool)
        vmask = np.asarray(nc["vmask"][0, :, :, :]).astype(bool)
    e3u, e3v = orca.e3uv_from_partial_steps(e3t_0, e3t_ps, tmask, umask, vmask)
    return e3u.astype(acc_dtype), e3v.astype(acc_dtype), umask, vmask


# -- one snapshot -----------------------------------------------------
def baroclinic_field(src_path, varname, e3, mask3d):
    r"""
    baroclinic_field: the baroclinic part of one daily velocity file.

    The velocity minus its thickness-weighted depth mean. NaN on land and below the sea floor.

        Parameters:
            src_path (str): Daily file
            varname  (str): 'vozocrtx' or 'vomecrty'
            e3     (float): (nz, ny, nx) thicknesses at the field's points
            mask3d  (bool): (nz, ny, nx) wet mask at the field's points

        Returns:
            bc (float): (nz, ny, nx) field minus its e3-weighted depth mean; NaN on land and
                        below the sea floor

        Formula:
            u' = u - \frac{\sum_k e_{3,k} u_k}{\sum_k e_{3,k}}

        See also:
            orcakit.orca.split_barotropic
    """
    with netcdf.Dataset(src_path) as nc:
        full = nc[varname][:, :, :].filled(0.0).astype(acc_dtype)
    _, bc = orca.split_barotropic(full, mask3d, e3=e3)          # (nz, y, x), 0 on land
    return np.where(mask3d, bc, np.nan).astype(acc_dtype)


def clone_with(src_path, dst_path, varname, new3d):
    r"""
    clone_with: copy a file exactly, replacing the data of one variable.

    Dimensions, coordinates, attributes, fill value, compression and chunking are kept. The copy is
    written under a temporary name and renamed, so a stopped run leaves no partial file.

        Parameters:
            src_path (str): File to copy
            dst_path (str): File to write
            varname  (str): Variable whose data are replaced
            new3d  (float): New data

        Returns:
            None

        Formula:
            \text{dst} = \text{src} \text{ with } \text{src}[v] \leftarrow x

        See also:
            main
    """
    tmp = f"{dst_path}.tmp{rank}"
    with netcdf.Dataset(src_path) as s, netcdf.Dataset(tmp, "w") as d:
        for dn, dd in s.dimensions.items():
            d.createDimension(dn, None if dd.isunlimited() else len(dd))
        for a in s.ncattrs():
            d.setncattr(a, s.getncattr(a))
        d.baroclinic_processing = (
            "barotropic (transport-weighted depth-mean) component removed: "
            "u_bc = u - sum_k(u e3u)/sum_k(e3u); e3u/e3v from GLORYS2V4_mesh_mask.nc "
            "via orca.e3uv_from_partial_steps; split via orca.split_barotropic")
        d.baroclinic_source = os.path.basename(src_path)
        d.baroclinic_date = time.strftime("%Y-%m-%dT%H:%M:%S")
        for vn, sv in s.variables.items():
            flt = sv.filters() or {}
            ch = sv.chunking()
            kw = dict(zlib=flt.get("zlib", False),
                      complevel=flt.get("complevel", 4),
                      shuffle=flt.get("shuffle", False))
            if isinstance(ch, (list, tuple)):
                kw["chunksizes"] = ch
            dv = d.createVariable(vn, sv.dtype, sv.dimensions,
                                  fill_value=getattr(sv, "_FillValue", None), **kw)
            for a in sv.ncattrs():
                if a != "_FillValue":
                    dv.setncattr(a, sv.getncattr(a))
            dv[...] = new3d if vn == varname else sv[...]
    os.replace(tmp, dst_path)


# -- main -----------------------------------------------------------------
def main():
    r"""
    main: command-line entry point: write the daily baroclinic velocities of GLORYS04, dates split across ranks.

    An older version of the barotropic and baroclinic split. Its output is wrong in the last column of u
    and the top row of v, where it does not obey the seam and the fold.

        Parameters:
            None; options are read from the command line, see parse_args

        Returns:
            None

        Formula:
            u' = u - \bar u^{z} \text{ for every day, written next to the raw files}

        See also:
            baroclinic_field, orcakit.orca.split_barotropic
    """
    args = parse_args()

    pairs = discover_pairs(args.years) if rank == 0 else None
    pairs = comm.bcast(pairs, root=0)
    if args.limit:
        pairs = pairs[:args.limit]
    n_days = len(pairs)
    if n_days == 0:
        log("no matching days -- nothing to do", all_ranks=True)
        return

    for sub, _ in VARS:
        os.makedirs(os.path.join(args.dst, sub), exist_ok=True)
    comm.Barrier()

    if rank == 0:
        d0 = udate_re.search(os.path.basename(pairs[0][0])).groups()
        d1 = udate_re.search(os.path.basename(pairs[-1][0])).groups()
        log(f"{n_days} paired days {'-'.join(d0)} .. {'-'.join(d1)}, "
            f"{size} rank(s), ~{-(-n_days // size)} days/rank -> {args.dst}"
            + ("   [DRY RUN]" if args.dry_run else ""))

    e3u_3d, e3v_3d, umask, vmask = build_metric()
    metric = {"vozocrtx": (e3u_3d, umask),
              "vomecrty": (e3v_3d, vmask)}
    log("e3u/e3v ready", all_ranks=(rank == 0))

    my = list(range(rank, n_days, size))                        # round-robin
    t0, written, skipped = time.time(), 0, 0
    for n, i in enumerate(my):
        uf, vf = pairs[i]
        outs = {"vozocrtx": os.path.join(args.dst, "vozocrtx", os.path.basename(uf)),
                "vomecrty": os.path.join(args.dst, "vomecrty", os.path.basename(vf))}
        if not args.overwrite and all(os.path.exists(p) for p in outs.values()):
            skipped += 1
            continue
        if args.dry_run:
            log(f"would write {os.path.basename(uf)} (+ vomecrty)", all_ranks=True)
            written += 1
            continue
        for src, (sub, var) in zip((uf, vf), VARS):
            e3, m3 = metric[var]
            clone_with(src, outs[var], var, baroclinic_field(src, var, e3, m3))
        written += 1
        if (n + 1) % 20 == 0 or (n + 1) == len(my):
            log(f"{n + 1}/{len(my)} done (wrote {written}, skipped {skipped}, "
                f"{time.time() - t0:.0f}s)", all_ranks=True)

    comm.Barrier()
    tot_w = comm.reduce(written, op=MPI.SUM, root=0)
    tot_s = comm.reduce(skipped, op=MPI.SUM, root=0)
    if rank == 0:
        log(f"ALL DONE: {tot_w} day(s) written, {tot_s} already present. tree: {args.dst}")


if __name__ == "__main__":
    main()
