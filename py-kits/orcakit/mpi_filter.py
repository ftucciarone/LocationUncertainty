#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mpi_filter.py
=============

An MPI-parallel spatial filter over a global-ocean 3-D field, using a 2-D
horizontal domain decomposition with FULL VERTICAL PENCILS -- every MPI rank
owns *all* nz vertical levels but only a rectangular sub-block of the horizontal
plane.

Why this shape of decomposition
-------------------------------
The motivating operation is a per-level horizontal Gaussian low-pass. Such a
filter:
  * couples neighbouring points horizontally (it IS a horizontal stencil), so a
    horizontal split needs halo exchange, but
  * does not couple vertical levels (a 2-D Gaussian, sigma_z = 0, is applied
    level by level).
Keeping each rank's vertical column whole ("pencil") therefore means the filter
never has to communicate in z; all the MPI traffic is a thin horizontal halo.

Two run modes
-------------
  * NO MPI  (`python mpi_filter.py ...`)  or  `mpirun -np 1 ...`
        one rank reads the whole (nz, ny, nx) field and filters it. This is not
        a separate code path -- it is the P = 1 corner of the parallel code
        (see `Tile.from_whole` / the `cart_comm is None` branches).

  * `mpirun -np N ...`,  N a perfect square
        P = sqrt(N) tiles per horizontal axis, so a P x P grid of ranks. Each
        rank reads only its own hyperslab  var[:, j0:j1, i0:i1]  (full depth),
        surrounds it with a halo of width  hw = ceil(4*sigma)  filled by
        point-to-point exchange with its Cartesian neighbours, filters locally,
        and (optionally) the tiles are gathered back to rank 0.

Process grid (example: N = 9, P = 3)
------------------------------------
          cx=0       cx=1       cx=2
        +----------+----------+----------+
   cy=2 |          |          |          |   <- top tile-row  = global NORTH edge
        |          |          |          |      (tripolar fold)
        +----------+----------+----------+
   cy=1 |          |          |          |   (interior tile-row: ordinary halo
        |          |          |          |    exchange on all four sides)
        +----------+----------+----------+
   cy=0 |          |          |          |   <- bottom tile-row = global SOUTH
        |          |          |          |      edge (closed)
        +----------+----------+----------+
        x is PERIODIC: the cx=0 and cx=2 columns are neighbours (E-W wrap).
        y is NOT periodic: the cy=0 / cy=2 rows have a real domain edge.
   (cy, cx) = cart_comm.Get_coords(rank).  With reorder=True the integer rank
   labels may be permuted by MPI, so we always work from the coords, never from
   a hand-computed rank number.

Per-rank padded array (`Tile.a`, shape (nz, jl + 2*hw, il + 2*hw))
-----------------------------------------------------------------
                       x  ------------------------------------->
            col: 0        hw                     hw+il        hw+il+hw
               +----------+----------------------+----------+
        row 0  |  NW      |      NORTH halo      |     NE   |
               |  corner  |    a[:, 0:hw,        |   corner |
               |          |        hw:hw+il]     |          |
            hw +----------+----------------------+----------+
               |  WEST    |                      |   EAST   |
               |  halo    |      INTERIOR        |   halo   |
               | a[:, :,  |   a[:, hw:hw+jl,     | a[:, :,  |
               |  0:hw]   |       hw:hw+il]      | hw+il: ] |
         hw+jl +----------+----------------------+----------+
               |  SW      |      SOUTH halo      |     SE   |
               |  corner  |    a[:, hw+jl:,      |   corner |
               |          |        hw:hw+il]     |          |
       hw+jl+hw+----------+----------------------+----------+
  INTERIOR  = this rank's own data (jl rows x il cols, all nz levels).
  halos     = width hw, filled by `_exchange` every time we filter.
  hw        = ceil(4*sigma): a Gaussian truncated at 4 sigma reaches exactly
              this far, so one halo of this width is all the smoother needs.

Halo / boundary rules (implemented in `Tile._exchange`)
------------------------------------------------------
  x  : PERIODIC, period nx-2: the last two columns repeat the first two.
       Cartesian dim 1 is periodic => every rank always has a valid east and west
       neighbour (they wrap); the outermost tiles shift what they send by those
       two columns.
  y  : NON-periodic.
         interior tile faces  : ordinary neighbour exchange.
         global SOUTH (cy = 0) : closed. South halo = interior rows reflected,
                                 the same rule scipy applies in `orca.smooth_2D`.
         global NORTH (cy=Py-1): T-pivot fold, as `orca.halo`: the ghosts come
                                 from the source row and column of `orca.FOLD`,
                                 which depend on the grid point, times fold_sign
                                 (-1 for either velocity component, +1 otherwise
                                 and for the wet mask). Because each rank holds
                                 only part of that global row, the top tile-row
                                 pools its strips with an Allgatherv on a row
                                 sub-communicator. The corners are filled too.

Expected input layout
---------------------
  * one netCDF file per (component, date); the netCDF variable name equals the
    component name passed on the command line;
  * the data variable is 3-D (level, y, x), e.g. (75, 1021, 1442), with no time
    axis;
  * _FillValue is NaN  ->  masked_invalid  ->  data filled with 0 + a separate
    0/1 "wet" mask;
  * `--meridional` marks a field whose sign flips across the tripolar fold
    (either velocity component); by default the fold sign is +1. `--point` says
    where the field lives (T, U, V or F).

Sketch-level simplifications -- read before trusting the numbers
-------------------------------------------------------------
  1. hw + 5 must be <= the smallest local tile dimension: the halo is filled from
     the IMMEDIATE neighbour only (no multi-hop forwarding), and the fold reaches
     5 rows further down. This is asserted collectively in `Tile.from_mpi`; if it
     trips, use fewer ranks or a smaller sigma.
  2. Each rank opens the netCDF independently (no parallel HDF5). That is fine
     for read-only access; the (optional) write goes through a gather to rank 0.

Every function docstring below carries a `Called:` line saying whether it is
user-facing (you call it from a script / the CLI) or automatic (another routine
here calls it), and by whom.

Usage
-----
    # serial, whole field, cross-check against the reference smoother
    python mpi_filter.py --component vozocrtx --date 2015-01-01 --sigma 8 --check

    # 16-way (4x4) tiling, gather the result on rank 0 and cross-check
    mpirun -np 16 python mpi_filter.py --sigma 8 --check

The gathered result is tiling-invariant and matches the whole-domain reference.
"""

import argparse
import os
import sys
import time

import numpy as np

from orcakit import paths
import netCDF4 as netcdf

# --------------------------------------------------------------------------
# A reference whole-domain smoother (`orca.smooth_2D`, with its own global
# halo + tripolar fold) is used ONLY by `--check`; it lives outside this file,
# so its directory is put on sys.path here. The parallel filter itself calls
# scipy directly and does not need it.
# --------------------------------------------------------------------------
from orcakit import orca

# --------------------------------------------------------------------------
# mpi4py is an OPTIONAL import. If it is missing, or if the job is launched
# with a single process, the script silently runs the whole-domain serial
# path. `_HAVE_MPI` records which world we are in; `MPI` is only dereferenced
# when `_HAVE_MPI` is True.
# --------------------------------------------------------------------------
try:
    from mpi4py import MPI
    _HAVE_MPI = True
except Exception:                       # ImportError, or a broken MPI runtime
    _HAVE_MPI = False

# -- Configuration --------------------------------------------------------
# `database_dir` is expected to contain one sub-directory per
# component, each holding files named  <component>_<tag>_<date>.nc .
acc_dtype = np.float32                  # working precision for everything on-node
database_dir = os.path.join(paths.get("glorys04"), "daily")
file_tag = "glorys04"                   # middle token of the file name


# ==========================================================================
#  Small helpers
# ==========================================================================

def default_path(component, date):
    r"""
    default_path: path of the daily file of one variable and date, in this script's data folder.

    The folder and file tag are set at the top of this script.

        Parameters:
            component (str): Variable name, also the folder and file prefix, e.g. 'vozocrtx'
            date      (str): Date in the file name, e.g. '2015-01-01'

        Returns:
            path (str): The path; no file is opened

        Formula:
            \text{path} = \text{database\_dir}/c/c\_\text{tag}\_\text{date}.nc

        See also:
            main
    """
    return os.path.join(database_dir, component, f"{component}_{file_tag}_{date}.nc")


def _split(n, p):
    r"""
    _split: cut range(n) into p contiguous chunks of nearly equal size.

    Places the tiles of a decomposition in the global grid. The first (n mod p) chunks get one extra cell,
    as numpy.array_split does: _split(1442, 3) -> [(0, 481), (481, 962), (962, 1442)].

        Parameters:
            n (int): Global size along one axis
            p (int): Number of chunks, the tiles along that axis

        Returns:
            out (list): p (start, stop) tuples

        Formula:
            n = q\,p + r,\qquad \text{chunk } k \text{ has } q + [k < r] \text{ cells},\qquad \text{start}_k = k\,q + \min(k, r)

        See also:
            mpi_orca._split
    """
    base, rem = divmod(n, p)
    out, s = [], 0
    for k in range(p):
        e = s + base + (1 if k < rem else 0)   # +1 for the first `rem` chunks
        out.append((s, e))
        s = e
    return out


def _read_slab(path, varname, jsl=slice(None), isl=slice(None)):
    r"""
    _read_slab: read all levels of a horizontal window of a variable from disk.

    The field and its wet mask are returned apart because the land-aware smoother averages over wet points
    only, smoothing the field times the mask and the mask with the same kernel.

        Parameters:
            path    (str): NetCDF file
            varname (str): 3D (level, y, x) variable
            jsl   (slice): Rows to read, all by default
            isl   (slice): Columns to read, all by default

        Returns:
            data (float): (nz, ny_sel, nx_sel) float32 field, NaN and fill as 0
            wet  (float): Same shape, 1 where the value read was finite, else 0

        Formula:
            d = f[:, \text{jsl}, \text{isl}] \text{ with NaN} \to 0,\qquad w = [f \text{ finite}]

        See also:
            mpi_orca.read_slab
    """
    with netcdf.Dataset(path) as nc:
        # netCDF4 returns a masked array; _FillValue here is NaN, and
        # masked_invalid also catches any stray Inf. `[:, jsl, isl]` reads only
        # the requested hyperslab off disk -- each rank touches just its tile.
        raw = np.ma.masked_invalid(nc[varname][:, jsl, isl])
    filled = np.ma.filled(raw, 0.0)
    data = np.ascontiguousarray(filled, dtype=acc_dtype)
    # wet = "was finite in the source". `getmaskarray` gives the mask even when
    # `raw.mask` is the scalar False; the `& isfinite` is belt-and-braces in
    # case a build of netCDF4 hands back a plain ndarray with raw NaNs.
    wet = np.ascontiguousarray(
        ~np.ma.getmaskarray(raw) & np.isfinite(filled), dtype=acc_dtype
    )
    return data, wet


# ==========================================================================
#  Tile: one rank's full-depth horizontal block, its halo, and the filter
# ==========================================================================

class Tile:
    """One process's data + halo + local filter.

    Construction is via the two classmethods, not __init__ directly:
        Tile.from_whole(...)  -> serial, one tile == the whole globe (P = 1)
        Tile.from_mpi(...)    -> parallel, one tile == a P x P sub-block

    Everything after construction (halo exchange, smoothing, gather) is written
    so that `cart_comm is None` (serial) is just the P = 1 special case -- that
    is what makes the `--check` comparison between the two meaningful.

    Index cheat-sheet for the padded array `self.a` (see the module docstring
    diagram), with hw = self.hw, jl = self.jl (interior rows), il = self.il
    (interior cols):
        interior     a[:, hw:hw+jl, hw:hw+il]
        west  halo   a[:, :,        0:hw]
        east  halo   a[:, :,        hw+il:]
        south halo   a[:, 0:hw,     :]
        north halo   a[:, hw+jl:,   :]
    """

    def __init__(self, data, wet, hw, is_meridional, point="T",
                 cart_comm=None, Px=1, Py=1, cx=0, cy=0, i0=0, j0=0, ny_g=None, nx_g=None):
        r"""
        __init__: low-level constructor; use Tile.from_whole or Tile.from_mpi instead.

        Allocates the padded buffers of the field and its mask and fills their interior.

            Parameters:
                data, wet  (float): (nz, jl, il) this tile's field and wet mask
                hw           (int): Halo width, ceil(4 sigma)
                is_meridional (bool): True if the field's sign flips across the fold
                point        (str): 'T', 'U', 'V' or 'F', where the field lives
                cart_comm   (Comm): Cartesian communicator, or None for serial
                Px, Py       (int): Tiles per axis, 1 when serial
                cx, cy       (int): This tile's column and row in the P x P grid
                i0, j0       (int): Global index of its first interior column and row
                ny_g, nx_g   (int): Global sizes, by default those of data (serial)

            Returns:
                None

            Formula:
                a[:, h:h+j_l,\ h:h+i_l] = d,\qquad w[:, h:h+j_l,\ h:h+i_l] = m,\qquad h = \lceil 4\sigma \rceil

            See also:
                Tile.from_whole, Tile.from_mpi
        """
        # ---- geometry of THIS tile ----
        self.nz, self.jl, self.il = data.shape       # interior sizes (jl rows, il cols)
        self.hw = int(hw)                             # halo width (each side)
        self.is_meridional = bool(is_meridional)      # True => sign flip across the fold
        self.point = point                            # T, U, V or F: picks the fold mapping
        # ---- geometry of the GLOBAL grid + where this tile sits in it ----
        self.cart_comm = cart_comm                   # Cartesian communicator, or None (serial)
        self.Px, self.Py = Px, Py                    # tiles per axis (x, y)
        self.cx, self.cy = cx, cy                    # this tile's (col, row) index in the P x P grid
        self.i0, self.j0 = i0, j0                    # global index of this tile's first interior col / row
        self.ny_g = data.shape[1] if ny_g is None else ny_g   # global ny (== jl when serial)
        self.nx_g = data.shape[2] if nx_g is None else nx_g   # global nx (== il when serial)
        self.rank = cart_comm.Get_rank() if cart_comm is not None else 0

        # Row sub-communicator: all ranks sharing this tile's `cy` (same
        # tile-row), ordered by `cx`. Only the top tile-row actually uses it,
        # to Allgather the strips it needs for the tripolar fold. `Sub([keep_y,
        # keep_x])` -> keep the x dimension, collapse y.
        self.row_comm = cart_comm.Sub([False, True]) if cart_comm is not None else None

        # ---- allocate the padded buffers and drop the interior in ----
        hw = self.hw
        self.wet = wet                               # UNpadded interior wet mask (nz, jl, il)
        pad = (self.nz, self.jl + 2 * hw, self.il + 2 * hw)
        self.a = np.zeros(pad, dtype=acc_dtype)      # padded field   (halo starts at 0)
        self.w = np.zeros(pad, dtype=acc_dtype)      # padded wet mask (halo starts at 0)
        self.a[:, hw:hw + self.jl, hw:hw + self.il] = data
        self.w[:, hw:hw + self.jl, hw:hw + self.il] = wet

    # ---------------------------------------------------------------- constructors
    @classmethod
    def from_whole(cls, path, varname, sigma, is_meridional, point="T"):
        r"""
        from_whole: a single tile holding the whole field, with its halo filled locally.

        The serial case, also the reference that a tiled run should reproduce.

            Parameters:
                path          (str): NetCDF file
                varname       (str): Variable to load
                sigma       (float): Gaussian lengthscale the tile will be smoothed with, in grid points
                is_meridional (bool): True if the field's sign flips across the fold
                point         (str): 'T', 'U', 'V' or 'F', where the field lives

            Returns:
                tile (Tile): Tile spanning the global grid, cart_comm None

            Formula:
                \text{one tile: } (j_0, i_0) = (0, 0),\ (j_l, i_l) = (n_y, n_x),\ h = \lceil 4\sigma \rceil

            See also:
                Tile.from_mpi
        """
        data, wet = _read_slab(path, varname)                 # whole (nz, ny, nx)
        hw = int(np.ceil(4.0 * sigma))                        # Gaussian reach at 4 sigma
        return cls(data, wet, hw, is_meridional, point, cart_comm=None)   # cart_comm=None => serial branch everywhere

    @classmethod
    def from_mpi(cls, comm, path, varname, sigma, is_meridional, point="T"):
        r"""
        from_mpi: this rank's tile of a P x P decomposition, read as one window of the file.

        Each rank reads only its own part of the field, so a grid too large for one process can be filtered.

        Exits the program if the number of ranks is not a perfect square, or if the halo would be wider than
        the thinnest tile.

            Parameters:
                comm         (Comm): MPI communicator of size P*P
                path          (str): NetCDF file
                varname       (str): Variable to load
                sigma       (float): Gaussian lengthscale, in grid points; sets hw = ceil(4 sigma)
                is_meridional (bool): True if the field's sign flips across the fold
                point         (str): 'T', 'U', 'V' or 'F', where the field lives

            Returns:
                tile (Tile): This rank's tile

            Formula:
                P = \sqrt{N},\ h = \lceil 4\sigma \rceil \le \min(j_l, i_l),\ \text{tile from \_split}

            See also:
                Tile.from_whole, mpi_orca.Decomp.from_comm
        """
        N = comm.Get_size()
        P = int(round(N ** 0.5))
        if P * P != N:
            # "sqrt(N) tiles per axis" only makes sense for a square N.
            if comm.Get_rank() == 0:
                sys.stderr.write(
                    f"mpi_filter: need a perfect-square process count, got {N}; "
                    f"try one of {[k * k for k in range(2, 12)]}\n")
            comm.Barrier()          # let rank 0's message flush before everyone exits
            sys.exit(1)

        # Cartesian topology over a [Py, Px] = [P, P] grid.
        #   dims    = [P, P]            axis 0 = y/lat, axis 1 = x/lon
        #   periods = [False, True]     y has real edges; x wraps (E-W)
        #   reorder = True             MPI may renumber ranks to fit the hardware
        cart_comm = comm.Create_cart([P, P], periods=[False, True], reorder=True)
        # Get_coords returns (coord_axis0, coord_axis1) = (cy, cx) for THIS rank.
        cy, cx = cart_comm.Get_coords(cart_comm.Get_rank())

        # Read the global shape (cheap: header only), then this tile's window.
        with netcdf.Dataset(path) as nc:
            nz, ny_g, nx_g = nc[varname].shape
        j0, j1 = _split(ny_g, P)[cy]        # this tile's global row range [j0, j1)
        i0, i1 = _split(nx_g, P)[cx]        # this tile's global col range [i0, i1)
        data, wet = _read_slab(path, varname, slice(j0, j1), slice(i0, i1))

        # The halo is filled from the IMMEDIATE neighbour only, and the fold
        # reaches 5 rows further down; if that is wider than the thinnest tile,
        # a single exchange cannot fill it. Check it collectively (MIN over all ranks) so every rank aborts together.
        hw = int(np.ceil(4.0 * sigma))
        mn = cart_comm.allreduce(min(j1 - j0, i1 - i0), op=MPI.MIN)
        if hw + 5 > mn:
            if cart_comm.Get_rank() == 0:
                sys.stderr.write(
                    f"mpi_filter: halo hw+5={hw + 5} exceeds the smallest local tile dim {mn} "
                    f"(P={P}); use fewer ranks or a smaller sigma.\n")
            cart_comm.Barrier()
            sys.exit(1)

        return cls(data, wet, hw, is_meridional, point, cart_comm=cart_comm, Px=P, Py=P, cx=cx, cy=cy,
                   i0=i0, j0=j0, ny_g=ny_g, nx_g=nx_g)

    # ---------------------------------------------------------------- halo exchange
    def _exchange(self, a, fold_sign):
        r"""
        _exchange: fill the halo of the tile's field or mask, in place.

        The field and its mask are exchanged with the same pattern, so their halos agree when one is divided
        by the other; they differ only in the sign at the fold. In serial one tile spans the globe; in
        parallel each edge comes from the neighbouring tile.

            Parameters:
                a        (float): self.a or self.w, (nz, jl+2hw, il+2hw), interior at a[:, hw:hw+jl, hw:hw+il]
                fold_sign (float): +1 for a field that does not flip and for the mask, -1 for one that does

            Returns:
                None; a's halo is written in place

            Formula:
                \text{as mpi\_orca.exchange\_halo, with sign } s \text{ at the fold}

            NEMO:
                src/OCE/LBC/lbc_lnk_pt2pt_generic.h90: SUBROUTINE lbc_lnk_pt2pt_/**/PRECISION( cdname, ptab, cd_nat, psgn, kfld, kfillmode, pfillval, lsend, lrecv, ld4only, ldfull )

            See also:
                mpi_orca.exchange_halo
        """
        hw, jl, il = self.hw, self.jl, self.il

        # ---------------- SERIAL (P = 1): one tile already spans the globe ----
        if self.cart_comm is None:
            # x periodic with period nx-2: copy the opposite interior edge,
            # skipping the two repeated columns, into each x-halo.
            a[:, :, 0:hw]        = a[:, :, il - 2:il + hw - 2]
            a[:, :, hw + il:]    = a[:, :, hw + 2:2 * hw + 2]
            # south closed: reflect the interior rows into the south halo.
            a[:, 0:hw, :]        = a[:, 2 * hw - 1:hw - 1:-1, :]
            # north: tripolar fold (the single tile is the whole global row).
            self._north_fold(a, fold_sign, my_i0=0)
            return

        # ---------------- PARALLEL: neighbour exchange ------------------------
        c = self.cart_comm

        # ===== x direction (Cartesian axis 1), PERIODIC =====================
        # cart_comm.Shift(direction, disp) -> (rank_source, rank_dest) where
        #   rank_dest   = neighbour at coord[dir] + disp   (I SEND toward here)
        #   rank_source = neighbour at coord[dir] - disp   (I RECEIVE from here)
        # x is periodic, so both are always real ranks (they wrap around). The
        # last two global columns repeat the first two, so across the seam the
        # outermost tiles send the columns just inside them.
        shift_e = 2 if self.cx == self.Px - 1 else 0
        shift_w = 2 if self.cx == 0 else 0

        # (a) send my EAST interior strip; receive into my WEST halo.
        src_w, dst_e = c.Shift(1, +1)                       # dst_e = east nbr, src_w = west nbr
        send_e = np.ascontiguousarray(a[:, :, il - shift_e:il + hw - shift_e])  # my last hw interior columns, full padded height
        recv_w = np.empty_like(send_e)
        c.Sendrecv(send_e, dest=dst_e, sendtag=0,
                   recvbuf=recv_w, source=src_w, recvtag=0)
        a[:, :, 0:hw] = recv_w                              # west neighbour's east strip -> my west halo

        # (b) send my WEST interior strip; receive into my EAST halo.
        src_e, dst_w = c.Shift(1, -1)                       # dst_w = west nbr, src_e = east nbr
        send_w = np.ascontiguousarray(a[:, :, hw + shift_w:2 * hw + shift_w])   # my first hw interior columns
        recv_e = np.empty_like(send_w)
        c.Sendrecv(send_w, dest=dst_w, sendtag=1,
                   recvbuf=recv_e, source=src_e, recvtag=1)
        a[:, :, hw + il:] = recv_e
        # distinct tags (0 vs 1) so the two exchanges to the SAME partner
        # (which happens at P = 2, where east nbr == west nbr) cannot cross.
        # The x strips span the FULL padded height on purpose: they carry the
        # rows that will become the N/S halo, so the corners get filled when
        # the y exchange (below) copies whole padded-width strips.

        # ===== y direction (Cartesian axis 0), NON-periodic ================
        # (c) send my NORTH interior strip; receive into my SOUTH halo.
        src_s, dst_n = c.Shift(0, +1)                       # dst_n = north nbr, src_s = south nbr
        send_n = np.ascontiguousarray(a[:, hw + jl - hw:hw + jl, :])   # my last hw interior rows, full padded width
        recv_s = np.empty_like(send_n)
        c.Sendrecv(send_n, dest=dst_n, sendtag=2,
                   recvbuf=recv_s, source=src_s, recvtag=2)
        if src_s != MPI.PROC_NULL:          # PROC_NULL => I am the global south edge; recv_s is garbage
            a[:, 0:hw, :] = recv_s

        # (d) send my SOUTH interior strip; receive into my NORTH halo.
        src_n, dst_s = c.Shift(0, -1)                       # dst_s = south nbr, src_n = north nbr
        send_s = np.ascontiguousarray(a[:, hw:2 * hw, :])   # my first hw interior rows
        recv_n = np.empty_like(send_s)
        c.Sendrecv(send_s, dest=dst_s, sendtag=3,
                   recvbuf=recv_n, source=src_n, recvtag=3)
        if src_n != MPI.PROC_NULL:          # PROC_NULL => I am the global north edge; fold fills it instead
            a[:, hw + jl:, :] = recv_n

        # ===== global y boundaries (only the outermost tile-rows) ==========
        if self.cy == 0:                    # bottom tile-row: closed south
            a[:, 0:hw, :] = a[:, 2 * hw - 1:hw - 1:-1, :]  # reflect the interior rows
        if self.cy == self.Py - 1:          # top tile-row: tripolar fold
            self._north_fold(a, fold_sign, my_i0=self.i0)

    def _north_fold(self, a, fold_sign, my_i0):
        r"""
        _north_fold: fill the north halo of a padded array from the T-pivot fold, in place.

        Beyond the top row the ORCA grid folds back on itself: each ghost row is a real row further south,
        read in reverse, with the sign flipped for a velocity component. In parallel the tiles of the top row
        pool their strips to find the mirrored columns.

            Parameters:
                a        (float): Padded array being filled
                fold_sign (float): -1 for a field that flips sign at the fold, else +1
                my_i0      (int): This tile's first global column, 0 when serial

            Returns:
                None; a's north halo rows are written in place

            Formula:
                a_{n_y+p,\,i} = s\;a_{n_y-r-p,\ \mathrm{fold}(i)},\qquad p = 0, \dots, h-1

            NEMO:
                src/OCE/LBC/lbc_nfd_generic.h90: ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)

            See also:
                mpi_orca._north_fold
        """
        hw, jl, il, nz = self.hw, self.jl, self.il, self.nz
        row = orca.FOLD[self.point][0]
        depth = hw + row                                     # block row s <-> global row ny_g-depth+s

        # This rank's contribution: its top `depth` interior rows, all columns.
        contrib = np.ascontiguousarray(
            a[:, hw + jl - depth:hw + jl, hw:hw + il], dtype=acc_dtype
        )  # shape (nz, depth, il)

        if self.cart_comm is None:
            # Serial: the single tile already spans the globe in x.
            strip = contrib                                   # (nz, depth, nx_g)
        else:
            # Parallel: gather every top-row tile's block along x.
            # Allgather"v" because tiles can have different `il` (non-divisible
            # nx_g). counts[k] / displ[k] describe rank k's slot in the flat
            # receive buffer; row_comm orders ranks by cx, so the blocks come
            # back already west-to-east.
            counts = np.array(self.row_comm.allgather(contrib.size), dtype=int)
            widths = np.array(self.row_comm.allgather(il), dtype=int)
            displ = np.concatenate(([0], np.cumsum(counts)[:-1]))
            recv = np.empty(int(counts.sum()), dtype=acc_dtype)
            self.row_comm.Allgatherv([contrib.ravel(), MPI.FLOAT],
                                     [recv, counts, displ, MPI.FLOAT])
            parts = [recv[d:d + n].reshape(nz, depth, w)
                     for n, d, w in zip(counts, displ, widths)]
            strip = np.concatenate(parts, axis=2)            # (nz, depth, nx_g)

        source = orca.fold_columns(np.arange(my_i0 - hw, my_i0 + il + hw), self.nx_g, self.point)
        for p in range(hw):
            a[:, hw + jl + p, :] = fold_sign * strip[:, depth - row - p, source]

    # ---------------------------------------------------------------- the filter
    def smooth(self, sigma):
        r"""
        smooth: land-aware Gaussian low-pass of this tile.

        The field times the mask and the mask are smoothed with the same kernel and divided, so land does not
        pull the result towards zero near the coasts.

        Where the smoothed mask is below 0.01, almost all land, the result is set to 0, as on land.

            Parameters:
                sigma (float): Gaussian lengthscale, in grid points; levels are not mixed

            Returns:
                field_out (float): (nz, jl, il) smoothed interior, exactly 0 on land

            Formula:
                \tilde f = \frac{G_\sigma * (f\,m)}{G_\sigma * m}\;m \quad\text{where } G_\sigma * m > 0.01,\ \text{else } 0

            See also:
                filters.lowpass_tile, filters.masked_smooth
        """
        from scipy.ndimage import gaussian_filter

        # Fill halos first. The field's north-fold sign depends on the variable;
        # the mask's never does.
        fold_sign = -1.0 if self.is_meridional else 1.0
        self._exchange(self.a, fold_sign)     # padded field
        self._exchange(self.w, 1.0)           # padded wet mask

        # sigma per axis: 0 in z (levels stay independent), `sigma` in y and x.
        sig = (0.0, sigma, sigma)
        # `mode="nearest"` only affects the outermost ring of the PADDED array,
        # which the kernel (radius 4*sigma = hw) barely reaches from any
        # interior point -- so the choice is practically irrelevant here.
        num = gaussian_filter(self.a * self.w, sigma=sig, mode="nearest")
        den = gaussian_filter(self.w,          sigma=sig, mode="nearest")
        with np.errstate(invalid="ignore", divide="ignore"):
            sm = np.where(den > 0.01, num / den, 0.0)

        # Strip the halo, re-apply the interior wet mask (land -> exactly 0).
        hw = self.hw
        inner = sm[:, hw:hw + self.jl, hw:hw + self.il]
        return np.ascontiguousarray(inner * self.wet, dtype=acc_dtype)

    # ---------------------------------------------------------------- gather
    def gather_to_root(self, local):
        r"""
        gather_to_root: assemble the tiles of a field into one global array on rank 0.

        Needed to write or check the whole field. It can be a few hundred MB on rank 0, so ask for it only
        when the whole field is needed.

            Parameters:
                local (float): This rank's (nz, jl, il) result, or (nz, ny, nx) when serial

            Returns:
                field (float): (nz, ny_g, nx_g) on rank 0, None on the other ranks; local itself when serial

            Formula:
                F[:, j_0^{(r)} : j_0^{(r)} + j_l^{(r)},\ i_0^{(r)} : i_0^{(r)} + i_l^{(r)}] = f^{(r)} \quad \text{on rank 0}

            See also:
                mpi_orca.gather
        """
        if self.cart_comm is None:
            return local                                    # serial: already global

        c = self.cart_comm
        # Every rank tells rank 0 which global slice it owns.
        meta = c.gather((self.j0, self.j0 + self.jl, self.i0, self.i0 + self.il), root=0)

        if self.rank != 0:
            # Uppercase Send = fast buffer transfer (no pickle). One tag is
            # enough: rank 0 receives in a fixed order.
            c.Send([np.ascontiguousarray(local), MPI.FLOAT], dest=0, tag=7)
            return None

        out = np.zeros((self.nz, self.ny_g, self.nx_g), dtype=acc_dtype)
        out[:, self.j0:self.j0 + self.jl, self.i0:self.i0 + self.il] = local   # rank 0's own tile
        for r, (a, b, cc, d) in enumerate(meta):
            if r == 0:
                continue
            buf = np.empty((self.nz, b - a, d - cc), dtype=acc_dtype)
            c.Recv([buf, MPI.FLOAT], source=r, tag=7)
            out[:, a:b, cc:d] = buf
        return out


# ==========================================================================
#  Reference smoother (only used by --check)
# ==========================================================================

def reference_smooth(path, varname, sigma, is_meridional, point="T"):
    r"""
    reference_smooth: whole-domain smoothing in one process, the reference a tiled run is checked against.

    The land-aware ratio is written out so that only the field carries the sign of a velocity across the
    fold, never the mask: a sign-flipped mask would blow up the ratio on the rows with little ocean just
    south of the fold. The result does not depend on any tiling.

        Parameters:
            path          (str): NetCDF file
            varname       (str): Variable to smooth
            sigma       (float): Gaussian lengthscale, in grid points
            is_meridional (bool): True if the field's sign flips across the fold
            point         (str): 'T', 'U', 'V' or 'F', where the field lives

        Returns:
            field_out (float): (nz, ny, nx) smoothed field

        Formula:
            \tilde f = \frac{G_\sigma * (f\,m)}{G_\sigma * m}\;m,\ \text{the field with sign } s \text{ at the fold, the mask without}

        See also:
            Tile.smooth, orca.smooth_2D
    """
    data, wet = _read_slab(path, varname)
    num = orca.smooth_2D(data * wet, sigma, isMeridional=is_meridional, point=point)
    den = orca.smooth_2D(wet,        sigma, isMeridional=False, point=point)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(den > 0.01, num / den, 0.0)
    return np.ascontiguousarray(out * wet, dtype=acc_dtype)


# ==========================================================================
#  Driver
# ==========================================================================

def main():
    r"""
    main: command-line entry point: filter one variable, serial or tiled under mpirun.

    Options: --file, --component, --date, --sigma, --meridional, --point, --gather, --check. With --check
    the tiled result is compared with reference_smooth.

        Parameters:
            None; options are read from the command line

        Returns:
            None; prints each rank's tile, the timing and, with --check, the largest and rms
            difference from reference_smooth

        Formula:
            \text{max and rms of } \tilde f_{\text{tiled}} - \tilde f_{\text{reference}} \text{ with --check}

        See also:
            Tile.smooth, reference_smooth
    """
    ap = argparse.ArgumentParser(
        description="MPI-tiled spatial filter sketch (full vertical pencils)."
    )
    ap.add_argument("--file", default=None,
                    help="explicit netCDF path (overrides --component/--date)")
    ap.add_argument("--component", default="vozocrtx",
                    help="component name; also the in-file variable name")
    ap.add_argument("--date", default="2015-01-01", help="date token in the file name")
    ap.add_argument("--sigma", type=float, default=8.0,
                    help="Gaussian sigma in gridpoints (halo width hw = ceil(4*sigma))")
    ap.add_argument("--meridional", action="store_true",
                    help="field's sign flips across the tripolar fold (either velocity component)")
    ap.add_argument("--point", choices=sorted(orca.FOLD), default="T",
                    help="where the field lives")
    ap.add_argument("--gather", action="store_true",
                    help="assemble the filtered field on rank 0")
    ap.add_argument("--check", action="store_true",
                    help="rank 0 also runs reference_smooth and reports max|delta| (implies --gather)")
    args = ap.parse_args()

    path = args.file or default_path(args.component, args.date)
    varname = args.component
    is_meridional = args.meridional
    fname = os.path.basename(path)

    # Serial iff mpi4py is absent OR we were launched with a single process.
    serial = (not _HAVE_MPI) or MPI.COMM_WORLD.Get_size() == 1

    # -------------------------------------------------------------- serial path
    if serial:
        t0 = time.time()
        tile = Tile.from_whole(path, varname, args.sigma, is_meridional, args.point)
        out = tile.smooth(args.sigma)
        print(f"[serial] {fname}  sigma={args.sigma}  hw={tile.hw}  out{out.shape}  "
              f"{time.time() - t0:.1f}s  mean|f|={np.abs(out).mean():.3e}")
        if args.check:
            ref = reference_smooth(path, varname, args.sigma, is_meridional, args.point)
            e = np.abs(out - ref)
            print(f"[check] vs reference_smooth: "
                  f"max={e.max():.3e}  rms={np.sqrt((e ** 2).mean()):.3e}")
        return

    # -------------------------------------------------------------- parallel path
    comm = MPI.COMM_WORLD
    tile = Tile.from_mpi(comm, path, varname, args.sigma, is_meridional, args.point)

    if tile.rank == 0:
        print(f"[mpi] {fname}  N={comm.Get_size()}  grid={tile.Py}x{tile.Px}  "
              f"sigma={args.sigma}  hw={tile.hw}  "
              f"global=({tile.nz},{tile.ny_g},{tile.nx_g})")

    tile.cart_comm.Barrier()                       # so the timing below is comparable across ranks
    t0 = time.time()
    local = tile.smooth(args.sigma)
    dt = tile.cart_comm.reduce(time.time() - t0, op=MPI.MAX, root=0)   # slowest rank's time

    # one line per rank: which tile it owned, and a cheap sanity scalar
    print(f"  rank {tile.rank:3d}  coords=({tile.cy},{tile.cx})  "
          f"tile[z,{tile.j0}:{tile.j0 + tile.jl},{tile.i0}:{tile.i0 + tile.il}]  "
          f"mean|f|={np.abs(local).mean():.3e}")

    if args.gather or args.check:
        g = tile.gather_to_root(local)       # (nz,ny,nx) on rank 0, None elsewhere
        if tile.rank == 0:
            print(f"[mpi] gathered {g.shape} in worst-rank {dt:.1f}s")
            if args.check:
                ref = reference_smooth(path, varname, args.sigma, is_meridional, args.point)
                e = np.abs(g - ref)
                print(f"[check] vs reference_smooth: "
                      f"max={e.max():.3e}  rms={np.sqrt((e ** 2).mean()):.3e}")


if __name__ == "__main__":
    main()
