#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mpi_orca.py
===========

ORCA / eORCA global-ocean domain decomposition + halo / communication layer.

This module is ONLY the infrastructure. It cuts a global `(nz, ny, nx)` grid
into z-pencils (the full water column over a rectangular horizontal tile),
moves the thin overlaps between tiles that a horizontally non-local operation
needs, and assembles per-tile results back into a global array. It contains NO
grid operations (smoothing, regridding, metric reconstruction, binning, ...) --
those belong in a separate module that imports these primitives.

Serial vs tiled
---------------
Every routine takes a `Decomp` as its first argument and branches once on
    `dec.serial`:

    dec = Decomp.for_file(path, varname)   # serial if no MPI or -np 1,
                                           # else a P x P Cartesian tiling,
                                           # P = sqrt(N)

When `dec.serial` is True the tiled routines fall straight through to the plain
whole-array path; when it is a tiled `Decomp` they add exactly the MPI step
(halo exchange, a gather, or an allreduce) that a caller's operation needs.

What's here
-----------
  partitioning
      _split                 -- load-balanced 1-D chunking
      Decomp                 -- topology, this tile's global window, neighbours
  I/O / assembly
      read_slab              -- read this tile's hyperslab (all nz levels)
      read_slab_with_halo    -- read the hyperslab AND its halo, straight off
                                 disk, no communication (only valid for
                                 static/read-only source data)
      gather                 -- per-tile -> global array on rank 0
      Decomp.allreduce_sum   -- sum a per-tile partial result across all ranks
  halo primitives  (E-W wrap + tripolar north fold + closed south)
      tiled  : add_halo / strip_halo / exchange_halo / _north_fold
      serial : halo / dehalo   -- whole-array counterpart, imported from orca

Halo / boundary rules (in `exchange_halo`)
-----------------------------------------
  x  : PERIODIC (east-west wrap), period nx-2: the last two columns repeat the
       first two. The Cartesian x dimension is periodic, so every rank always
       has a real east and west neighbour; the outermost tiles shift what they
       send by those two columns.
  y  : NON-periodic.
         interior tile faces  : ordinary neighbour exchange.
         global SOUTH (cy = 0) : closed -- the interior rows reflected, the
                                 same rule scipy applies in `orca.smooth_2D`.
         global NORTH (cy=Py-1): T-pivot fold, as `orca.halo`: the ghosts come
                                 from the source row and column of `orca.FOLD`,
                                 which depend on the grid point, times fold_sign
                                 (-1 for either velocity component, else +1, and
                                 always +1 for a wet mask). The top tile-row
                                 pools its strips with an Allgatherv on a row
                                 sub-communicator, since the fold spans the full
                                 global width. The corners are filled too.

Padded local array (what `add_halo` returns)
-------------------------------------------
  shape  = leading_dims + (jl + 2*hw, il + 2*hw)     (jl, il = interior sizes)
  interior     a[..., hw:hw+jl, hw:hw+il]     -- this rank's own data
  west  halo   a[..., :,        0:hw]
  east  halo   a[..., :,        hw+il:]
  south halo   a[..., 0:hw,     :]
  north halo   a[..., hw+jl:,   :]
  NB `add_halo` uses CENTRED padding (a south halo exists); the whole-array
  `halo()` primitive does not. For exact `halo()` semantics call `halo()`.

Sketch-level simplifications
---------------------------
  1. hw + 5 <= the smallest local tile dimension: the halo is filled from the
     immediate neighbour only (no multi-hop), and the fold reaches 5 rows further
     down. Checked in `add_halo`; raise if violated -- use fewer ranks or a
     narrower halo.
  2. Each rank opens the netCDF independently (no parallel HDF5) -- fine for
     read-only; a caller that writes a result back does so through `gather`
     to rank 0.

Usage
-----
    dec       = Decomp.for_file(path, "vozocrtx")
    sub, wet  = read_slab(dec, path, "vozocrtx")     # this tile (or whole grid)
    a         = add_halo(dec, sub, hw)               # padded, halo filled
    # ... run a horizontally non-local op on `a` (in another module) ...
    full      = gather(dec, strip_halo(result, hw))  # (nz, ny, nx) on rank 0
"""

import sys

import numpy as np
import netCDF4 as netcdf

from .orca import FOLD, dehalo, fold_columns, halo, wrap_columns

# mpi4py is OPTIONAL. Missing import, or a single-process launch, => serial.
try:
    from mpi4py import MPI
    _HAVE_MPI = True
except Exception:
    _HAVE_MPI = False

# Accumulator dtype for field data read off disk (see `read_slab`).
acc_dtype = np.float32

# north-fold halo directions, kept for callers that use the whole-array `halo`
whalo = np.array([1, 0, 0])
nhalo = np.array([0, 1, 0])
ehalo = np.array([0, 0, 1])


def _split(n, p):
    r"""
    _split: cut range(n) into p contiguous chunks of nearly equal size.

    Places the tiles of a decomposition in the global grid, as NEMO splits its domain. The first
    (n mod p) chunks get one extra cell, as numpy.array_split does: _split(1021, 3) -> [(0, 341),
    (341, 681), (681, 1021)].

        Parameters:
            n (int): Global size along one axis
            p (int): Number of chunks, the tiles along that axis

        Returns:
            out (list): p (start, stop) tuples

        Formula:
            n = q\,p + r,\qquad \text{chunk } k \text{ has } q + [k < r] \text{ cells},\qquad \text{start}_k = k\,q + \min(k, r)

        NEMO:
            src/OCE/LBC/mppini.F90: SUBROUTINE mpp_basesplit( kiglo, kjglo, khls, knbi, knbj, kimax, kjmax, kimppt, kjmppt, klci, klcj)

        See also:
            Decomp.from_comm
    """
    base, rem = divmod(n, p)
    out, s = [], 0
    for k in range(p):
        e = s + base + (1 if k < rem else 0)
        out.append((s, e))
        s = e
    return out


# ==========================================================================
#  DECOMPOSITION
# ==========================================================================

class Decomp:
    """The domain decomposition: MPI Cartesian topology + this tile's place in
    the global grid. Field-agnostic -- build it once, pass it to every routine.

    Build one with a constructor below; never call `Decomp(...)` directly:
        Decomp.for_file(path, var)  -- recommended. Reads the grid size, then
                                       serial if no MPI / -np 1, else P x P.
        Decomp.auto(ny_g, nx_g)     -- same choice, sizes given explicitly.
        Decomp.from_comm(comm, ...) -- force a tiling over a given communicator.
        Decomp.serial_decomp(...)   -- force the whole-grid (serial) case.

    Attributes
        serial            : True  => one tile == the whole global grid
        cart_comm         : Cartesian communicator (None when serial)
        Px, Py            : tiles per axis (x, y); both = sqrt(N)
        cx, cy            : this tile's (col, row) index in the P x P grid
        i0, j0            : global index of this tile's first interior col / row
        il, jl            : this tile's interior sizes (cols, rows)
        nx_g, ny_g        : global horizontal sizes
        rank, nranks      : within `cart_comm`
        row_comm          : sub-communicator of the tile-row (varies cx), or None
        min_tile          : smallest interior dim over ALL ranks (halo cap)
    """

    def __init__(self, cart_comm, Px, Py, cx, cy, i0, j0, jl, il, ny_g, nx_g):
        r"""
        __init__: low-level constructor; use for_file, auto, from_comm or serial_decomp instead.

        Records where this tile lies in the global grid and in the P x P layout.

            Parameters:
                cart_comm (Comm): Cartesian communicator, or None for serial
                Px, Py     (int): Tiles per axis, both sqrt(N), or 1 when serial
                cx, cy     (int): This tile's column and row in the P x P grid
                i0, j0     (int): Global index of its first interior column and row
                jl, il     (int): Its interior sizes, rows and columns
                ny_g, nx_g (int): Global horizontal sizes

            Returns:
                None; sets the attributes listed in the class docstring

            Formula:
                \text{tile } (c_x, c_y) \text{ covers } j \in [j_0, j_0 + j_l),\ i \in [i_0, i_0 + i_l)

            See also:
                Decomp.from_comm, Decomp.serial_decomp
        """
        self.cart_comm = cart_comm
        self.serial = cart_comm is None
        self.Px, self.Py = Px, Py
        self.cx, self.cy = cx, cy
        self.i0, self.j0 = i0, j0
        self.jl, self.il = jl, il
        self.ny_g, self.nx_g = ny_g, nx_g
        self.rank = 0 if self.serial else cart_comm.Get_rank()
        self.nranks = 1 if self.serial else cart_comm.Get_size()
        self.row_comm = None if self.serial else cart_comm.Sub([False, True])
        self.min_tile = min(jl, il) if self.serial else cart_comm.allreduce(min(jl, il), op=MPI.MIN)

    # ---- constructors ----
    @classmethod
    def serial_decomp(cls, ny_g, nx_g):
        r"""
        serial_decomp: the whole-grid case, one tile spanning the entire grid.

        Code written for tiles then runs unchanged in one process: every routine of this module falls back to
        its whole-array form.

            Parameters:
                ny_g, nx_g (int): Global horizontal sizes

            Returns:
                dec (Decomp): Serial decomposition, .serial True and .cart_comm None

            Formula:
                P_x = P_y = 1,\qquad (j_0, i_0) = (0, 0),\qquad (j_l, i_l) = (n_y, n_x)

            See also:
                Decomp.auto
        """
        return cls(None, 1, 1, 0, 0, 0, 0, ny_g, nx_g, ny_g, nx_g)

    @classmethod
    def from_comm(cls, comm, ny_g, nx_g):
        r"""
        from_comm: a P x P tiling of the grid over the ranks of a communicator.

        The grid is cut into P strips in each direction, one tile per rank. The x axis is periodic and the y
        axis is not, as on the ORCA grids. MPI may renumber ranks, so each tile's position comes from its
        Cartesian coordinates.

        Exits the program if the number of ranks is not a perfect square.

            Parameters:
                comm      (Comm): MPI communicator of size P*P
                ny_g, nx_g (int): Global horizontal sizes

            Returns:
                dec (Decomp): Decomposition for this rank's tile

            Formula:
                P = \sqrt{N},\qquad [j_0, j_0 + j_l) = \text{\_split}(n_y, P)[c_y],\qquad [i_0, i_0 + i_l) = \text{\_split}(n_x, P)[c_x]

            NEMO:
                src/OCE/LBC/mppini.F90: SUBROUTINE mpp_basesplit( kiglo, kjglo, khls, knbi, knbj, kimax, kjmax, kimppt, kjmppt, klci, klcj)

            See also:
                _split, Decomp.auto
        """
        N = comm.Get_size()
        P = int(round(N ** 0.5))
        if P * P != N:
            if comm.Get_rank() == 0:
                sys.stderr.write(
                    f"mpi_orca: need a perfect-square process count, got {N}; "
                    f"try one of {[k * k for k in range(2, 12)]}\n")
            comm.Barrier()
            sys.exit(1)
        # dims [Py, Px]; axis 0 = y (real edges), axis 1 = x (periodic wrap);
        # reorder => MPI may renumber ranks, so always read Get_coords.
        cart_comm = comm.Create_cart([P, P], periods=[False, True], reorder=True)
        cy, cx = cart_comm.Get_coords(cart_comm.Get_rank())
        j0, j1 = _split(ny_g, P)[cy]
        i0, i1 = _split(nx_g, P)[cx]
        return cls(cart_comm, P, P, cx, cy, i0, j0, j1 - j0, i1 - i0, ny_g, nx_g)

    @classmethod
    def auto(cls, ny_g, nx_g):
        r"""
        auto: the decomposition suited to how the program was started.

        Serial if mpi4py is missing or the job has one process, otherwise a P x P tiling of all ranks, so the
        same script runs with or without mpirun.

            Parameters:
                ny_g, nx_g (int): Global horizontal sizes

            Returns:
                dec (Decomp): Serial or tiled decomposition, chosen from the run

            Formula:
                \text{serial if } N = 1 \text{ or no MPI, else from\_comm}(\mathrm{COMM\_WORLD})

            See also:
                Decomp.from_comm, Decomp.serial_decomp
        """
        if (not _HAVE_MPI) or MPI.COMM_WORLD.Get_size() == 1:
            return cls.serial_decomp(ny_g, nx_g)
        return cls.from_comm(MPI.COMM_WORLD, ny_g, nx_g)

    @classmethod
    def for_file(cls, path, varname):
        r"""
        for_file: Decomp.auto with the grid size read from a file; the usual entry point.

        Saves passing the grid size by hand.

            Parameters:
                path    (str): NetCDF file
                varname (str): Name of a 3D (z, y, x) variable in it

            Returns:
                dec (Decomp): Serial or tiled decomposition for that grid

            Formula:
                (n_y, n_x) = \text{shape of varname}[1:],\ \text{then auto}(n_y, n_x)

            See also:
                Decomp.auto
        """
        with netcdf.Dataset(path) as nc:
            _, ny_g, nx_g = nc[varname].shape
        return cls.auto(ny_g, nx_g)

    # ---- geometry helpers ----
    def slices(self):
        r"""
        slices: the window of a global array that belongs to this tile.

        For reading this tile's part of a file or of a global array.

            Parameters:
                None

            Returns:
                jsl (slice): Rows of the tile; slice(None) when serial
                isl (slice): Columns of the tile; slice(None) when serial

            Formula:
                \text{jsl} = [j_0, j_0 + j_l),\qquad \text{isl} = [i_0, i_0 + i_l)

            See also:
                read_slab
        """
        if self.serial:
            return slice(None), slice(None)
        return slice(self.j0, self.j0 + self.jl), slice(self.i0, self.i0 + self.il)

    @property
    def is_north_edge(self):
        """
        is_north_edge: whether this tile touches the global north edge, the fold.

        Tiles on that edge fill their north halo from the fold rather than from a neighbour. Always true
        when serial.

            Parameters:
                None

            Returns:
                edge (bool): True on the last row of tiles, and always when serial

            Formula:
                c_y = P_y - 1

            See also:
                Decomp.is_south_edge
        """
        return self.serial or self.cy == self.Py - 1

    @property
    def is_south_edge(self):
        """
        is_south_edge: whether this tile touches the global south edge, which is closed.

        Tiles on that edge reflect their south halo. Always true when serial.

            Parameters:
                None

            Returns:
                edge (bool): True on the first row of tiles, and always when serial

            Formula:
                c_y = 0

            See also:
                Decomp.is_north_edge
        """
        return self.serial or self.cy == 0

    def barrier(self):
        r"""
        barrier: wait until every rank reaches this point.

        For lining up ranks, for example before one of them reads a file another has just written. Does
        nothing when serial.

            Parameters:
                None

            Returns:
                None

            Formula:
                \text{all ranks wait until every rank arrives}

            NEMO:
                src/OCE/LBC/lib_mpp.F90: CALL mpi_barrier( mpi_comm_oce, ierror )

            See also:
                Decomp.allreduce_sum
        """
        if not self.serial:
            self.cart_comm.Barrier()

    def allreduce_sum(self, arr):
        r"""
        allreduce_sum: element-wise sum of an array over all ranks.

        When several tiles contribute to the same cells, such as the numerator and denominator of an average
        onto a coarser grid, their partial results must be summed before use. The array has the same shape
        on every rank. Does nothing when serial.

            Parameters:
                arr (array): Partial result, same shape on every rank

            Returns:
                total (array): Sum over all ranks, reduced in place (on a contiguous copy if needed);
                               arr unchanged when serial

            Formula:
                a \leftarrow \sum_{r=0}^{N-1} a^{(r)}

            See also:
                gather
        """
        if self.serial:
            return arr
        buf = np.ascontiguousarray(arr)
        self.cart_comm.Allreduce(MPI.IN_PLACE, buf, op=MPI.SUM)
        return buf


# ==========================================================================
#  I/O
# ==========================================================================

def read_slab(dec, path, varname, level=None):
    r"""
    read_slab: read this rank's window of a variable from disk.

    Each rank opens the file on its own and reads only its tile; no parallel file system library is
    needed.

    Missing and fill values are returned as 0, with the mask of valid points beside them.

        Parameters:
            dec  (Decomp): Decomposition; serial reads the whole grid
            path    (str): NetCDF file
            varname (str): 3D (z, y, x) variable
            level   (int): One level, keeping a leading axis of size 1, or None for all levels

        Returns:
            data (float): (nz, jl, il), or (nz, ny, nx) when serial; float32, NaN and fill as 0
            wet  (float): Same shape, 1 where the value read was finite, else 0

        Formula:
            d = f[:, \text{jsl}, \text{isl}] \text{ with NaN} \to 0,\qquad w = [f \text{ finite}]

        See also:
            read_slab_with_halo, Decomp.slices
    """
    lvl = slice(None) if level is None else slice(level, level + 1)
    jsl, isl = dec.slices()
    with netcdf.Dataset(path) as nc:
        raw = np.ma.masked_invalid(nc[varname][lvl, jsl, isl])
    filled = np.ma.filled(raw, 0.0)
    data = np.ascontiguousarray(filled, dtype=acc_dtype)
    wet = np.ascontiguousarray(
        ~np.ma.getmaskarray(raw) & np.isfinite(filled), dtype=acc_dtype)
    return data, wet


def read_slab_with_halo(dec, path, varname, hw, is_meridional=False, level=None, point="T"):
    r"""
    read_slab_with_halo: read this rank's tile and its halo from disk, with no communication between ranks.

    When the data are a file that every rank can open, a slightly wider window can be read in place of
    exchanging the halo with the neighbours, which is faster in a loop over many files. The same boundary
    rules as exchange_halo apply: columns wrap with period nx - 2 across the seam, rows below 0 are
    reflected, and rows beyond the top come from the fold, sign-flipped for a velocity component.

    Unlike add_halo, the halo is not limited by the tile size. Missing values are returned as 0, with
    the mask beside them.

        Parameters:
            dec        (Decomp): Decomposition, serial or tiled
            path          (str): NetCDF file with a 3D (z, y, x) variable
            varname       (str): Variable to read
            hw            (int): Halo width; unlike add_halo, not limited by the tile size
            is_meridional (bool): True for a field whose sign flips across the fold, either
                                  velocity component; only the fold halo is negated
            level         (int): One level, keeping a leading axis of size 1, or None for all
            point         (str): 'T', 'U', 'V' or 'F', where the field lives

        Returns:
            a (float): (nz, jl+2hw, il+2hw) float32 field with halo, or (nz, ny_g+2hw, nx_g+2hw)
                       when serial; NaN and fill as 0
            w (float): Same shape, 1 where finite, else 0; never sign-flipped

        Formula:
            a_{j,i} = f_{j,\ i \bmod (n_x-2)} \text{ in x};\qquad a_{-1-p,\,i} = a_{p,\,i} \text{ south (reflection)};\qquad a_{n_y+p,\,i} = s\,f_{n_y-r-p,\ \mathrm{fold}(i)} \text{ north}

        NEMO:
            src/OCE/LBC/lbc_nfd_generic.h90: ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)

        See also:
            read_slab, exchange_halo, orca.halo
    """
    j0, jl, i0, il = dec.j0, dec.jl, dec.i0, dec.il
    ny_g, nx_g = dec.ny_g, dec.nx_g
    lvl = slice(None) if level is None else slice(level, level + 1)  # keeps a leading axis of size 1

    with netcdf.Dataset(path) as nc:
        var = nc[varname]

        # ---- x: periodic, possibly wrapping ----
        columns = np.arange(i0 - hw, i0 + il + hw)   # halo'd column window, global index space
        idx = wrap_columns(columns, nx_g)
        runs = np.split(idx, np.where(np.diff(idx) != 1)[0] + 1)   # contiguous runs

        def xread(rows):
            r"""
            xread: read the halo'd columns of some rows, in the contiguous runs that wrap across the seam.

            A helper of read_slab_with_halo: a window crossing the seam is read in up to three contiguous parts.

                Parameters:
                    rows (slice): Rows to read

                Returns:
                    part (MaskedArray): (nz, rows, il+2hw) data

                Formula:
                    \text{columns } [i_0 - h, i_0 + i_l + h) \bmod (n_x - 2), \text{ read as contiguous runs}

                See also:
                    read_slab_with_halo
            """
            parts = [var[lvl, rows, r[0]:r[-1] + 1] for r in runs]
            return np.ma.concatenate(parts, axis=-1)

        # ---- y: south closed / interior / north fold ----
        y_lo, y_hi = j0 - hw, j0 + jl + hw
        y_read_lo, y_read_hi = max(y_lo, 0), min(y_hi, ny_g)
        core = np.ma.masked_invalid(xread(slice(y_read_lo, y_read_hi)))
        nz = core.shape[0]
        south_pad = y_read_lo - y_lo             # >0 only for the tile-row on the global south edge
        north_pad = y_hi - y_read_hi              # >0 only for the tile-row on the global north edge

        a = np.ma.masked_all((nz, jl + 2 * hw, il + 2 * hw), dtype=acc_dtype)
        a[:, south_pad:south_pad + core.shape[1], :] = core
        if south_pad > 0:
            a[:, :south_pad, :] = a[:, 2 * south_pad - 1:south_pad - 1:-1, :]   # reflect
        if north_pad > 0:
            # ascending, contiguous disk read; reversing puts block row p at
            # global row ny_g-row-p, the fold's source row for ghost p.
            top = ny_g - FOLD[point][0] + 1
            fold_asc = np.ma.masked_invalid(var[lvl, top - north_pad:top, :])   # (nz, north_pad, nx_g)
            fold_sel = fold_asc[:, ::-1, :][:, :, fold_columns(columns, nx_g, point)]
            sign = -1.0 if is_meridional else 1.0
            a[:, jl + 2 * hw - north_pad:, :] = sign * fold_sel

    filled = np.ma.filled(a, 0.0)
    data = np.ascontiguousarray(filled, dtype=acc_dtype)
    wet = np.ascontiguousarray(
        ~np.ma.getmaskarray(a) & np.isfinite(filled), dtype=acc_dtype)
    return data, wet


def gather(dec, local, root=0):
    r"""
    gather: assemble the tiles of a field into one global array on one rank.

    Needed to write a result or check it against a whole-grid computation. The whole field can be
    hundreds of MB, so gather only when it is needed; choosing a different root for each field spreads
    the writing across ranks.

        Parameters:
            dec  (Decomp): Decomposition the tiles were made under
            local (float): This rank's (..., jl, il) tile, or (..., ny, nx) when serial
            root    (int): Rank that receives the field

        Returns:
            field (float): (..., ny_g, nx_g) on root, None on the other ranks; local itself when serial

        Formula:
            F[\ldots, j_0^{(r)} : j_0^{(r)} + j_l^{(r)},\ i_0^{(r)} : i_0^{(r)} + i_l^{(r)}] = f^{(r)} \quad \text{for every rank } r

        See also:
            Decomp.allreduce_sum
    """
    if dec.serial:
        return local

    c = dec.cart_comm
    lead = local.shape[:-2]
    meta = c.gather((dec.j0, dec.j0 + dec.jl, dec.i0, dec.i0 + dec.il), root=root)
    mt = MPI.FLOAT if local.dtype == np.float32 else MPI.DOUBLE

    if dec.rank != root:
        c.Send([np.ascontiguousarray(local), mt], dest=root, tag=7)
        return None

    out = np.zeros(lead + (dec.ny_g, dec.nx_g), dtype=local.dtype)
    out[..., dec.j0:dec.j0 + dec.jl, dec.i0:dec.i0 + dec.il] = local
    for r, (a, b, cc, d) in enumerate(meta):
        if r == root:
            continue
        buf = np.empty(lead + (b - a, d - cc), dtype=local.dtype)
        c.Recv([buf, mt], source=r, tag=7)
        out[..., a:b, cc:d] = buf
    return out


# ==========================================================================
#  HALO EXCHANGE  (the core of the tiled mode)
# ==========================================================================

def _north_fold(dec, a, hw, fold_sign, point="T"):
    r"""
    _north_fold: fill the north halo of a padded array from the T-pivot fold, in place.

    Beyond the top row the ORCA grid folds back on itself: each ghost row is a real row further south,
    read in reverse, with the sign flipped for a velocity component, as NEMO's north-fold condition does.
    In parallel the tiles of the top row pool their strips to find the mirrored columns.

        Parameters:
            dec     (Decomp): Decomposition, serial or tiled
            a        (float): Padded array, interior at a[..., hw:hw+jl, hw:hw+il]
            hw         (int): Halo width
            fold_sign (float): +1 for scalars, variances and masks, -1 for velocity components
            point      (str): 'T', 'U', 'V' or 'F', where the field lives

        Returns:
            None; a's north halo rows are written in place

        Formula:
            a_{n_y+p,\,i} = s\;a_{n_y-r-p,\ \mathrm{fold}(i)},\qquad p = 0, \dots, h-1

        NEMO:
            src/OCE/LBC/mpp_nfd_generic.h90: SUBROUTINE mpp_nfd_/**/PRECISION( cdname, ptab, cd_nat, psgn, kfillmode, pfillval, kfld, ldfull )
            src/OCE/LBC/lbc_nfd_generic.h90: ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)

        See also:
            exchange_halo, orca.fold_columns
    """
    lead = a.shape[:-2]
    jl = a.shape[-2] - 2 * hw
    il = a.shape[-1] - 2 * hw
    row = FOLD[point][0]
    depth = hw + row                                       # block row s <-> global row ny_g-depth+s
    contrib = np.ascontiguousarray(a[..., hw + jl - depth:hw + jl, hw:hw + il])  # lead+(depth, il)

    if dec.serial:
        strip = contrib                                    # already global in x
    else:
        rc = dec.row_comm
        mt = MPI.FLOAT if a.dtype == np.float32 else MPI.DOUBLE
        counts = np.array(rc.allgather(contrib.size), dtype=int)
        widths = np.array(rc.allgather(il), dtype=int)
        displ = np.concatenate(([0], np.cumsum(counts)[:-1]))
        recv = np.empty(int(counts.sum()), dtype=a.dtype)
        rc.Allgatherv([contrib.ravel(), mt], [recv, counts, displ, mt])
        parts = [recv[d:d + n].reshape(lead + (depth, w))
                 for n, d, w in zip(counts, displ, widths)]
        strip = np.concatenate(parts, axis=-1)             # lead+(depth, nx_g)

    source = fold_columns(np.arange(dec.i0 - hw, dec.i0 + il + hw), dec.nx_g, point)
    for p in range(hw):
        a[..., hw + jl + p, :] = fold_sign * strip[..., depth - row - p, source]


def exchange_halo(dec, a, hw, is_meridional, point="T"):
    r"""
    exchange_halo: fill the halo of a padded array, in place.

    A tile needs a margin of its neighbours' values for any operation that looks at nearby points, such
    as a filter or a gradient. This fills it as NEMO's halo exchange does: east-west wrap across the seam,
    copies from the neighbouring tiles, reflection at the closed south edge and the fold at the north.

    The same call serves a field or its mask; only the sign at the fold differs. add_halo allocates the
    padded array and calls this.

        Parameters:
            dec        (Decomp): Serial (one tile spans the globe) or tiled (exchange with
                                 neighbours)
            a           (float): Padded array, interior at a[..., hw:hw+jl, hw:hw+il]
            hw            (int): Halo width
            is_meridional (bool): True for a field whose sign flips across the fold, either
                                  velocity component
            point         (str): 'T', 'U', 'V' or 'F', where the field lives

        Returns:
            None; a's halo is written in place

        Formula:
            \text{west/east: wrap with period } n_x - 2;\quad \text{south: } a_{-1-p} = a_p;\quad \text{north: fold};\quad \text{between tiles: copy of the neighbour's edge}

        NEMO:
            src/OCE/LBC/lbc_lnk_pt2pt_generic.h90: SUBROUTINE lbc_lnk_pt2pt_/**/PRECISION( cdname, ptab, cd_nat, psgn, kfld, kfillmode, pfillval, lsend, lrecv, ld4only, ldfull )

        See also:
            add_halo, _north_fold
    """
    lead_ok = a.ndim >= 2
    assert lead_ok
    jl = a.shape[-2] - 2 * hw
    il = a.shape[-1] - 2 * hw
    fold_sign = -1.0 if is_meridional else 1.0

    # ---- SERIAL: the single tile already spans the globe ----
    if dec.serial:
        a[..., :, 0:hw]     = a[..., :, il - 2:il + hw - 2]     # west  <- east interior (period nx-2)
        a[..., :, hw + il:] = a[..., :, hw + 2:2 * hw + 2]      # east  <- west interior (period nx-2)
        a[..., 0:hw, :]     = a[..., 2 * hw - 1:hw - 1:-1, :]   # south: reflect
        _north_fold(dec, a, hw, fold_sign, point)               # north: tripolar fold
        return

    # ---- PARALLEL: point-to-point exchange ----
    c = dec.cart_comm
    mt = MPI.FLOAT if a.dtype == np.float32 else MPI.DOUBLE

    # x (Cartesian axis 1), PERIODIC. Shift(dir, disp) -> (source, dest):
    #   dest   = neighbour at coord[dir] + disp   (I send toward here)
    #   source = neighbour at coord[dir] - disp   (I receive from here)
    src_w, dst_e = c.Shift(1, +1)
    # The last two global columns repeat the first two, so across the seam the
    # outermost tiles send the columns just inside them.
    shift_e = 2 if dec.cx == dec.Px - 1 else 0
    shift_w = 2 if dec.cx == 0 else 0
    send_e = np.ascontiguousarray(a[..., :, il - shift_e:il + hw - shift_e])   # my last hw interior cols, full height
    recv_w = np.empty_like(send_e)
    c.Sendrecv([send_e, mt], dest=dst_e, sendtag=0, recvbuf=[recv_w, mt], source=src_w, recvtag=0)
    a[..., :, 0:hw] = recv_w

    src_e, dst_w = c.Shift(1, -1)
    send_w = np.ascontiguousarray(a[..., :, hw + shift_w:2 * hw + shift_w])    # my first hw interior cols
    recv_e = np.empty_like(send_w)
    c.Sendrecv([send_w, mt], dest=dst_w, sendtag=1, recvbuf=[recv_e, mt], source=src_e, recvtag=1)
    a[..., :, hw + il:] = recv_e
    # distinct tags so the two exchanges to the SAME partner (P = 2) cannot
    # cross. x strips span the full padded height => the corners get filled
    # when the y exchange below copies whole padded-width strips.

    # y (Cartesian axis 0), NON-periodic.
    src_s, dst_n = c.Shift(0, +1)
    send_n = np.ascontiguousarray(a[..., hw + jl - hw:hw + jl, :])  # my last hw interior rows, full width
    recv_s = np.empty_like(send_n)
    c.Sendrecv([send_n, mt], dest=dst_n, sendtag=2, recvbuf=[recv_s, mt], source=src_s, recvtag=2)
    if src_s != MPI.PROC_NULL:            # PROC_NULL => global south edge; recv_s is garbage
        a[..., 0:hw, :] = recv_s

    src_n, dst_s = c.Shift(0, -1)
    send_s = np.ascontiguousarray(a[..., hw:2 * hw, :])            # my first hw interior rows
    recv_n = np.empty_like(send_s)
    c.Sendrecv([send_s, mt], dest=dst_s, sendtag=3, recvbuf=[recv_n, mt], source=src_n, recvtag=3)
    if src_n != MPI.PROC_NULL:            # PROC_NULL => global north edge; fold fills it instead
        a[..., hw + jl:, :] = recv_n

    # global y boundaries (outermost tile-rows only)
    if dec.cy == 0:
        a[..., 0:hw, :] = a[..., 2 * hw - 1:hw - 1:-1, :]   # closed south: reflect
    if dec.cy == dec.Py - 1:
        _north_fold(dec, a, hw, fold_sign, point)           # tripolar fold


def add_halo(dec, field, hw, is_meridional=False, point="T"):
    r"""
    add_halo: a padded copy of a tile, with its halo filled from the neighbours.

    The building block of operations on tiles that look at nearby points, such as a filter: pad, apply,
    then strip_halo.

    The halo comes from one neighbour, and the fold reaches 5 rows further down, so the halo width plus
    5 must fit in the smallest tile.

        Parameters:
            dec        (Decomp): Decomposition
            field       (float): (..., y, x) tile, or the whole grid when serial
            hw            (int): Halo width; hw + 5 must be at most dec.min_tile, since the halo
                                 comes from one neighbour and the fold reaches 5 rows further down,
                                 else ValueError
            is_meridional (bool): True for a field whose sign flips across the fold, either
                                  velocity component
            point         (str): 'T', 'U', 'V' or 'F', where the field lives

        Returns:
            a (float): (..., jl+2hw, il+2hw) padded array, interior at a[..., hw:hw+jl, hw:hw+il]

        Formula:
            a[\ldots, h:h+j_l,\ h:h+i_l] = f,\ \text{halo from exchange\_halo}

        See also:
            exchange_halo, strip_halo
    """
    if hw + 5 > dec.min_tile:
        raise ValueError(
            f"add_halo: hw+5={hw + 5} exceeds the smallest local tile dim "
            f"{dec.min_tile}; use fewer ranks or a narrower halo.")
    lead = field.shape[:-2]
    ny, nx = field.shape[-2:]
    a = np.zeros(lead + (ny + 2 * hw, nx + 2 * hw), dtype=field.dtype)
    a[..., hw:hw + ny, hw:hw + nx] = field
    exchange_halo(dec, a, hw, is_meridional, point)
    return a


def strip_halo(a, hw):
    r"""
    strip_halo: the interior of a padded array, the inverse of add_halo's padding.

    Returns a contiguous copy, ready to be gathered or written.

        Parameters:
            a (float): Padded array from add_halo
            hw  (int): Halo width used

        Returns:
            interior (float): a[..., hw:-hw, hw:-hw], contiguous

        Formula:
            f = a[\ldots, h:-h,\ h:-h]

        See also:
            add_halo
    """
    return np.ascontiguousarray(a[..., hw:-hw, hw:-hw])

