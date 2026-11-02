"""Shared file search and reading. Use it through the per-dataset modules."""
import glob
import os
import re

import netCDF4
import numpy as np

DATE = re.compile(r"_(\d{4})-?(\d{2})(?:-(\d{2}))?(?:_|\.nc$)")
# Common names of the fields. Source names keep working.
NAMES = {"u": "vozocrtx", 
         "v": "vomecrty", 
         "w": "vovecrtz", 
         "ssh": "sossheig", 
         "sst": "sosstsst", 
         "temp": "votemper",
         "ubt": "vobtcrtx", 
         "vbt": "vobtcrty"}
# NEMO's barotropic velocity, stored in the file of the full-depth variable.
BAROTROPIC = {"vozocrtx": "vobtcrtx", 
              "vomecrty": "vobtcrty"}
# Where each field lives on the C grid, and whether it flips sign at the north fold.
FIELDS = {"vozocrtx": ("U", True), 
          "vomecrty": ("V", True), 
          "vobtcrtx": ("U", True), 
          "vobtcrty": ("V", True),
          "vovecrtz": ("T", False), 
          "uv_tvcov": ("T", False), 
          "sossheig": ("T", False), 
          "sosstsst": ("T", False), 
          "votemper": ("T", False)}


def file_date(path):
    r"""
    file_date: the date written in a file name.

    Files of every dataset carry their date in the name, daily or monthly, in a few spellings. This reads
    it in one form, so files of different variables can be matched by date.

    Recognises _YYYY-MM-DD, _YYYY-MM and _YYYYMM before an underscore or .nc.

        Parameters:
            path    (str): File path; the name holds _YYYY-MM-DD, _YYYY-MM or _YYYYMM

        Returns:
            date    (str): 'YYYY-MM-DD' for daily files, 'YYYY-MM' for monthly ones

        Formula:
            \text{date} = \texttt{YYYY-MM-DD} \text{ or } \texttt{YYYY-MM} \text{ from } \texttt{\_YYYY[-]MM[-DD]} \text{ in the name}

        See also:
            find_files, Source.groups
    """
    y, m, d = DATE.search(os.path.basename(path)).groups()
    return f"{y}-{m}-{d}" if d else f"{y}-{m}"


def find_files(root, var, start=None, end=None):
    r"""
    find_files: the files of one variable in a dataset folder, between two dates.

    Every dataset stores one file per variable and date under <root>/<var>/. This lists those of one
    variable in date order, accepting common names such as u or ssh as well as NEMO's, and a barotropic
    name, whose values sit in the file of the full-depth variable.

    Dates compare on their common length, so '2015' selects a whole year of daily or monthly files.
    Raises FileNotFoundError if the folder holds none.

        Parameters:
            root    (str): Dataset folder
            var     (str): Common or NEMO name of the variable
            start   (str): First date, inclusive, or None for no limit
            end     (str): Last date, inclusive, or None for no limit

        Returns:
            files  (list): Paths sorted by date. Raises FileNotFoundError if the folder has none

        Formula:
            \{\, f \in \text{root}/v/v\_*.nc \ :\ \text{start} \le \text{date}(f) \le \text{end} \,\}, \text{ sorted by date}

        See also:
            file_date, Source.files
    """
    var = NAMES.get(var, var)
    var = {bt: full for full, bt in BAROTROPIC.items()}.get(var, var)
    files = glob.glob(os.path.join(root, var, f"{var}_*.nc"))
    if not files:
        raise FileNotFoundError(f"no {var} files in {root}")
    out = []
    for f in files:
        date = file_date(f)
        if start and date[:len(start)] < start[:len(date)]:
            continue
        if end and date[:len(end)] > end[:len(date)]:
            continue
        out.append(f)
    return sorted(out, key=file_date)


def clean(field, limit=10.0, fill=np.nan):
    r"""
    clean: mark unusable values, such as fill values left in raw files.

    Raw files sometimes keep a huge fill value, or infinities, where a NaN was meant. Any value that is
    not finite or larger than limit in size is replaced, by NaN unless another fill is asked for.

        Parameters:
            field (float): Array of any shape
            limit (float): Largest absolute value kept
            fill  (float): Value put in place of the others, NaN by default

        Returns:
            field_out (float): field, with non-finite values and values beyond limit replaced by fill

        Formula:
            \text{clean}(x) = x \text{ if } |x| \le \text{limit}, \text{ else fill}

        See also:
            read_one
    """
    field = np.asarray(field)
    with np.errstate(invalid="ignore"):
        return np.where(np.isfinite(field) & (np.abs(field) <= limit), field, fill)


def read_one(path, var, level=None, box=None, dtype=np.float32):
    r"""
    read_one: one snapshot of a variable from one file.

    The datasets differ in how they lay out a file: with or without a time axis, and with depth called
    deptht, depthu, depthv or depthw. This reads any of them into the same shape, land as NaN, with one
    level or a box of grid points if asked.

    The time axis, if any, is reduced to its first entry.

        Parameters:
            path    (str): NetCDF file
            var     (str): Common or NEMO name of the variable
            level   (int): One level, or None for all levels; ignored for 2D fields
            box   (tuple): (j0, j1, i0, i1) index box, or None for the whole grid
            dtype (dtype): Type of the result

        Returns:
            field (float): (ny, nx) or (nz, ny, nx), with masked points (land, fill) as NaN

        Formula:
            x = \text{file}[v][t{=}0,\ \text{level},\ j_0{:}j_1,\ i_0{:}i_1],\ \text{masked} \to \text{NaN}

        See also:
            read_stack, Source.snapshot
    """
    var = NAMES.get(var, var)
    j0, j1, i0, i1 = box if box else (None, None, None, None)
    with netCDF4.Dataset(path) as nc:
        v = nc[var]
        idx = []
        for name in v.dimensions:
            if name == "time_counter":
                idx.append(0)
            elif name.startswith("dep"):
                idx.append(slice(None) if level is None else level)
            elif name == "y":
                idx.append(slice(j0, j1))
            else:
                idx.append(slice(i0, i1))
        data = v[tuple(idx)]
    return np.ma.filled(np.ma.asarray(data, dtype=dtype), np.nan)


def stream_files(file_lists, names, level=0, box=None):
    r"""
    stream_files: read several variables snapshot by snapshot, checking their dates match.

    Long records do not fit in memory; reading one date at a time does. The variables are read together
    so each step gets, for example, u and v of the same day.

    Raises ValueError if the lists differ in length or a step's files differ in date.

        Parameters:
            file_lists (list): One sorted list of files per variable, all of the same length
            names      (list): Variable name of each list
            level       (int): See read_one
            box       (tuple): See read_one

        Returns:
            generator: (date, field, field, ...) per snapshot. Raises ValueError if the lists
                       differ in length or a snapshot's files differ in date

        Formula:
            (d, x_1, \dots, x_n) \text{ for each date } d, \text{ with } \text{date}(f_k) = d \text{ for every variable } k

        See also:
            Source.stream, read_one
    """
    if len({len(f) for f in file_lists}) != 1:
        raise ValueError(f"file counts differ: {[len(f) for f in file_lists]}")
    for group in zip(*file_lists):
        dates = {file_date(f) for f in group}
        if len(dates) != 1:
            raise ValueError(f"dates differ: {group}")
        yield (dates.pop(), *(read_one(f, n, level, box) for f, n in zip(group, names)))


def read_stack(files, var, level=0, box=None):
    r"""
    read_stack: read files into one array, one snapshot per file.

    For records small enough to hold in memory, such as one level over a year.

        Parameters:
            files  (list): Files, in the order wanted
            var     (str): Common or NEMO name of the variable
            level   (int): One level, or None for all levels; ignored for 2D fields
            box   (tuple): (j0, j1, i0, i1) index box, or None for the whole grid

        Returns:
            stack (float): float32 (nt, ny, nx), or (nt, nz, ny, nx) with level=None on a 3D field;
                           land as NaN

        Formula:
            X[t] = \text{read\_one}(f_t)

        See also:
            read_one, Source.read
    """
    return np.stack([read_one(f, var, level, box) for f in files])


def read_mesh(paths, *names):
    r"""
    read_mesh: grid fields by name, from whichever mesh file holds them.

    A grid may be described by one mesh file or by several (horizontal, vertical, masks). Fields are looked
    up by name in each file in turn, with their singleton dimensions removed.

    Raises KeyError if a name is in no file.

        Parameters:
            paths  (list): Mesh files, searched in order
            *names  (str): Names of the fields, e.g. 'e1t', 'tmask'

        Returns:
            fields (array): One array per name, singleton dimensions removed; a single array if
                            one name is given. Raises KeyError if a name is in no file

        Formula:
            \text{field}(n) = \text{first file holding } n,\ \text{squeezed}

        See also:
            Source.mesh
    """
    out = []
    for name in names:
        for p in paths:
            with netCDF4.Dataset(p) as nc:
                if name in nc.variables:
                    out.append(np.squeeze(nc[name][:].data))
                    break
        else:
            raise KeyError(f"{name} not in {paths}")
    return out[0] if len(out) == 1 else out
