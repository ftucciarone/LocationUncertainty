"""Write fields in the common file convention: one variable and date per file, NEMO names, time in hours since 1950-01-01, NaN on land, and the history of how each file was made."""
import os
import subprocess
import sys
import time

import netCDF4
import numpy as np

from ._read import NAMES
from ._source import TIME_UNITS

ATTRS = {"vozocrtx": ("Zonal velocity", "m s-1"),
         "vomecrty": ("Meridional velocity", "m s-1"),
         "vovecrtz": ("Vertical velocity", "m s-1"),
         "sossheig": ("Sea surface height", "m"),
         "vobtcrtx": ("Barotropic zonal velocity", "m s-1"),
         "vobtcrty": ("Barotropic meridional velocity", "m s-1")}


def history(inputs):
    r"""
    history: the history a new file inherits from what it read.

    A file records every step that made it in its history attribute, newest first, so this collects the
    history of its inputs. A file written by write_field passes on its own; a raw file passes on its name
    and, if it has one, its own history.

    Several inputs whose histories may differ, such as u and v, are listed one after the other, each
    under its label and indented.

        Parameters:
            inputs (str): Path read, or {label: path} for several inputs whose histories may
                          differ; each is then listed under 'inputs: <label>', indented

        Returns:
            history (str): Lines, newest first

        Formula:
            H = \text{history}(\text{input}),\ \text{or}\ \bigoplus_\ell \big[\ \text{inputs: } \ell \ \oplus\ \text{indent}(H_\ell)\ \big]

        See also:
            provenance
    """
    def one(path):
        r"""
        one: the history one input passes on.

            Parameters:
                path (str): File read

            Returns:
                history (str): Its history attribute if a pipeline step wrote it; otherwise its name,
                               then its own history if any

            Formula:
                H = \text{history attribute, or ``source file: name'' and its own history}

            See also:
                history
        """
        with netCDF4.Dataset(path) as nc:
            own = nc.getncattr("history") if "history" in nc.ncattrs() else ""
            if "command" in nc.ncattrs():
                return own
        return "\n".join(line for line in (f"source file: {os.path.basename(path)}", own) if line)

    if isinstance(inputs, dict):
        return "\n".join(f"inputs: {label}\n" + "\n".join("  " + line for line in one(path).splitlines())
                         for label, path in inputs.items())
    return one(inputs)


def provenance(step=None, inputs=None, **extra):
    r"""
    provenance: the global attributes saying how a file was made.

    A result is only useful if it can be traced back: the commit and command that made it, when, from
    what. With a step, the history attribute gets one line, the date, the step and the command, above the
    history of the inputs, so a file carries its whole chain of steps down to the raw data.

        Parameters:
            step     (str): e.g. 'PP5 coarse_grain: bin onto domain_cfg.nc', or None
            inputs   (str): File read, or {label: path}; used with step
            **extra       : Further attributes, e.g. source='cmems'

        Returns:
            attrs (dict): git_commit, command, created; with a step also processing, history and
                          source_file; then extra

        Formula:
            \text{history} = \text{date}\ \text{step}\ |\ \text{command} \ \oplus\ \text{history}(\text{inputs})

        See also:
            history, write_field
    """
    here = os.path.dirname(os.path.abspath(__file__))
    try:
        commit = subprocess.run(["git", "-C", here, "describe", "--always", "--dirty"],
                                capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "unknown"
    created = time.strftime("%Y-%m-%dT%H:%M:%S")
    command = " ".join([os.path.basename(sys.argv[0]), *sys.argv[1:]])
    out = {"git_commit": commit,
           "command": " ".join(sys.argv),
           "created": created}
    if step:
        out["processing"] = step
        out["history"] = f"{created} {step} | {command}"
        if inputs:
            out["history"] += "\n" + history(inputs)
            paths = inputs.values() if isinstance(inputs, dict) else [inputs]
            out["source_file"] = ", ".join(os.path.basename(p) for p in paths)
    return {**out,
            **extra}


def write_field(root, var, date, field, lat, lon, time_hours, tag, depth=None, dtype=np.float32, attrs=None, extra=None, depth_name="deptht"):
    r"""
    write_field: write one snapshot in the common file convention.

    Every result written this way reads back like a raw dataset: root/<var>/<var>_<tag>_<date>.nc, NEMO's
    variable names, a scalar time in hours since 1950-01-01, NaN as fill value, compressed and chunked.
    Extra variables, such as a barotropic part, can share the file.

    The file is written under a temporary name and then renamed, so a stopped run never leaves a partial
    file behind.

        Parameters:
            root        (str): Folder of the dataset written
            var         (str): Common or NEMO name of the variable
            date        (str): Date in the file name
            field     (float): (ny, nx), or (nz, ny, nx) with depth; land as NaN
            lat       (float): (ny, nx) latitude of the field's points
            lon       (float): (ny, nx) longitude of the field's points
            time_hours(float): Time, in hours since 1950-01-01
            tag         (str): Tag in the file name
            depth     (float): (nz,) levels, for a 3D field
            dtype     (dtype): Stored type
            attrs      (dict): Global attributes, see provenance
            extra      (dict): {name: field} also stored in this file, e.g. the barotropic part
            depth_name  (str): 'deptht', or 'depthw' for fields at the top of the T-cells

        Returns:
            path (str): The file written

        Formula:
            \text{root}/v/v\_\text{tag}\_\text{date}.nc \leftarrow \{v:\ x,\ \text{extra},\ \varphi,\ \lambda,\ t,\ z\}

        See also:
            derived, provenance
    """
    var = NAMES.get(var, var)
    folder = os.path.join(root, var)
    os.makedirs(folder, exist_ok=True)
    path = os.path.join(folder, f"{var}_{tag}_{date}.nc")
    tmp = f"{path}.tmp{os.getpid()}"
    ny, nx = field.shape[-2:]
    with netCDF4.Dataset(tmp, "w") as nc:
        for name, value in (attrs or {}).items():
            nc.setncattr(name, value)
        dims, chunks = ("y", "x"), (min(256, ny), min(361, nx))
        nc.createDimension("y", ny)
        nc.createDimension("x", nx)
        if field.ndim == 3:
            nc.createDimension(depth_name, field.shape[0])
            dims, chunks = (depth_name, *dims), (min(19, field.shape[0]), *chunks)
            d = nc.createVariable(depth_name, "f4", (depth_name,))
            d[:] = depth
            d.units, d.long_name, d.standard_name = "m", "Nominal levels", "depth"
        for name, size in (("x", nx), ("y", ny)):
            nc.createVariable(name, "i4", (name,))[:] = np.arange(size)
        t = nc.createVariable("time_counter", "f8", ())
        t[...] = time_hours
        t.units, t.calendar, t.standard_name = TIME_UNITS, "gregorian", "time"
        for name, values, units, standard in (("nav_lat", lat, "degrees_north", "latitude"), ("nav_lon", lon, "degrees_east", "longitude")):
            c = nc.createVariable(name, "f4", ("y", "x"), zlib=True, complevel=4)
            c[:] = values
            c.units, c.standard_name = units, standard
        for name, values in ((var, field), *(extra or {}).items()):
            ndim = np.ndim(values)
            v = nc.createVariable(name, np.dtype(dtype), dims[-ndim:], zlib=True, complevel=4, shuffle=True, chunksizes=chunks[-ndim:], fill_value=np.nan)
            v[:] = np.asarray(values, dtype=dtype)
            v.coordinates = "time_counter nav_lat nav_lon"
            if name in ATTRS:
                v.long_name, v.units = ATTRS[name]
    os.replace(tmp, path)
    return path
