#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module netcdf

Generic NetCDF file utilities not tied to any particular grid or variable
naming convention.
"""

import os.path

import netCDF4
import numpy as np


def nc_copy(src_file, dst_file, atts=False, dims=False, vars=False, exclude=["TO_REMOVE"]):
    r"""
    nc_copy: start a new file from parts of another.

    Copies the global attributes, dimensions or variables of a template, typically a domain file, so that
    new fields can be added with the same dimensions.

    The new file is NETCDF3_CLASSIC and replaces any file of that name.

        Parameters:
            src_file (str): File to copy from
            dst_file (str): File to create
            atts    (bool): Copy the global attributes
            dims    (bool): Copy the dimensions
            vars    (bool): Copy the variables and their data
            exclude (list): Variables not to copy

        Returns:
            dst (Dataset): The new file, open for writing; the caller closes it

        Formula:
            \text{dst} \leftarrow \{\text{attributes},\ \text{dimensions},\ \text{variables} \setminus \text{exclude}\} \text{ of src, as selected}

        See also:
            write_fields
    """
    if os.path.exists(dst_file):
        os.remove(dst_file)

    src = netCDF4.Dataset(src_file)
    dst = netCDF4.Dataset(dst_file, "w", format="NETCDF3_CLASSIC")

    # copy attributes
    if atts:
        for name in src.ncattrs():
            dst.setncattr(name, src.getncattr(name))

    # copy dimensions
    if dims:
        for name, dimension in src.dimensions.items():
            dst.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))

    # copy all file data except for the excluded
    if vars:
        for name, variable in src.variables.items():
            if name not in exclude:
                dst.createVariable(name, variable.datatype, variable.dimensions)
                dst.variables[name][:] = src.variables[name][:]

    return dst


def write_fields(grid_file, dst_file, fields, dims=("z", "y", "x"), fill_value=0.0):
    r"""
    write_fields: write arrays to a new file with the dimensions of a grid file.

    Each array takes the last of the given dimension names that match its rank, so 3D and 2D fields can
    share one file.

    Masked and NaN values are stored as the fill value.

        Parameters:
            grid_file   (str): File whose dimensions are copied
            dst_file    (str): File to create
            fields     (dict): {name: array}; each array takes the last array.ndim names of dims
            dims      (tuple): Dimension names, outermost first
            fill_value (float): Value stored for masked and NaN points

        Returns:
            None

        Formula:
            \text{dst}[\text{name}] = \text{fields}[\text{name}],\ \text{masked and NaN} \to \text{fill\_value}

        See also:
            nc_copy, write_modes
    """
    dst = nc_copy(grid_file, dst_file, dims=True)
    try:
        for name, data in fields.items():
            data = np.ma.filled(np.ma.masked_invalid(data), fill_value)
            var = dst.createVariable(name, "f4", dims[-data.ndim:], fill_value=fill_value)
            var[...] = data
    finally:
        dst.close()


def write_modes(dst_file, means, modes, omega):
    r"""
    write_modes: write DMD modes in the format of NEMO's stochastic noise input (noise_DMDstatistics.nc).

    The mean fields, the real and imaginary parts of each mode, and their angular frequencies, with the
    variable names the model reads.

        Parameters:
            dst_file (str): File to create
            means   (dict): {name: (z, y, x)} mean fields, e.g. {'uavg': ...}
            modes   (dict): {prefix: complex (modes, z, y, x)}, written as <prefix>_real_r and
                            <prefix>_imag_r
            omega  (float): (modes,) angular frequencies, in rad/s

        Returns:
            None

        Formula:
            \text{prefix\_real\_r} = \mathrm{Re}\,\phi,\qquad \text{prefix\_imag\_r} = \mathrm{Im}\,\phi,\qquad \text{omega} = \omega

        See also:
            write_fields
    """
    shape = next(iter(means.values())).shape
    with netCDF4.Dataset(dst_file, "w") as nc:
        nc.createDimension("modes", len(omega))
        for name, size in zip(("z", "y", "x"), shape):
            nc.createDimension(name, size)
        nc.createVariable("modes", "i4", ("modes",))[:] = np.arange(len(omega))
        for name, size in zip(("z", "y", "x"), shape):
            nc.createVariable(name, "f8", (name,))[:] = np.arange(size)
        nc.createVariable("omega", "f8", ("modes",))[:] = omega
        for name, data in means.items():
            nc.createVariable(name, "f8", ("z", "y", "x"))[...] = data
        for prefix, data in modes.items():
            nc.createVariable(f"{prefix}_real_r", "f8", ("modes", "z", "y", "x"))[...] = data.real
            nc.createVariable(f"{prefix}_imag_r", "f8", ("modes", "z", "y", "x"))[...] = data.imag
