"""
GLORYS04 (GLORYS2V4, 1/4 degree, 1021 x 1442, 75 levels), daily, 2015-2016.
sossheig stops at 2015-12-30. baroclinic=True reads the barotropic-removed u and v.

Variables: sossheig, vozocrtx, vomecrty
Dates are 'YYYY-MM-DD' or 'YYYY-MM', inclusive.
"""
import os

from orcakit import paths

from ._read import find_files, read_mesh, read_stack, stream_files
from ._source import Source

ROOT = paths.get("glorys04")
MESH = [os.path.join(ROOT, "static", "GLORYS2V4_mesh_mask.nc")]
SOURCE = Source("glorys04", os.path.join(ROOT, "daily"), MESH, "glorys04", 86400.0)


def files(var, start=None, end=None, baroclinic=False):
    r"""
    files: the files of one GLORYS04 variable between two dates.

    The baroclinic option reads an older folder of baroclinic velocities and is due to be removed.

        Parameters:
            var     (str): Common or NEMO name of the variable
            start   (str): First date, inclusive, or None
            end     (str): Last date, inclusive, or None
            baroclinic (bool): Read the old daily-baroclinic/ tree; to be removed

        Returns:
            files  (list): Paths sorted by date

        Formula:
            \text{find\_files}(\text{root}, v, \text{start}, \text{end})

        See also:
            datakit._source.Source.files
    """
    return find_files(os.path.join(ROOT, "daily-baroclinic" if baroclinic else "daily"), var, start, end)


def read(var, start=None, end=None, level=0, box=None, baroclinic=False):
    r"""
    read: read a GLORYS04 variable between two dates into one array.

        Parameters:
            var     (str): Common or NEMO name of the variable
            start   (str): First date, inclusive, or None
            end     (str): Last date, inclusive, or None
            level   (int): One level, or None for all levels
            box   (tuple): (j0, j1, i0, i1) index box, or None for the whole grid
            baroclinic (bool): Read the old daily-baroclinic/ tree; to be removed

        Returns:
            stack (float): (nt, ny, nx), or (nt, nz, ny, nx) with level=None; land as NaN

        Formula:
            X[t] = \text{read\_one}(f_t)

        See also:
            datakit._source.Source.read
    """
    return read_stack(files(var, start, end, baroclinic), var, level, box)


def stream(*variables, start=None, end=None, level=0, box=None, baroclinic=False):
    r"""
    stream: read several GLORYS04 variables one snapshot at a time.

        Parameters:
            *variables (str): Common or NEMO names
            start      (str): First date, inclusive, or None
            end        (str): Last date, inclusive, or None
            level      (int): One level, or None for all levels
            box      (tuple): (j0, j1, i0, i1) index box, or None for the whole grid
            baroclinic (bool): Read the old daily-baroclinic/ tree; to be removed

        Returns:
            generator: (date, field, ...) per snapshot

        Formula:
            (d, x_1, \dots, x_n) \text{ for each date } d

        See also:
            datakit._source.Source.stream
    """
    return stream_files([files(v, start, end, baroclinic) for v in variables], variables, level, box)


def mesh(*names):
    r"""
    mesh: GLORYS04 grid fields by name.

        Parameters:
            *names (str): Names of the fields

        Returns:
            fields (array): One array per name, or a single array for one name

        Formula:
            \text{read\_mesh}(\text{mesh files}, n_1, \dots)

        See also:
            datakit._source.Source.mesh
    """
    return read_mesh(MESH, *names)
