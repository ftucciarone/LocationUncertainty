"""
GLORYS12 (1/12 degree, 3059 x 4322), monthly. Heavy: read a box or one level.
sossheig and sosstsst cover 2015; 3D fields cover 1993-2002.

Variables: sossheig, sosstsst, vozocrtx, vomecrty, vovecrtz, votemper
Dates are 'YYYY-MM-DD' or 'YYYY-MM', inclusive.
"""
import os

from orcakit import paths

from ._read import find_files, read_mesh, read_stack, stream_files
from ._source import MONTH, Source

ROOT = paths.get("glorys12")
MESH = [os.path.join(ROOT, "static", f) for f in
        ("glorys12v1_mesh_hgr.nc", "glorys12v1_mesh_zgr.nc", "glorys12v1_mask.nc", "glorys12v1_ORCA12_bathymetry.nc")]
SOURCE = Source("glorys12", os.path.join(ROOT, "monthly"), MESH, "glorys12", MONTH)


def files(var, start=None, end=None):
    r"""
    files: the files of one GLORYS12 variable between two dates.

        Parameters:
            var     (str): Common or NEMO name of the variable
            start   (str): First date, inclusive, or None
            end     (str): Last date, inclusive, or None

        Returns:
            files  (list): Paths sorted by date

        Formula:
            \text{find\_files}(\text{root}, v, \text{start}, \text{end})

        See also:
            datakit._source.Source.files
    """
    return find_files(os.path.join(ROOT, "monthly"), var, start, end)


def read(var, start=None, end=None, level=0, box=None):
    r"""
    read: read a GLORYS12 variable between two dates into one array.

    GLORYS12 is heavy: read one level or a box rather than the whole grid where possible.

        Parameters:
            var     (str): Common or NEMO name of the variable
            start   (str): First date, inclusive, or None
            end     (str): Last date, inclusive, or None
            level   (int): One level, or None for all levels
            box   (tuple): (j0, j1, i0, i1) index box, or None for the whole grid

        Returns:
            stack (float): (nt, ny, nx), or (nt, nz, ny, nx) with level=None; land as NaN

        Formula:
            X[t] = \text{read\_one}(f_t)

        See also:
            datakit._source.Source.read
    """
    return read_stack(files(var, start, end), var, level, box)


def stream(*variables, start=None, end=None, level=0, box=None):
    r"""
    stream: read several GLORYS12 variables one snapshot at a time.

        Parameters:
            *variables (str): Common or NEMO names
            start      (str): First date, inclusive, or None
            end        (str): Last date, inclusive, or None
            level      (int): One level, or None for all levels
            box      (tuple): (j0, j1, i0, i1) index box, or None for the whole grid

        Returns:
            generator: (date, field, ...) per snapshot

        Formula:
            (d, x_1, \dots, x_n) \text{ for each date } d

        See also:
            datakit._source.Source.stream
    """
    return stream_files([files(v, start, end) for v in variables], variables, level, box)


def mesh(*names):
    r"""
    mesh: GLORYS12 grid fields by name.

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
