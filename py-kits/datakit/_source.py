"""One description of a dataset, raw or derived, so that all sources are read the same way."""
import importlib
import os

import netCDF4
import numpy as np

from orcakit import orca

from ._read import NAMES, file_date, find_files, read_mesh, read_one, read_stack, stream_files

MONTH = 31 * 86400.0
TIME_UNITS = "hours since 1950-01-01"


class Source:
    """
    Fields under root/<var>/<var>_*.nc, with the grid of the mesh files.

        name, tag   : dataset name, and the tag in the names of files written from it
        dt          : time between snapshots, in seconds
        boundary    : 't-fold' is the only north boundary orca.halo handles
    """

    def __init__(self, name, root, mesh, tag, dt, boundary="t-fold"):
        r"""
        __init__: describe a dataset, raw or derived.

        A Source holds what is needed to read any dataset the same way: where its files are, which mesh
        describes its grid, how often it is sampled, and the tag of the files written from it.

            Parameters:
                name      (str): Dataset name, used in messages and provenance
                root      (str): Folder holding <var>/<var>_*.nc
                mesh     (list): Mesh files of the grid
                tag       (str): Tag in the names of the files written from this dataset
                dt      (float): Time between snapshots, in seconds
                boundary  (str): North boundary; 't-fold' is the only one orca.halo handles

            Returns:
                None

            Formula:
                \text{Source} = (\text{name},\ \text{root},\ \text{mesh},\ \text{tag},\ \Delta t,\ \text{boundary})

            See also:
                source, derived
        """
        self.name, self.root, self.mesh_files, self.tag, self.dt, self.boundary = name, root, mesh, tag, dt, boundary

    def files(self, var, start=None, end=None):
        r"""
        files: the files of one variable, between two dates.

            Parameters:
                var     (str): Common or NEMO name of the variable
                start   (str): First date, inclusive, or None
                end     (str): Last date, inclusive, or None

            Returns:
                files  (list): Paths sorted by date

            Formula:
                \text{find\_files}(\text{root}, v, \text{start}, \text{end})

            See also:
                find_files
        """
        return find_files(self.root, var, start, end)

    def read(self, var, start=None, end=None, level=0, box=None):
        r"""
        read: read a variable between two dates into one array.

        For records that fit in memory. Use stream for longer ones.

            Parameters:
                var     (str): Common or NEMO name of the variable
                start   (str): First date, inclusive, or None
                end     (str): Last date, inclusive, or None
                level   (int): One level, or None for all levels
                box   (tuple): (j0, j1, i0, i1) index box, or None for the whole grid

            Returns:
                stack (float): (nt, ny, nx), or (nt, nz, ny, nx) with level=None; land as NaN

            Formula:
                X[t] = \text{read\_one}(f_t),\quad t = \text{start}, \dots, \text{end}

            See also:
                Source.stream, read_stack
        """
        return read_stack(self.files(var, start, end), var, level, box)

    def stream(self, *variables, start=None, end=None, level=0, box=None):
        r"""
        stream: read several variables one snapshot at a time.

        Keeps one date in memory, so any length of record can be processed, with the dates of all variables
        checked to match.

            Parameters:
                *variables (str): Common or NEMO names
                start      (str): First date, inclusive, or None
                end        (str): Last date, inclusive, or None
                level      (int): One level, or None for all levels
                box      (tuple): (j0, j1, i0, i1) index box, or None for the whole grid

            Returns:
                generator: (date, field, ...) per snapshot, see stream_files

            Formula:
                (d, x_1, \dots, x_n) \text{ for each date } d

            See also:
                stream_files, Source.groups
        """
        return stream_files([self.files(v, start, end) for v in variables], variables, level, box)

    def groups(self, *variables, start=None, end=None):
        r"""
        groups: the dates on which all the variables have a file, with those files.

        Processing that needs several variables per date, such as u and v, can only use dates where all of
        them exist. This lists those dates with their files, and the dates where some are missing, so a run
        can report them.

            Parameters:
                *variables (str): Common or NEMO names
                start      (str): First date, inclusive, or None
                end        (str): Last date, inclusive, or None

            Returns:
                groups  (list): Sorted (date, path, path, ...), one path per variable
                missing (list): Sorted dates on which only some of the variables have a file

            Formula:
                D = \bigcap_k \text{dates}(v_k),\qquad \text{missing} = \bigcup_k \text{dates}(v_k) \setminus D

            See also:
                Source.files, Source.stream
        """
        by_date = [{file_date(f): f for f in self.files(v, start, end)} for v in variables]
        common = set.intersection(*(set(d) for d in by_date))
        missing = sorted(set.union(*(set(d) for d in by_date)) - common)
        return [(date, *(d[date] for d in by_date)) for date in sorted(common)], missing

    def mesh(self, *names):
        r"""
        mesh: grid fields by name, from the dataset's mesh files.

            Parameters:
                *names (str): Names of the fields, e.g. 'e1t', 'tmask'

            Returns:
                fields (array): One array per name, or a single array for one name

            Formula:
                \text{read\_mesh}(\text{mesh files}, n_1, \dots)

            See also:
                read_mesh, Source.grid
        """
        return read_mesh(self.mesh_files, *names)

    def point_mask(self, point):
        r"""
        point_mask: the wet mask of one kind of C-grid point.

        T, U, V and F points have different masks: a U point is wet only when both its T-cells are. The
        mask must match the point a field lives on.

            Parameters:
                point (str): 'T', 'U', 'V' or 'F'

            Returns:
                mask (bool): (nz, ny, nx), True on water

            Formula:
                m_p = p\text{mask from the mesh},\ p \in \{T, U, V, F\}

            See also:
                Source.point_coords, orcakit.orca.create_umask_from_tmask
        """
        return self.mesh(f"{point.lower()}mask").astype(bool)

    def grid(self):
        r"""
        grid: the horizontal grid: T-point positions and the cell widths at every point type.

        Differences, divergences and areas on the C grid all need the widths e1 and e2 of the cells at the
        right points. They are returned together, with the type of north boundary.

            Parameters:
                None

            Returns:
                grid (dict): 'lat' and 'lon' of T-points, 'boundary', and the scale factors e1t, e2t,
                             e1u, e2u, e1v, e2v, e1f, e2f in metres, each (ny, nx)

            Formula:
                \{\varphi_T,\ \lambda_T,\ e_{1p},\ e_{2p} : p \in \{t, u, v, f\}\}

            See also:
                Source.vertical, Source.point_coords
        """
        names = ("e1t", "e2t", "e1u", "e2u", "e1v", "e2v", "e1f", "e2f")
        lat, lon, *metrics = self.mesh("gphit", "glamt", *names)
        return {"lat": lat,
                "lon": lon,
                "boundary": self.boundary,
                **dict(zip(names, metrics))}

    def vertical(self, dtype=np.float32):
        r"""
        vertical: cell thicknesses and wet masks at U and V points.

        Depth means and transports need the thickness of every cell, including the partial cells at the sea
        floor. These are the mesh's thicknesses, fixed in time, so they do not follow the free surface.

        Raises NotImplementedError for a north boundary other than a T-point fold.

            Parameters:
                dtype (dtype): Type of the thicknesses

            Returns:
                e3u   (float): (nz, ny, nx) U-cell thicknesses, in metres, 0 on land
                e3v   (float): (nz, ny, nx) V-cell thicknesses, in metres, 0 on land
                umask  (bool): (nz, ny, nx) U-point wet mask
                vmask  (bool): (nz, ny, nx) V-point wet mask

            Formula:
                e_{3u},\ e_{3v} = \text{orca.e3uv\_from\_partial\_steps}(e_{3t}^0,\ e_{3t}^{ps},\ t, u, v)

            See also:
                orcakit.orca.e3uv_from_partial_steps, Source.point_mask
        """
        if self.boundary != "t-fold":
            raise NotImplementedError(f"orca.halo does not handle boundary {self.boundary}")
        e3t_0, e3t_ps, tmask, umask, vmask = self.mesh("e3t_0", "e3t_ps", "tmask", "umask", "vmask")
        tmask, umask, vmask = tmask.astype(bool), umask.astype(bool), vmask.astype(bool)
        e3u, e3v = orca.e3uv_from_partial_steps(e3t_0, e3t_ps, tmask, umask, vmask)
        return e3u.astype(dtype), e3v.astype(dtype), umask, vmask

    def point_coords(self, point):
        r"""
        point_coords: latitude and longitude of one kind of C-grid point, from the mesh.

        U and V points lie half a cell from the T points. Files often store the coordinates of one point type
        only; these are the right ones for a field at any point.

            Parameters:
                point (str): 'T', 'U', 'V' or 'F'

            Returns:
                lat (float): (ny, nx) latitude, in degrees
                lon (float): (ny, nx) longitude, in degrees

            Formula:
                (\varphi_p,\ \lambda_p) = (\text{gphi}p,\ \text{glam}p)

            See also:
                Source.coords, Source.point_mask
        """
        return self.mesh(f"gphi{point.lower()}", f"glam{point.lower()}")

    def coords(self, path):
        r"""
        coords: nav_lat and nav_lon as stored in one file.

        They may belong to another point than the file's field. point_coords gives those of a chosen point.

            Parameters:
                path (str): NetCDF file

            Returns:
                lat (float): (ny, nx) nav_lat
                lon (float): (ny, nx) nav_lon

            Formula:
                (\varphi, \lambda) = (\text{nav\_lat},\ \text{nav\_lon})

            See also:
                Source.point_coords
        """
        with netCDF4.Dataset(path) as nc:
            nc.set_auto_mask(False)
            return nc["nav_lat"][:], nc["nav_lon"][:]

    def depth(self, path):
        r"""
        depth: the depth levels of one file.

            Parameters:
                path (str): NetCDF file

            Returns:
                depth (float): (nz,) values of its dep* variable, or None for a 2D field

            Formula:
                z = \text{dep*}[:]

            See also:
                Source.snapshot
        """
        with netCDF4.Dataset(path) as nc:
            name = next((n for n in nc.variables if n.startswith("dep")), None)
            return nc[name][:].data if name else None

    def time_hours(self, path):
        r"""
        time_hours: the time of one file, in hours since 1950-01-01.

        Datasets count time from different origins and in different units. Converting all of them to one
        axis lets files from different sources be compared and written alike.

            Parameters:
                path (str): NetCDF file

            Returns:
                hours (float): Time of its first snapshot

            Formula:
                t = \text{date2num}(\text{num2date}(t_{\text{file}}),\ \text{hours since 1950-01-01})

            See also:
                write_field
        """
        with netCDF4.Dataset(path) as nc:
            t = nc["time_counter"]
            date = netCDF4.num2date(float(np.ravel(t[:])[0]), t.units, t.calendar)
            return float(netCDF4.date2num(date, TIME_UNITS, t.calendar))

    def snapshot(self, path, var, dtype=np.float32):
        r"""
        snapshot: all levels of one variable from one file.

            Parameters:
                path    (str): NetCDF file
                var     (str): Common or NEMO name of the variable
                dtype (dtype): Type of the result

            Returns:
                field (float): (ny, nx) or (nz, ny, nx), land as NaN

            Formula:
                x = \text{read\_one}(\text{path}, v)

            See also:
                read_one
        """
        return read_one(path, NAMES.get(var, var), None, None, dtype)


def source(name):
    r"""
    source: open a raw dataset by name.

        Parameters:
            name (str): 'glorys04', 'glorys12' or 'cmems'

        Returns:
            source (Source): The dataset, with its folder, mesh, tag and time step

        Formula:
            \text{datakit.}\langle\text{name}\rangle\text{.SOURCE}

        See also:
            derived, Source.__init__
    """
    return importlib.import_module(f"datakit.{name}").SOURCE


def derived(root, like, name=None):
    r"""
    derived: open data written with write_field, on the grid of a raw dataset.

    Fields computed from a dataset are on its grid and at its time step, so they read like it. This
    opens such a folder as a Source of its own, with the mesh, tag and time step of the original.

        Parameters:
            root  (str): Folder written by the step, holding <var>/<var>_*.nc
            like  (str): Raw dataset, or its Source, whose mesh, tag and time step the data share
            name  (str): Dataset name; the folder's name by default

        Returns:
            source (Source): The derived dataset

        Formula:
            \text{Source}(\text{name},\ \text{root},\ \text{like.mesh},\ \text{like.tag},\ \text{like.}\Delta t)

        See also:
            source, write_field
    """
    like = source(like) if isinstance(like, str) else like
    return Source(name or os.path.basename(os.path.normpath(root)), root, like.mesh_files, like.tag, like.dt, like.boundary)
