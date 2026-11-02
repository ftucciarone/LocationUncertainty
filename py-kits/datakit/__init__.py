"""
Read the ocean reanalyses, and the fields computed from them, all in the same way.

The datasets this library works from (GLORYS04 daily, GLORYS12 and the CMEMS ORAS5 reanalysis
monthly) differ in file layout, variable names, time units and grids. datakit hides those
differences: every dataset, raw or derived, is a Source with the same functions, returns arrays of
the same shape with land as NaN, and gives its grid, masks and cell thicknesses on request. Results
are written back in one file convention, with their history, so they read like a raw dataset and can
feed the next step.

Modules:
    _source   Source, the description of a dataset; source() opens a raw one, derived() a computed one
    _read     finding files by variable and date, and reading snapshots, records and mesh fields
    _write    write_field, the common file convention, and provenance and history attributes
    glorys04  GLORYS2V4, 1/4 degree, daily, 2015-2016
    glorys12  GLORYS12V1, 1/12 degree, monthly
    cmems     ORAS5 from CMEMS, 1/4 degree, monthly

Data locations come from paths.toml (orcakit.paths), so the same code runs on any machine.

Conventions: arrays are float32, land is NaN, and north is at the high y index, with no flip;
only plotting may flip. Time is axis 0, then depth, then y and x. Fields can be asked for by
common name (u, v, w, ssh, sst, temp, ubt, vbt) or by NEMO's name.

Example:
    from datakit import source, derived, clean

    src = source("glorys04")
    ssh = src.read("ssh", "2015-01", "2015-03")                     # (nt, ny, nx)
    u3d = src.read("u", "2015-01-01", "2015-01-01", level=None)     # (1, nz, ny, nx)
    box = source("glorys12").read("ssh", box=(1500, 1800, 1000, 1400))
    e1t, e2t = src.mesh("e1t", "e2t")
    e3u, e3v, umask, vmask = src.vertical()

    for date, u, v in src.stream("u", "v", start="2015-01", end="2015-02"):
        u = clean(u)                                                # one snapshot at a time

    groups, missing = src.groups("u", "v", start="2015")            # (date, u file, v file)
    split = derived("/path/to/out/split", src)                      # data written with write_field
"""
from ._read import BAROTROPIC, FIELDS, NAMES, clean, read_one  # noqa: E402,F401
from ._source import Source, derived, source  # noqa: E402,F401
from ._write import history, provenance, write_field  # noqa: E402,F401
