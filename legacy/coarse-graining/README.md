# Coarse graining

Standalone tool that filters high-resolution GYRE velocities (R27) and coarse-grains them (R3), with the fields for LU noise. Not built or tested here.

**Authors:** F.L. Tucciarone and L. Li. It is a complete reshaping by F.L. Tucciarone of L. Li's original code (2020, f77).

## `srcF90/`
From `crsgrain-base`, the later version: two filter widths (band), difference of Gaussians, time-mean removal, command-line options (`-w1`, `-w2`, `-units`, `-baroclinic`, `-timeavgrem`), SSH, and separate high- and low-resolution output. A run of this version finished normally over 690 snapshots.

To build, set the netCDF install folder: `make NCBASE=/path/to/netcdf/installs`. Only the personal paths in `make.macro` were changed.

## `srcPY/`
EOF, curl and divergence-free scripts from the old `develop` branch, commit `9cc82cc`. Their absolute paths were emptied, so set `base_dir` and the file names before running. `modredkit/pod.py` supersedes `snapshot_POD.py`.

The earlier Fortran (2022) is in tag `develop-2022`.

## Check before reuse
- `main.F90` builds input paths with `indir(1:25)`, so it only works with a directory prefix of exactly that length.
- The paths in `namelist_ioNetCDF.nml` are relative to an old folder layout.
