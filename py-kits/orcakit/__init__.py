"""
Operations on NEMO's ORCA grids: the C-grid stencils, the seam and fold, masks, filters and parallel tiles.

The ORCA family of global grids is curvilinear and staggered (Arakawa C: tracers at cell centres,
velocities on the faces), periodic east-west with two repeated columns, and folded onto itself at the
north. Any operation that looks at neighbours, a gradient or a filter, must respect all of that, or it
gives wrong values along the seam and near the pole. orcakit gathers those operations once, written
to match NEMO's own discretisation, so the rest of the library never handles a boundary by hand.

Modules:
    orca        the seam and the north fold (halo), Gaussian smoothing across them, masks, partial-step
                thicknesses, the barotropic split, nearest-point search and averaging onto a coarser grid
    grid_c      C-grid operators: averages between points, differences, divergence, vorticity,
                vertical velocity by continuity and geostrophic velocity
    filters     land-aware low-pass and high-pass in space, on the whole grid or on a tile, and the
                deseasonalising and Butterworth filters in time
    mpi_orca    the domain decomposition into tiles, halo exchange between ranks, and gathering
    mpi_filter  a stand-alone parallel filter script built on the same ideas, with a check against
                a whole-domain reference
    netcdf      writing fields and modes to NetCDF files with the dimensions of a grid file
    paths       data locations, read from paths.toml

Every routine takes and returns plain arrays, so it can be used at any step of any computation.

Conventions: north stays at the high y index, as NEMO stores it; flipping the grid, as some older
code did, puts the fold on the wrong edge. The seam and fold rules (orca.halo, orca.FOLD) are those
of the fine grids with a T-point pivot and two overlap columns, GLORYS04, GLORYS12 and CMEMS; eORCA1's
domain file has an F-point pivot and no overlap columns, and must not go through them.
"""
