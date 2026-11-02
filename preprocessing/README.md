# Preprocessing

Steps that read a source and write a new one in the same format, numbered PP1 to PP7. Each script reads, calls the routines of `py-kits` and writes; the computation itself lives in the kits. Data locations come from `py-kits/paths.toml`.

| Step | Script | What it writes | Status |
|---|---|---|---|
| PP1 | `baroclinic.py` | Barotropic and baroclinic parts of u and v | Validated |
| PP2 | `vertical_velocity.py` | Vertical velocity by continuity | Validated |
| PP3 | `geostrophic.py` | Surface geostrophic velocity from SSH | Validated |
| PP4 | `spatial_filter.py` | Gaussian low-pass or high-pass in space | Validated |
| PP5 | `coarse_grain.py` | Fields averaged or sampled onto a target grid | Validated |
| PP6 | `time_filter.py` | Deseasonalised or Butterworth-filtered fields | **Untested** |
| PP7 | `running_statistics.py` | Running variance, std, mean or u-v covariance | Validated |

`check_baroclinic.py` and `check_derived.py` check the outputs of PP1 to PP3; `_common.py` holds the options, resuming and splitting of dates that the scripts share. "Validated" means the results were looked at and judged right, on the data listed under each step.

Each file records its chain of steps in the global attribute `history`, newest first: the step, what it did and its command, then the history of the file it read, down to the raw source file. `ncdump -h` shows it. Steps that read two inputs, such as the covariance, keep both histories.

## Baroclinic and barotropic split (PP1)

    export PYTHONPATH=/path/to/LocationUncertainty/py-kits
    mpirun -np 16 python baroclinic.py --source glorys12 --start 1993-01 --end 1993-12 --dst /path/to/out
    python check_baroclinic.py --source glorys12 --dst /path/to/out --start 1993-01 --end 1993-12

| Option | Meaning |
|---|---|
| `--source` | `glorys04`, `glorys12` or `cmems` |
| `--start`, `--end` | Dates, inclusive, `YYYY-MM-DD` or `YYYY-MM` |
| `--dst` | Writes `split/<var>/<var>_<tag>_<date>.nc` here |
| `--dtype` | `float32` (default) or `float64`, for both the computation and the file |
| `--overwrite` | Recompute existing files. Without it they are skipped, and a run stops if one has another precision or lacks a part |
| `--dry-run` | List what each rank would write |

Each file holds the baroclinic part under `vozocrtx` or `vomecrty`, and the barotropic part under NEMO's names `vobtcrtx` or `vobtcrty`.

Dates are split round-robin across ranks. Read the results back with `datakit.derived("/path/to/out/split", "glorys12")`: `read("u", level=None)` gives the baroclinic part and `read("ubt")` the barotropic part.

The thicknesses `e3u` and `e3v` come from the mesh, fixed in time, with no dependence on SSH.

Measured in float32:

| Source | Time per date | Peak memory per rank | Output per date |
|---|---|---|---|
| GLORYS12 | 110 s | 29 GB | 2.2 GB |
| CMEMS | 15 s | not measured | 360 MB |
| GLORYS04 | 18 s on 3 ranks | not measured | 370 MB |

Mind the memory on GLORYS12: keep the ranks per node below the node's memory divided by 29 GB.

## Vertical velocity (PP2)

    mpirun -np 16 python vertical_velocity.py --source cmems --start 2015-01 --end 2015-12 --dst /path/to/out
    python check_derived.py --step w --source cmems --dst /path/to/out --start 2015-01 --end 2015-12

Options as above. Writes `w/vovecrtz/vovecrtz_<tag>_<date>.nc` on `depthw`, from the full velocity, by continuity with `e3u` and `e3v` from the mesh. `w(k)` is at the top of T-cell `k`. Read back with `datakit.derived("/path/to/out/w", "cmems").read("w", level=None)`.

## Geostrophic velocity (PP3)

    mpirun -np 16 python geostrophic.py --source cmems --start 2015-01 --end 2015-12 --dst /path/to/out
    python check_derived.py --step geostrophic --source cmems --dst /path/to/out --start 2015-01 --end 2015-12

Options as above, plus `--min-lat` (default 5): NaN for `|lat|` below it, where geostrophy does not hold. Writes `geostrophic/<var>/<var>_<tag>_<date>.nc`, 2D, under `vozocrtx` and `vomecrty`.

Known limits: no ice-cover mask, so velocities under ice are kept; within about two cells of land the wet-only SSH average biases the velocity.

## Spatial filter (PP4)

    mpirun -np 16 python spatial_filter.py --source glorys04 --keep low --sigma 8 --dst /path/to/out

Options as above, plus `--keep low` or `high`, `--sigma` in grid points of the source, and `--vars` (default `u v`). Writes `lowpass/<var>/` or `highpass/<var>/`. Land is NaN, as are the few wet points with under 1% of the filter footprint valid. On GLORYS12 with `sigma` 8 a rank needs about 31 GB.

## Regridding (PP5)

    mpirun -np 16 python coarse_grain.py --source /path/to/out/lowpass --like glorys04 --dst /path/to/out

Options as above, plus `--target` (a `domain_cfg`, default eORCA1), `--method bin` (area-weighted average) or `nearest`, and `--vars`. Writes `coarse/<var>/` with the target's coordinates. The target's mask is used only at the surface, and the levels stay those of the source.

## Time filter (PP6)

**Untested.** It passes its synthetic tests but has never been checked on real data: check its output before using it.

    mpirun -np 16 python time_filter.py --source glorys04 --method deseason --dst /path/to/out
    mpirun -np 16 python time_filter.py --source cmems --method highpass --cutoff 365 --space low --sigma 8 --order time-space --dst /path/to/out

Options as above, plus:

| Option | Meaning |
|---|---|
| `--method` | `deseason`, or Butterworth `lowpass` or `highpass` |
| `--harmonics` | `deseason`: harmonics of the year, 2 by default for annual and semiannual |
| `--keep-mean` | `deseason`: add the time mean back, which is otherwise removed with the trend |
| `--cutoff`, `--order-butter` | Butterworth: cutoff period in days, and filter order |
| `--space`, `--sigma` | Also a Gaussian `low` or `high` in space, with a halo exchanged between tiles |
| `--order` | `space-time` or `time-space` |
| `--vars` | Default `u v` |

Writes `<name>/<var>/`, one file per date, under the same names, where `<name>` says what was applied, such as `deseason` or `time_highpass_lowpass`.

Runs parallel in space, because it needs the whole record of each point: the rank count must be a perfect square (1, 4, 9, 16, ...). Each rank holds its tile of the whole record, so the memory needed in total is about the size of the record. With `--space`, `4*sigma + 5` must fit in the smallest tile. A run is all or nothing: if every output exists it is skipped, otherwise everything is recomputed. A Butterworth filter is unreliable near the start and end of the record, and `deseason` needs a few years to separate the trend from the seasons.

## Running statistic (PP7)

    mpirun -np 16 python running_statistics.py --source glorys04 --window 21 --stride 7 --stat std --dst /path/to/out

Options as above, plus `--window`, `--stride`, in snapshots, `--stat var`, `std`, `mean` or `cov`, and `--vars`. Writes `running_<stat>/<var>/`, dated by the newest snapshot of the window. `--stat cov` takes two `--vars`, u and v by default, and writes `uv_tvcov` (or `--cov-name`) on T-points: U is first averaged in x and V in y onto T-points, with land faces as 0, and land T-points are NaN.

## Chaining

Any step can read the output of another: give its folder as `--source` and the raw source as `--like`. For example, low-pass then coarse-grain:

    mpirun -np 16 python spatial_filter.py --source glorys04 --keep low --sigma 8 --dst /path/to/out
    mpirun -np 16 python coarse_grain.py --source /path/to/out/lowpass --like glorys04 --dst /path/to/out

The computations are primitives in the kits (`orcakit.filters`, `orcakit.orca`, `orcakit.grid_c`, `noisekit.running`), so they can be called directly at any other step. 

Every output carries the latitude and longitude of the point its field lives on, from the mesh: T for `w` and the covariance, U and V for the velocity components, and the target's points for regridded fields.
