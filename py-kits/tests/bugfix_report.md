## Bugfix report (2026-10-31)

All 589 checks in `tests/test_all_routines.ipynb` pass (87 test groups, 252 routines, 9 command-line entry points left out). `tests/test_kits.py` and `tests/test_datakit.py` pass (28 tests).

Results made before these fixes, including `outputs/*.nc` and the original folders, still carry the bugs. New results are in `outputs/fixed/`. No file NEMO reads has been replaced.

### Effect on real data

Old kit and new kit against the existing reference files. V3 without low-pass is the fair comparison for `noise_stat.nc`, which has no low-pass. The new V3 differs through the sampling map only, the one fix on that path.

| Comparison | Variable | Correlation | Largest difference |
|---|---|---|---|
| Old V3 vs `noise_stat.nc` | all | 0.99999 or more | 0.3 to 0.4% of the reference maximum |
| New V3 vs `noise_stat.nc` | `umean` | 0.996 | 0.90 m/s |
| New V3 vs `noise_stat.nc` | `vmean` | 0.975 | 0.75 m/s |
| New V3 vs `noise_stat.nc` | `uvvar` | 0.950 | 0.65 m²/s² |
| Old V4 vs `noise_BcBtS_stats.nc` | checked: `vbt_mean`, `ubt_mean` | exact to 4 decimals | below 0.0001 m/s |
| New V4 vs `noise_BcBtS_stats.nc` | `ubc_mean` | 0.9987 | 0.36 m/s |
| New V4 vs `noise_BcBtS_stats.nc` | `vbt_mean` | 0.922 | 0.22 m/s |

These correlations and largest differences count points that changed between land and water as 0. On points wet in both files, new V3 against the reference correlates at 0.982 (`vmean`) to 0.998 (`umean`).

The large differences are expected where the sampled point moved in a fine, unsmoothed field, but they are large, and they are concentrated at few points. In V4 `vbt_mean` the median difference is 0, the 99th percentile 0.002 m/s and the largest values are in eORCA1 rows 329 and 330.

### Fixes

| Kit | Routine | Bug | Fix | Measured effect |
|---|---|---|---|---|
| orcakit | `halo`, `smooth_2D`, `mpi_orca`, `mpi_filter` | North fold used the wrong geometry for U, V and F points; corners and the south edge were not filled | `FOLD` rule per grid point, checked on the mesh, with a `point` argument throughout | Result with 1 and 4 ranks now agrees to 6e-8 m/s, float rounding. `build_noise_sigma --check` showed 8e-3 m/s before |
| orcakit | `find_nearest_LatLon`, `find_nearest_grid` | Distance taken in latitude and longitude, not on the sphere | KD-tree on the sphere; the two repeated fine rows left out | 24590 of 65087 points changed. Median distance 10.61 to 10.37 km; 12 to 40 km down to 3 to 5 km at the points checked near the fold |
| orcakit | `create_umask_from_tmask`, `create_vmask_from_tmask` | Not implemented | Implemented | Surface masks match NEMO's, 0 differences |
| orcakit | `masked_smooth`, `scale_split` | The mask was smoothed with the field's sign flip at the fold | Separate `smooth_mask`, no flip | Applies to v near the fold |
| orcakit | `grid_c.interp_1D`, `fdiff_1D`, `hDiv_Tgrid`, `vort_Vgrid` | Wrong exchange sign or grid point; column 0 of the divergence wrong; vorticity used the wrong metrics | `point` argument, west halo, F-point metrics | Analytic tests |
| noisekit | `uv_from_ssh_dA` | Geostrophic components wrong | `Ug=-(g/f)deta/dy`, `Vg=+(g/f)deta/dx` | Analytic test |
| noisekit | `variance`, `v4_bcbts`, `build_noise_BcResid_stats`, `runvar_dmd`, `build_noise_sigma`, `build_noise_lscale` | u and v means smoothed with no sign flip at the fold | Means flip at U or V; variances and masks do not; covariance at T (`variance`) or U (`build_noise_sigma`) | With the map: V4 `vbt_mean` at eORCA1 row 330, column 294, 0.255 to 0.031 m/s. Row 328 nearby is 0.038 m/s |
| gyrekit | `coarse_grain_P` | Wrong divisor | `sum(shape)+1` | Analytic test |
| diagkit | `init_3Dwvnb`, `iso_spec` | Wrong axis order and 3D factor; no 3D input | Order (nz, ny, nx), 3D factor, 3D input | Analytic test |
| modredkit | `RelativeInformationContent` | Rank off by one when `r` is an integer or None | No +1 | Analytic test |
| modredkit | `getConjugatePairsMask`, `DMD_dynamics_noPairs` | Wrong pair selection and weights; decay left out | Im >= 0; weight 2 for complex modes, 1 for real; decay included | Analytic test |
| modredkit | `RemGtT`, `SortByFrequency`, mrDMD | Wrong threshold and order | Keeps abs(lambda) <= threshold + 1e-3; sorts by abs(Im); mrDMD uses `RemGtT` | Analytic test |
| modredkit | `DMD_vecmodes`, `DMD_phyModes` | Unused classes | Deleted | None |

### Not fixed
- **Old `daily-baroclinic/` tree.** Its interior matches a new computation exactly on the day compared, but `u` in the last column and `v` in the top row are wrong, by up to 0.34 m/s, from the old halo. They are duplicate points, about 21000 per variable and day. It has not been regenerated.
- **`build_noise_BcResid_stats.py`** has no `main`, so importing it starts the full computation.

- **Wet points sampling land.** 1086 of 65087 wet eORCA1 points sample a fine-grid land point at the surface for T, 1453 for U and 1486 for V. They are stored as 0, which is also what the old kit and the reference do (1496 for V).
- **50 points beyond the fine grid.** They lie at 77.9 to 78.6°S, south of the end of GLORYS04 at 77.0°S, up to 171.8 km from their sample. No map can improve them.
- **Old V3 against the reference.** The 0.4% difference is unexplained.
- **`RelativeInformationContent`** still uses the sum of `S`, not `S**2`, and says so.
- **Fold rule on the mesh.** It holds to 0.35 km at the 99th percentile. The largest errors, 3.4 km (T) to 35.5 km (V), sit at the pivot and the poles.
- **`grid_c.ke_compute`** still returns 4 times the kinetic energy; use `diagkit.stats.kinetic_energy`.
