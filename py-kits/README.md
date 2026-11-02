# py-kits

Python tools for Location Uncertainty in NEMO: read ocean data, operate on ORCA grids, estimate the
statistics and structure of the stochastic noise, reduce records to a few modes, and diagnose the
results. Fortran is out of scope.

| Kit | What it is for |
|---|---|
| `datakit` | Reading the reanalyses (GLORYS04, GLORYS12, CMEMS), and fields computed from them, all in the same way |
| `orcakit` | Operations on NEMO's ORCA grids: C-grid stencils, the seam and fold, masks, filters, parallel tiles |
| `noisekit` | Statistics and structure of the stochastic noise; `noisekit/pipelines` holds the older noise methods |
| `modredkit` | POD, DMD, kernel DMD and multi-level DMD |
| `diagkit` | Metrics, flow quantities, statistics, spectra and maps |
| `gyrekit` | What is specific to the GYRE configurations |

Every routine takes and returns arrays, so it can be used at any step; `preprocessing/` and
`noisekit/pipelines` only read, call routines and write.

## Reference
`docs/py-kits.html` has a page for every kit, module and routine: what it does and why, its formula, an
example from the tests, where the library uses it, and the NEMO lines it follows, linked to the NEMO
forge. It is built from the docstrings and the executed test notebook:

    jupyter nbconvert --to notebook --execute --inplace tests/test_all_routines.ipynb
    python docs/make_html.py
    jupyter nbconvert --clear-output --inplace tests/test_all_routines.ipynb

`docs/figures/` holds the diagram of the noise pipeline and its script.

## Setup
Copy `paths.example.toml` to `paths.toml` and set the data locations; `nemo_src` lists NEMO source
trees, needed only to build the reference. Install with `pip install -e .`, or put `py-kits` on
`PYTHONPATH`. Parallel scripts need `mpi4py`, maps need `cartopy`.

## Tests
| File | Checks |
|---|---|
| `tests/test_all_routines.ipynb` | Every routine against analytic or hand results |
| `tests/test_kits.py`, `tests/test_datakit.py` | Synthetic checks, and reading the real datasets |
| `tests/datakit_*.ipynb`, `tests/modredkit_*.ipynb` | Worked examples on real data |
| `tests/bugfix_report.md` | The bugs fixed when the code was brought here, and their effect on real data |

    python -B tests/test_kits.py
    python -B tests/test_datakit.py
