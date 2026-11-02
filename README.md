# Transport under Location Uncertainty in NEMO

This branch is the latest working version. It is being rebuilt: the Python tools and preprocessing are in; the current Fortran and the install documents are still to be added.

## Layout
| Folder | Contents |
|---|---|
| `py-kits/` | Python tools: data reading, ORCA grid operations, noise statistics, model reduction, diagnostics. Its reference is `py-kits/docs/py-kits.html` |
| `preprocessing/` | Steps PP1 to PP7 that turn reanalysis data into the fields the noise is built from |
| `noise-models/` | The noise models under test, with the rules a model must follow; empty for now |
| `legacy/` | Routines no longer used, kept because they would be costly to rewrite |

## Branches
- `main`: latest working version.
- `develop`: new features under test, merged into `main` when they work.
- `JAMES-2024`, tag `frozen-JAMES2025`: frozen code and configurations of the JAMES paper, including `nemo-config`, `nemo-install`, `nemo-diagnostic` and `python-preprocessing`.

## Licence
Fortran derived from NEMO is under CeCILL 2.0, as NEMO is; everything else is under the MIT licence.
See `LICENSE.md`.
