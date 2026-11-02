"""
What is specific to the GYRE configurations (R3, R9, R27): runs, masks and coarse-graining.

GYRE is NEMO's idealised double-gyre basin, used to test methods at several resolutions before the
global ocean. gyrekit holds only what depends on it: where the runs are and how to read them, its
masks and the passage between resolutions, and a standard summary of a run. Everything generic comes
from the other kits.

Modules:
    runs   the experiments of each resolution and their files, from config/paths_<R>.yaml
    grids  the jet mask, matching grids between resolutions and coarse-graining of pressure fields
    diag   the standard velocity diagnostics of a run
"""
