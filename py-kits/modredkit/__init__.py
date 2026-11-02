"""
Model reduction: the decompositions that turn a long record of fields into a few patterns.

A record of ocean fields is huge, but most of its variability lives in a few spatial patterns. Proper
orthogonal decomposition (POD) finds the patterns that carry the most variance; dynamic mode
decomposition (DMD) finds patterns that each oscillate at one frequency and grow or decay at one
rate, and its kernel form extends it to nonlinear dynamics. These modes are used to build reduced
models and, in this library, to give the stochastic noise a structure in space and time.

Modules:
    pod    snapshot POD, weighted or not, for one field or several fields together
    dmd    exact DMD, kernel DMD, windowed multi-level DMD, amplitudes, eigenvalue filters and plots
    modes  working with snapshot matrices and modes: packing fields, selecting modes, rebuilding records

Conventions: snapshot matrices are (space, time), one snapshot per column, and the decompositions
work in float64. A rank given as an int is a number of modes; as a float it is a share of
information to keep (100.0 keeps all).
"""
