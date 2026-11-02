"""
Statistics and structure of the stochastic noise of Location Uncertainty, from ocean data.

In the Location Uncertainty (LU) framework a model is forced by a random velocity whose statistics
represent the scales it cannot resolve. noisekit estimates those statistics from high-resolution data:
running means, variances and covariances, their change over time, correlation lengths and times,
and random fields with a given spatial structure. The pipelines in noisekit.pipelines combine these
with the other kits to write the files NEMO reads.

Modules:
    welford      one-pass mean, variance and covariance, sample by sample, for records too long for memory
    running      the same statistics over a sliding window, to follow how variability changes in time
    correlation  correlation lengths and times by lagged correlation, and their e-folding lags
    structure    spatially correlated random fields, from an anisotropic diffusion (SPDE) operator
    noise        geostrophic velocity from sea level, in its earlier form
    pipelines    the noise methods, each a script from the data to NEMO's input file
"""
