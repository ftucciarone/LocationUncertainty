"""
The noise methods of the catalogue, each a script from the data to NEMO's noise input.

Each pipeline reads a dataset, computes one kind of noise statistic with the kits and writes it in
the format tlu_gss.F90 reads. They are kept as they were ported and are being reworked, one by one,
into thin wrappers of the kits' routines.

Modules:
    variance                 V1 to V3: constant mean, variance and covariance of u and v
    v4_bcbts                 V4: statistics of the barotropic and baroclinic parts
    runvar_dmd               D1 and D2: DMD of the running variance, exact or kernel
    w_modes                  W1: DMD modes of the vertical velocity
    build_noise_sigma        statistics of the small-scale baroclinic residual, by high-pass
    build_noise_lscale       correlation lengths of the same residual
    build_noise_BcResid_stats  statistics of the residual, by coarse-graining
    build_baroclinic_daily   the daily baroclinic velocities these read
"""
