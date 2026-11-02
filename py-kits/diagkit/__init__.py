"""
Diagnostics of model output and of data: metrics, flow quantities, statistics, spectra and maps.

Once a model has run, or a reduced model has been built, it has to be judged against a reference.
diagkit holds the measures used for that: error and pattern metrics and a relative entropy for
ensembles, strain, vorticity and the Okubo-Weiss parameter, kinetic energy and anomalies, spectra in
space and time, conservation checks, and the helpers for maps of the study regions.

Modules:
    metrics       weighted rms, rms error, pattern correlation and Gaussian relative entropy
    flow          strain, vorticity and the Okubo-Weiss parameter on a uniform C grid
    stats         anomalies, kinetic energy, horizontal means and the kinetic energy spectrum
    spectra       isotropic spectra in space and power spectra in time
    conservation  volume means of tracers and layer volumes, to check conservation
    fancy_plots   colormaps, the study regions and maps that locate them
"""
