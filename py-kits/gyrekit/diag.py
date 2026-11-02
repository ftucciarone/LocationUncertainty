#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module diag

The velocity part of the old R3/R9 diagnostic, as one function over the other
kits. Replaces the DIAG_NEMO class, which computed and wrote everything on creation.
"""
import numpy as np

from diagkit.spectra import temporal_psd
from diagkit.stats import fluctuations, hovmoller, ke_spectrum, kinetic_energy


def velocity_summary(u, v, dx, dy, dt, delay=0, nb=13):
    r"""
    velocity_summary: standard diagnostics of GYRE velocities.

    One call gives what a GYRE run is usually judged by: mean flow, mean and eddy kinetic energy, their
    evolution, and their spectra in space and time. The spin-up is left out of means and spectra, and the
    points near the boundaries out of the time spectrum.

        Parameters:
            u   (float): (time, z, y, x) or (time, y, x) velocity
            v   (float): Same shape
            dx  (float): Spacing in x
            dy  (float): Spacing in y
            dt  (float): Time between snapshots
            delay (int): Snapshots of spin-up left out of means and spectra
            nb    (int): Points near the boundaries left out of the time spectrum

        Returns:
            summary (dict): umean, vmean, mke, eke, hov_ke, hov_eke, k, spec_ke, freq, psd_eke

        Formula:
            \text{MKE} = K(\bar u, \bar v),\qquad \text{EKE} = \overline{K(u', v')},\qquad \text{spectra of } K \text{ and } \tfrac12(P_u + P_v)

        See also:
            diagkit.stats.kinetic_energy, diagkit.stats.ke_spectrum, diagkit.spectra.temporal_psd
    """
    uf, vf = fluctuations(u, delay), fluctuations(v, delay)
    umean, vmean = u[delay:].mean(axis=0), v[delay:].mean(axis=0)
    weights = np.zeros(u.shape[-2:])
    weights[nb:-nb, nb:-nb] = 1. / np.prod(weights[nb:-nb, nb:-nb].shape)
    k, spec_ke = ke_spectrum(u[delay:], v[delay:], dx, dy)
    freq, psd_u = temporal_psd(u[delay:], 1./dt, weights)
    _, psd_v = temporal_psd(v[delay:], 1./dt, weights)
    return {"umean": umean,
            "vmean": vmean,
            "mke": kinetic_energy(umean, vmean),
            "eke": kinetic_energy(uf[delay:], vf[delay:]).mean(axis=0),
            "hov_ke": hovmoller(kinetic_energy(u, v)),
            "hov_eke": hovmoller(kinetic_energy(uf, vf)),
            "k": k,
            "spec_ke": spec_ke.mean(axis=0),
            "freq": freq,
            "psd_eke": 0.5*(psd_u + psd_v)}
