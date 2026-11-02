#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module stats

Time-mean and fluctuation diagnostics. Time is axis 0, y and x the last two.
Extracted from the DIAG_NEMO class of the earlier NEMO diagnostic code.
"""
import numpy as np
from orcakit.grid_c import ke_compute

from .spectra import init_2Dwvnb, iso_spec


def fluctuations(q, delay=0):
    r"""
    fluctuations: anomalies in time, removing the means of the spin-up and of the rest separately.

    A model's first months drift towards equilibrium; removing one mean over the whole run would mix that
    drift into the anomalies.

        Parameters:
            q   (float): (nt, ...) record
            delay (int): Snapshots of spin-up, 0 for none

        Returns:
            q_prime (float): (nt, ...) anomalies

        Formula:
            q' = q - \overline{q}^{[0, d)} \text{ on } t < d,\qquad q' = q - \overline{q}^{[d, n)} \text{ on } t \ge d

        See also:
            hovmoller
    """
    out = np.zeros_like(q)
    if delay:
        out[:delay] = q[:delay] - q[:delay].mean(axis=0, keepdims=True)
    out[delay:] = q[delay:] - q[delay:].mean(axis=0, keepdims=True)
    return out


def kinetic_energy(u, v):
    r"""
    kinetic_energy: kinetic energy on interior T-points, from velocities on the C grid.

        Parameters:
            u (float): (..., ny, nx) U-point velocity
            v (float): (..., ny, nx) V-point velocity

        Returns:
            ke (float): (..., ny, nx), 0 on the outer rows and columns

        Formula:
            K = \tfrac14\big(u_{j,i-1}^2 + u_{j,i}^2 + v_{j-1,i}^2 + v_{j,i}^2\big) = \tfrac12\big(\overline{u^2} + \overline{v^2}\big)

        See also:
            orcakit.grid_c.ke_compute, ke_spectrum
    """
    return 0.25 * ke_compute(u, v)


def hovmoller(q):
    r"""
    hovmoller: horizontal mean, leaving time and any other leading axes.

    A Hovmoller diagram of the domain mean shows how a quantity evolves, for example the spin-up of
    kinetic energy.

        Parameters:
            q (float): (..., ny, nx) field

        Returns:
            q_mean (float): (...) mean over the last two axes

        Formula:
            \overline{q}(t) = \frac{1}{n_x n_y}\sum_{j,i} q_{t,j,i}

        See also:
            fluctuations
    """
    return np.mean(q, axis=(-2, -1))


def ke_spectrum(u, v, dx, dy):
    r"""
    ke_spectrum: isotropic kinetic energy spectrum on a uniform grid.

    How kinetic energy spreads across scales: its slope separates regimes of turbulence, and its fall-off
    shows the effective resolution of a model.

        Parameters:
            u  (float): (..., ny, nx) velocity; nx and ny even
            v  (float): (..., ny, nx) velocity
            dx (float): Spacing in x
            dy (float): Spacing in y

        Returns:
            kr  (float): Isotropic wavenumbers
            spec (float): (..., len(kr)) spectrum

        Formula:
            E(k_r) = \text{iso\_spec}\Big(\tfrac12\big(|\hat u|^2 + |\hat v|^2\big)\Big)

        See also:
            iso_spec, kinetic_energy
    """
    ny, nx = u.shape[-2:]
    wv, kr, efc = init_2Dwvnb(nx, ny, dx, dy)
    fh = 0.5*(np.abs(np.fft.fft2(u, axes=(-2, -1)))**2 + np.abs(np.fft.fft2(v, axes=(-2, -1)))**2)
    return kr, iso_spec(fh, wv, kr, efc)
