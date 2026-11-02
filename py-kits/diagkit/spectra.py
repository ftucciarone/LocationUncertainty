#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module spectra

Spatial and temporal power spectral density (PSD) tools.

Spatial: isotropic spectrum of a 2D/3D field, following the procedure in
Durran et al. (2017), "Practical considerations for computing dimensional
spectra from gridded data", http://doi.org/10.1175/MWR-D-17-0056.1

Temporal: periodogram-based PSD along the leading (time) axis of a field,
optionally reduced over the trailing spatial axes by a weight array.

Author: Long Li, long.li@inria.fr
        Francesco Tucciarone, francesco.tucciarone@inria.fr
"""

import numpy as np
from scipy import signal


def init_2Dwvnb(nx, ny, dx, dy):
    r"""
    init_2Dwvnb: the wavenumbers of a 2D discrete Fourier transform, and of its isotropic spectrum.

    An isotropic spectrum sums the energy of all Fourier coefficients with the same wavenumber modulus.
    This gives that modulus for every coefficient, the rings to sum over and the scaling of the result.

        Parameters:
            nx (int): Points in x, even
            ny (int): Points in y, even
            dx (float): Spacing in x
            dy (float): Spacing in y

        Returns:
            wv  (float): (nx*ny,) wavenumber modulus of each Fourier coefficient, in the order of
                         a flattened (ny, nx) array
            kr  (float): Isotropic wavenumbers, up to min(nx, ny)/2 steps of max(dkx, dky)
            efac (float): Scaling of the spectrum for the DFT and the integration

        Formula:
            k_x = \frac{2\pi}{n_x \Delta x}\,\{0, \dots, \tfrac{n_x}{2}-1, -\tfrac{n_x}{2}, \dots, -1\},\qquad |k| = \sqrt{k_x^2 + k_y^2},\qquad k_r = \Delta k_r\,\{1, \dots, \tfrac{\min(n_x, n_y)}{2}\}

        See also:
            iso_spec, init_3Dwvnb
    """

    assert nx % 2 == 0, 'Wavenumbers length in x-axis is not even'
    assert ny % 2 == 0, 'Wavenumbers length in y-axis is not even'

    # Create horizontal wavenumbers
    dkx = 2. * np.pi / (nx * dx)
    dky = 2. * np.pi / (ny * dy)
    kx = dkx * np.append(np.arange(0., nx / 2), np.arange(-nx / 2, 0.))
    ky = dky * np.append(np.arange(0., ny / 2), np.arange(-ny / 2, 0.))
    kxx, kyy = np.meshgrid(kx.astype('float64'), ky.astype('float64'))
    wv = np.sqrt(kxx**2 + kyy**2).flatten()

    # Create isotropic wavenumbers
    nmax = np.ceil(np.sqrt(2.) * np.maximum(nx, ny) / 2.)
    nmax = np.minimum(nx, ny) / 2.
    dkr = np.maximum(dkx, dky)
    kr = dkr * np.arange(1., nmax + 1, dtype=np.float64)

    # Scaling factor for spectrum due to DFT and integration
    efac = dx * dy * np.minimum(dkx, dky) / (4. * np.pi**2 * nx * ny)
    return wv, kr, efac


def init_3Dwvnb(nx, ny, nz, dx, dy, dz):
    r"""
    init_3Dwvnb: the wavenumbers of a 3D discrete Fourier transform, and of its isotropic spectrum.

        Parameters:
            nx (int): Points in x, even
            ny (int): Points in y, even
            nz (int): Points in z, even
            dx (float): Spacing in x
            dy (float): Spacing in y
            dz (float): Spacing in z

        Returns:
            wv  (float): (nz*ny*nx,) wavenumber modulus of each Fourier coefficient, in the order
                         of a flattened (nz, ny, nx) array
            kr  (float): Isotropic wavenumbers
            efac (float): Scaling of the spectrum for the DFT and the integration

        Formula:
            |k| = \sqrt{k_x^2 + k_y^2 + k_z^2},\ \text{in the order of a flattened } (n_z, n_y, n_x) \text{ array}

        See also:
            iso_spec, init_2Dwvnb
    """

    assert nx % 2 == 0, 'Wavenumbers length in x-axis is not even'
    assert ny % 2 == 0, 'Wavenumbers length in y-axis is not even'
    assert nz % 2 == 0, 'Wavenumbers length in z-axis is not even'

    # Create horizontal wavenumbers
    dkx = 2. * np.pi / (nx * dx)
    dky = 2. * np.pi / (ny * dy)
    dkz = 2. * np.pi / (nz * dz)

    kx = dkx * np.append(np.arange(0., nx / 2), np.arange(-nx / 2, 0.))
    ky = dky * np.append(np.arange(0., ny / 2), np.arange(-ny / 2, 0.))
    kz = dkz * np.append(np.arange(0., nz / 2), np.arange(-nz / 2, 0.))

    # Flattened in the axis order of a (nz, ny, nx) field
    kzz, kyy, kxx = np.meshgrid(kz.astype('float64'),
                                ky.astype('float64'),
                                kx.astype('float64'), indexing='ij')

    wv = np.sqrt(kxx**2 + kyy**2 + kzz**2).flatten()

    # Create isotropic wavenumbers
    nmax = np.ceil(np.sqrt(2.) * max(nx, ny, nz) / 2.)
    dkr = np.maximum(dkx, dky)
    kr = dkr * np.arange(1., nmax + 1, dtype=np.float64)

    # Scaling factor for spectrum due to DFT and integration
    efac = dx * dy * dz * np.minimum(np.minimum(dkx, dky), dkz)**2 / (8. * np.pi**3 * nx * ny * nz)
    return wv, kr, efac


def iso_spec(fh, wv, kr, efc):
    r"""
    iso_spec: the isotropic spectrum, summing Fourier coefficients over rings of wavenumber.

    The trailing two or three axes are spatial, three when their size matches the wavenumbers.

        Parameters:
            fh  (float): (..., ny, nx) or (..., nz, ny, nx) spectral density, e.g. |FFT|^2
            wv  (float): From init_2Dwvnb or init_3Dwvnb
            kr  (float): Isotropic wavenumbers, same source
            efc (float): Scaling factor, same source

        Returns:
            psdr (float): (..., len(kr)) isotropic spectrum

        Formula:
            E(k_r) = c\sum_{|\,|k| - k_r| \le \Delta k_r / 2} \hat f(k)

        See also:
            init_2Dwvnb, diagkit.stats.ke_spectrum
    """

    nk = kr.shape[0]  # number of wavenumbers
    nspace = 3 if fh.ndim >= 3 and np.prod(fh.shape[-3:]) == wv.size else 2
    fh = fh.reshape((fh.shape[:len(fh.shape) - nspace] + (-1,)))
    psdr = np.zeros((fh.shape[:len(fh.shape) - 1] + (nk,)), dtype=np.float64)

    # Summing over annular rings
    dkr = kr[1] - kr[0]
    for p in range(nk):
        # Integration between lower and upper bounds of annular rings
        idk = (wv >= kr[p] - dkr / 2) & (wv <= kr[p] + dkr / 2)
        psdr[..., p] = efc * np.sum(fh[..., idk], axis=-1, dtype=np.float64)
    return psdr


def temporal_psd(q, fs, weights=None):
    r"""
    temporal_psd: power spectral density in time, by periodogram, optionally summed over space with weights.

    The spectrum in time shows at which periods a field varies; summed over a region with weights it gives
    one spectrum for the whole region.

        Parameters:
            q         (float): Field, time along axis 0
            fs        (float): Sampling frequency [1/s]
            weights   (float): Optional 2D weight array broadcast against
                                q's trailing two axes; if given, the result
                                is summed over those axes.

        Returns:
            freq      (float): Frequency axis, shape (q.shape[0]//2 + 1,)
            psd       (float): PSD; same shape as q with axis 0 replaced by
                                freq if weights is None, else weights applied
                                and the trailing two axes summed away.

        Formula:
            P(f) = \text{periodogram}(q)(f),\qquad \text{or}\ \sum_{j,i} w_{j,i}\, P_{j,i}(f)

        See also:
            iso_spec
    """
    freq, psd = signal.periodogram(q, fs, axis=0)
    if weights is None:
        return freq, psd
    return freq, np.sum(psd * weights[None], axis=(-2, -1))
