#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
        Module Spectrums

        Computes the spectrum of a given field based on the procedure
        illustrated in Durran et al. (2017) 

          title:  " Practical considerations for computing  
                    dimensional spectra from gridded data "
            url:  " http://doi.org/10.1175/MWR-D-17-0056.1 "


         Author: Long Li,                           long.li@inria.fr
                 Francesco Tucciarone, francesco.tucciarone@inria.fr

        Created: 15 April 2022

""" 


import numpy as np

def init_2Dwvnb(nx, ny, dx, dy):
    """Initialize horizaontal and isotropic wavenumbers `wv` and `kr`"""

    assert nx % 2 == 0, 'Wavenumbers length in x-axis is not even'
    assert ny % 2 == 0, 'Wavenumbers length in y-axis is not even'

    # Create horizontal wavenumbers
    dkx = 2. * np.pi/(nx * dx)
    dky = 2. * np.pi/(ny * dy)
    kx = dkx * np.append(np.arange(0., nx/2), np.arange(-nx/2, 0.))
    ky = dky * np.append(np.arange(0., ny/2), np.arange(-ny/2, 0.))
    kxx, kyy = np.meshgrid(kx.astype('float64'), ky.astype('float64'))
    wv = np.sqrt(kxx**2 + kyy**2).flatten()

    # Create isotropic wavenumbers
    nmax = np.ceil(np.sqrt(2.) * np.maximum(nx, ny)/2.)
    nmax = np.minimum(nx,ny)/2.
    dkr = np.maximum(dkx, dky)
    kr = dkr * np.arange(1., nmax+1, dtype=np.float64)

    # Scaling factor for spectrum due to DFT and integration
    efac = dx * dy * np.minimum(dkx, dky)/(4. * np.pi**2 * nx * ny)
    return wv, kr, efac


def init_3Dwvnb(nx, ny, nz, dx, dy, dz):
    """Initialize horizaontal and isotropic wavenumbers `wv` and `kr`"""

    assert nx % 2 == 0, 'Wavenumbers length in x-axis is not even'
    assert ny % 2 == 0, 'Wavenumbers length in y-axis is not even'
    assert nz % 2 == 0, 'Wavenumbers length in z-axis is not even'

    # Create horizontal wavenumbers
    dkx = 2. * np.pi/(nx * dx)
    dky = 2. * np.pi/(ny * dy)
    dkz = 2. * np.pi/(nz * dz)

    kx = dkx * np.append(np.arange(0., nx/2), np.arange(-nx/2, 0.))
    ky = dky * np.append(np.arange(0., ny/2), np.arange(-ny/2, 0.))
    kz = dkz * np.append(np.arange(0., nz/2), np.arange(-nz/2, 0.))

    kxx, kyy, kzz = np.meshgrid(kx.astype('float64'),
                                ky.astype('float64'),
                                kz.astype('float64'))

    wv = np.sqrt(kxx**2 + kyy**2 + kzz**2).flatten()

    # Create isotropic wavenumbers
    nmax = np.ceil(np.sqrt(2.) * np.maximum(nx, ny, nz)/2.)
    # nmax = np.minimum(nx,ny)/2.
    dkr = np.maximum(dkx, dky)
    kr = dkr * np.arange(1., nmax+1, dtype=np.float64)

    # Scaling factor for spectrum due to DFT and integration
    efac = dx * dy * np.minimum(dkx, dky)/(4. * np.pi**2 * nx * ny)
    return wv, kr, efac




def iso_spec(fh, wv, kr, efc):
    """Compute isotropic spectrum `psdr` of 2D Fourier coef. `fh`"""
    
    nk = kr.shape[0]  # number of wavenumbers
    fh = fh.reshape((fh.shape[:len(fh.shape)-2] + (-1,)))
    psdr = np.zeros((fh.shape[:len(fh.shape)-1] + (nk,)), dtype=np.float64)

    # Summing over annular rings
    dkr = kr[1] - kr[0]
    for p in range(nk):
        # Integration between lower and upper bounds of annular rings
        idk = (wv >= kr[p]-dkr/2) & (wv <= kr[p]+dkr/2)
        psdr[..., p] = efc * np.sum(fh[..., idk], axis=-1, dtype=np.float64)
    return psdr

