#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module filters

Filters as primitives: arrays in, arrays out, for any field at any step.
Space: lowpass, highpass, lowpass_tile, and masked_smooth and scale_split, which keep land as 0.
Time: deseason, time_lowpass.

Scale separation via spatial low-pass filtering: a mask-aware smoothing
that isolates the large-scale field, and (field - large_scale) as the
small-scale residual. Grid-generic -- the actual smoothing kernel (e.g.
ORCA's Gaussian smooth_2D with north-fold handling) is injected as a
callable, same pattern as grid_c.py's `exchange`.

    smooth signature: smooth(field, sigma) -> smoothed_field
"""

import numpy as np


def masked_smooth(field_in, mask, sigma, smooth, smooth_mask=None):
    r"""
    masked_smooth: low-pass a field without letting masked points, such as land, bias it near the mask edge.

    Smoothing a field with zeros on land pulls the result towards zero along every coast, since the
    kernel averages those zeros in. Here the field times the mask and the mask itself are smoothed with
    the same kernel and divided, which is a local average over valid points only.

    Where under 1% of the kernel's footprint is valid the ratio divides a tiny number by a tiny number:
    near the ORCA north fold this once turned v = -0.006 m/s into -325 m/s. Those points are set to 0,
    as are masked points. lowpass does the same with NaN in place of 0.

        Parameters:
            field_in    (float): Field to be smoothed
            mask        (float): 1 for valid (e.g. ocean) points, 0
                                  elsewhere; broadcastable against
                                  field_in's trailing two axes
            sigma       (float): Filter lengthscale, in gridpoints
            smooth   (callable): smooth(field, sigma) -> smoothed_field,
                                  e.g. functools.partial(orca.smooth_2D,
                                  mask=None) -- this routine handles
                                  masking itself, so the injected smooth
                                  should not also mask internally
            smooth_mask (callable): smooth for the mask, default `smooth`.
                                  Needed when `smooth` flips the sign at the
                                  north fold: a mask must never flip

        Returns:
            field_out   (float): Smoothed field, zeroed both on masked
                                  points and wherever the filter footprint
                                  found no valid neighbour at all

        Formula:
            \tilde f = \frac{G_\sigma * (f\,m)}{G_\sigma * m}\;m \quad\text{where } G_\sigma * m > 0.01,\ \text{else } 0

        See also:
            lowpass, scale_split
    """
    num = smooth(field_in * mask, sigma)
    den = (smooth_mask or smooth)(mask, sigma)
    # Below this, fewer than ~1% of the filter footprint was valid data --
    # not exactly 0, but not enough to trust the ratio either.
    den_floor = 0.01
    with np.errstate(invalid="ignore", divide="ignore"):
        field_out = num / den
    field_out = np.where(den > den_floor, field_out, 0.0)
    return field_out * mask


def scale_split(field_in, mask, sigma, smooth, smooth_mask=None):
    r"""
    scale_split: split a field into a large-scale part and the small-scale residual.

    A field is often studied as a resolved part plus what a coarser model cannot see. The large-scale
    part is masked_smooth of the field and the residual is the rest, so the two add back to the field
    on valid points.

    Both parts are 0 on masked points, as in masked_smooth.

        Parameters:
            field_in    (float): Field to split
            mask        (float): 1 for valid points, 0 elsewhere
            sigma       (float): Filter lengthscale, in gridpoints --
                                  the large/small-scale cutoff
            smooth   (callable): see masked_smooth
            smooth_mask (callable): see masked_smooth

        Returns:
            large       (float): Low-pass (large-scale) field, masked
            small       (float): field_in - large, masked (small-scale
                                  residual -- the noise-generating part)

        Formula:
            f = \tilde f + f',\qquad \tilde f = \text{masked\_smooth}(f),\qquad f' = (f - \tilde f)\,m

        See also:
            masked_smooth, highpass
    """
    large = masked_smooth(field_in, mask, sigma, smooth, smooth_mask)
    small = (field_in - large) * mask
    return large, small


def gaussian(field, sigma):
    r"""
    gaussian: Gaussian smoothing in the last two axes of a field that already carries its halo.

    This is the plain kernel, with no knowledge of the ORCA seam or fold: it relies on a halo filled
    beforehand, from the neighbouring tiles or from disk. Pass it as `smooth` to lowpass to filter a
    padded tile.

    Edges are extended with their nearest value, which matters only if the halo is narrower than the
    kernel, 4 sigma.

        Parameters:
            field (float): (..., ny, nx) padded field
            sigma (float): Filter lengthscale, in grid points

        Returns:
            field_out (float): Smoothed field, same shape

        Formula:
            (G_\sigma * f)_{i,j} = \sum_{m,n} \frac{e^{-(m^2+n^2)/2\sigma^2}}{\sum e^{-(m^2+n^2)/2\sigma^2}}\, f_{i+m,\,j+n}

        See also:
            lowpass, lowpass_tile
    """
    from scipy.ndimage import gaussian_filter
    sig = np.zeros(np.ndim(field))
    sig[-2:] = sigma
    return gaussian_filter(field, sigma=sig, mode="nearest")


def lowpass(field, mask, sigma, point="T", vector=False, smooth=None):
    r"""
    lowpass: Gaussian low-pass in space over valid points only, for fields or statistics.

    A Gaussian low-pass keeps the scales larger than sigma and removes the rest. On an ocean grid it must
    ignore land, or every coast is pulled towards zero: the field and the mask of valid points are
    smoothed with the same kernel and divided. The same call filters a velocity, a mean or a variance;
    only the sign at the north fold differs.

    Where under 1% of the footprint is valid the ratio is not trusted, and the result is NaN there, as on
    land. By default the ORCA Gaussian handles the seam and the fold; a custom `smooth` can replace it,
    for example on a tile that already has its halo.

        Parameters:
            field  (float): (..., ny, nx); NaN counts as invalid
            mask    (bool): True where valid, broadcastable to field
            sigma  (float): Filter lengthscale, in grid points
            point    (str): 'T', 'U', 'V' or 'F', where the field lives
            vector  (bool): True for a velocity component, which flips sign at the north fold
            smooth (callable): smooth(field, sigma); default the ORCA Gaussian (orca.smooth_2D)
                               with this point and sign. The mask is smoothed with it too, so a
                               custom one must not flip sign

        Returns:
            field_out (float): Low-passed field, NaN on invalid points and where under 1% of the
                               footprint was valid

        Formula:
            \tilde f = \frac{G_\sigma * (f\,m)}{G_\sigma * m} \quad\text{where } m = 1 \text{ and } G_\sigma * m > 0.01,\ \text{else NaN}

        See also:
            highpass, lowpass_tile, masked_smooth, orca.smooth_2D
    """
    field = np.asarray(field)
    valid = np.asarray(mask, dtype=bool) & np.isfinite(field)
    if smooth is None:
        from . import orca
        num = orca.smooth_2D(np.where(valid, field, 0.0), sigma, isMeridional=vector, point=point)
        den = orca.smooth_2D(valid.astype(float), sigma, point=point)
    else:
        num, den = smooth(np.where(valid, field, 0.0), sigma), smooth(valid.astype(float), sigma)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(valid & (den > 0.01), num/den, np.nan)


def highpass(field, mask, sigma, point="T", vector=False, smooth=None):
    r"""
    highpass: Gaussian high-pass in space, the field minus its low-pass.

    The high-pass keeps the scales smaller than sigma, such as eddies and fronts, and is what is left
    once lowpass has removed the large scales. lowpass and highpass add back to the field on valid points.

    It is NaN wherever lowpass is.

        Parameters:
            field  (float): (..., ny, nx); NaN counts as invalid
            mask    (bool): True where valid, broadcastable to field
            sigma  (float): Filter lengthscale, in grid points
            point    (str): 'T', 'U', 'V' or 'F', where the field lives
            vector  (bool): True for a velocity component, which flips sign at the north fold
            smooth (callable): See lowpass

        Returns:
            field_out (float): High-passed field, NaN where lowpass is NaN

        Formula:
            f' = f - \tilde f,\qquad \tilde f = \text{lowpass}(f)

        See also:
            lowpass
    """
    return np.asarray(field) - lowpass(field, mask, sigma, point, vector, smooth)


def lowpass_tile(dec, field, mask, sigma, point="T", vector=False):
    r"""
    lowpass_tile: lowpass on one rank's tile of a domain decomposition.

    A whole high-resolution grid with many levels may not fit in one process. Each rank then filters its
    own tile: a halo of width ceil(4 sigma) is filled from the neighbouring tiles, the padded tile is
    filtered, and the halo is stripped, which gives the same result as filtering the whole grid.

    With a serial decomposition it filters the whole grid. The halo, plus 5 rows for the fold, must fit
    in the smallest tile.

        Parameters:
            dec    (Decomp): mpi_orca decomposition; 4 sigma + 5 must fit in its smallest tile
            field   (float): (..., jl, il) this rank's tile; NaN counts as invalid
            mask     (bool): True where valid, (jl, il) or the shape of field
            sigma   (float): Filter lengthscale, in grid points
            point     (str): 'T', 'U', 'V' or 'F', where the field lives
            vector   (bool): True for a velocity component, which flips sign at the north fold

        Returns:
            field_out (float): (..., jl, il) low-passed tile, NaN as in lowpass

        Formula:
            \tilde f = \frac{G_\sigma * (f\,m)}{G_\sigma * m} \text{ on the tile, with a halo of } \lceil 4\sigma \rceil \text{ from the neighbours}

        See also:
            lowpass, mpi_orca.add_halo
    """
    from . import mpi_orca
    hw = int(np.ceil(4*sigma))
    field = np.asarray(field)
    valid = np.broadcast_to(np.asarray(mask, dtype=bool), field.shape).astype(field.dtype)
    a = mpi_orca.add_halo(dec, np.where(np.isfinite(field), field, 0.0).astype(field.dtype), hw, is_meridional=vector, point=point)
    w = mpi_orca.add_halo(dec, np.where(np.isfinite(field), valid, 0.0).astype(field.dtype), hw, point=point)
    return mpi_orca.strip_halo(lowpass(a, w > 0.5, sigma, smooth=gaussian), hw)


def deseason(field, time_days, harmonics=2, period=365.25):
    r"""
    deseason: remove, point by point, the mean, a linear trend and the harmonics of a period.

    The seasonal cycle and a slow drift often dominate a record but are not the variability of interest.
    This fits them by least squares at every point and returns what is left. Unlike a broad high-pass
    in time, it removes only those components and leaves every other time scale as it was.

    All points are fitted in one small least-squares system. A point with any NaN in time stays NaN. The
    record must be longer than the number of coefficients, and a few years are needed to separate the
    trend from the seasons.

        Parameters:
            field     (float): (nt, ...) record
            time_days (float): (nt,) time, in days
            harmonics   (int): Number of harmonics; 2 for annual and semiannual
            period    (float): Period, in days

        Returns:
            residual (float): (nt, ...) field minus the fit
            coef     (float): (2 + 2*harmonics, ...) mean, trend over the record, then the cosine and
                              sine of each harmonic. Raises ValueError if nt is too short for them

        Formula:
            f(t) \approx a_0 + a_1\,\frac{t - t_0}{t_{\mathrm{end}} - t_0} + \sum_{k=1}^{K} \big[b_k \cos(k\omega t) + c_k \sin(k\omega t)\big],\qquad \omega = \frac{2\pi}{P}
            \text{residual} = f - \text{fit, by least squares at each point}

        See also:
            time_lowpass
    """
    field, t = np.asarray(field), np.asarray(time_days, dtype=float)
    omega = 2*np.pi/period
    columns = [np.ones_like(t), (t - t[0])/max(np.ptp(t), 1.0)]
    for k in range(1, harmonics + 1):
        columns += [np.cos(k*omega*t), np.sin(k*omega*t)]
    design = np.stack(columns, axis=1)
    if design.shape[0] <= design.shape[1]:
        raise ValueError(f"{design.shape[0]} snapshots cannot fit {design.shape[1]} coefficients")
    good = np.isfinite(field).all(axis=0)
    coef = np.full((design.shape[1],) + field.shape[1:], np.nan)
    residual = np.full(field.shape, np.nan)
    coef[:, good] = np.linalg.lstsq(design, field[:, good], rcond=None)[0]
    residual[:, good] = field[:, good] - design @ coef[:, good]
    return residual, coef


def time_lowpass(field, dt, cutoff_period, order=5):
    r"""
    time_lowpass: zero-phase Butterworth low-pass in time, along axis 0.

    A Butterworth filter has a flat response below its cutoff and falls off quickly above it. Run forward
    and then backward, it shifts no event in time, which matters when a filtered series is compared with
    the original. The high-pass is the field minus the result.

    The start and end of the record are less reliable, since the filter needs data on both sides. A point
    with any NaN in time stays NaN.

        Parameters:
            field         (float): (nt, ...) record; points with any NaN in time stay NaN
            dt            (float): Time between snapshots
            cutoff_period (float): Shortest period kept, in the unit of dt
            order           (int): Order of the filter

        Returns:
            field_out (float): (nt, ...) low-passed record

        Formula:
            |H(\omega)|^2 = \Big[1 + (\omega/\omega_c)^{2N}\Big]^{-1},\qquad \omega_c = \frac{2\pi}{T_c}
            \text{applied forward and backward, so the gain is } |H|^2 \text{ and there is no phase shift}

        See also:
            deseason, highpass
    """
    from scipy import signal
    field = np.asarray(field)
    b, a = signal.butter(order, 2*dt/cutoff_period, "low")
    good = np.isfinite(field).all(axis=0)
    out = np.full(field.shape, np.nan)
    out[:, good] = signal.filtfilt(b, a, field[:, good], axis=0)
    return out
