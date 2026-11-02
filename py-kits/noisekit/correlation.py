#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module correlation

Correlation lengths and times, by time-ensemble lagged correlation. Replaces the
efold, efold_frac, efold_pooled, efolding_lag and directional variants in the
O1 notebooks. V5's build_noise_lscale.py keeps its own MPI-tiled accumulator.
"""
import numpy as np


def efolding_lag(rho, kmax):
    r"""
    efolding_lag: the lag at which a correlation first drops below 1/e.

    The e-folding lag is a simple measure of how far, in space or time, a field stays correlated with
    itself. The crossing is interpolated linearly between the two lags around it, with lag 0 taken as 1.

    Where the correlation never drops below 1/e within kmax, the true lag is at least kmax: it is
    returned as kmax and flagged as censored rather than given a false precise value. NaN never counts
    as a crossing.

        Parameters:
            rho  (float): (kmax, ...) correlation at lags 1 to kmax
            kmax   (int): Largest lag, returned where there is no crossing

        Returns:
            L        (float): (...) e-folding lag, in lag units; kmax where censored
            censored  (bool): (...) True where rho never crossed 1/e

        Formula:
            L = (k-1) + \frac{\rho_{k-1} - e^{-1}}{\rho_{k-1} - \rho_k},\qquad k = \min\{k : \rho_k < e^{-1}\},\quad \rho_0 = 1

        See also:
            spatial_correlation, temporal_correlation
    """
    below = rho < (1.0 / np.e)
    censored = ~below.any(axis=0)
    first_k = np.argmax(below, axis=0) + 1
    rho_k = np.take_along_axis(rho, (first_k - 1)[None], axis=0)[0]
    rho_km1 = np.where(first_k == 1, 1.0, np.take_along_axis(rho, (first_k - 2).clip(min=0)[None], axis=0)[0])
    denom = rho_km1 - rho_k
    frac = np.where(denom > 1e-12, (rho_km1 - 1.0 / np.e) / denom, 0.0)
    return np.where(censored, float(kmax), (first_k - 1) + frac), censored


def spatial_correlation(fields, kmax):
    r"""
    spatial_correlation: lagged correlation in space, along x and along y.

    The correlation between a point and the point k cells away, over time, measures the spatial scale of
    a field's variability. Time serves as the ensemble. Several fields, such as u' and v', can be pooled:
    their lagged products are summed before dividing by their pooled variances.

    The fields are taken as anomalies; their mean is not removed. Lags that run off the array are NaN,
    and points with no variance get 0. Each pair is normalised by the variances at both its ends, which
    keeps the correlation between -1 and 1.

        Parameters:
            fields (list): (nt, ny, nx) arrays
            kmax    (int): Largest lag, in grid points

        Returns:
            rho_x (float): (kmax, ny, nx) correlation with the point k columns east; NaN where the
                           lag runs off the array, 0 where the variance vanishes
            rho_y (float): (kmax, ny, nx) correlation with the point k rows north; same conventions

        Formula:
            \rho_x(k) = \frac{\sum_f \langle f_{i,j}\, f_{i+k,j} \rangle_t}{\sqrt{\sum_f \langle f_{i,j}^2 \rangle_t \ \sum_f \langle f_{i+k,j}^2 \rangle_t}},\qquad \rho_y \text{ likewise along } j

        See also:
            efolding_lag, temporal_correlation
    """
    P0 = sum((f**2).mean(axis=0) for f in fields)
    ny, nx = P0.shape
    rho_x, rho_y = np.full((kmax, ny, nx), np.nan), np.full((kmax, ny, nx), np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        for k in range(1, kmax+1):
            Px = sum((f[:, :, :-k]*f[:, :, k:]).mean(axis=0) for f in fields)
            d = np.sqrt(P0[:, :-k]*P0[:, k:])
            rho_x[k-1, :, :-k] = np.where(d > 1e-20, Px/d, 0.0)
            Py = sum((f[:, :-k, :]*f[:, k:, :]).mean(axis=0) for f in fields)
            d = np.sqrt(P0[:-k, :]*P0[k:, :])
            rho_y[k-1, :-k, :] = np.where(d > 1e-20, Py/d, 0.0)
    return rho_x, rho_y


def temporal_correlation(field, kmax):
    r"""
    temporal_correlation: lagged correlation in time, point by point.

    How long a field stays correlated with itself sets the time scale of its variability, for example the
    memory of a noise process. The time mean is removed first.

    Points with no variance get 0.

        Parameters:
            field (float): (nt, ...) record
            kmax    (int): Largest lag, in snapshots

        Returns:
            rho (float): (kmax, ...) correlation at lags 1 to kmax; 0 where the variance vanishes

        Formula:
            \rho(k) = \frac{\langle a_t\, a_{t+k} \rangle}{\langle a_t^2 \rangle},\qquad a = f - \bar f

        See also:
            efolding_lag, spatial_correlation
    """
    a = field - field.mean(axis=0, keepdims=True)
    P0 = (a**2).mean(axis=0)
    out = np.empty((kmax,) + a.shape[1:])
    with np.errstate(invalid="ignore", divide="ignore"):
        for k in range(1, kmax+1):
            out[k-1] = np.where(P0 > 1e-20, (a[:-k]*a[k:]).mean(axis=0)/P0, 0.0)
    return out
