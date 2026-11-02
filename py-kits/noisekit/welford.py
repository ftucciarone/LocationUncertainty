#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module welford

One-pass (online) mean/variance/covariance accumulation via Welford's
algorithm. Works array-wise: every aggregate entry is a full field
(e.g. shape (nlev, nlat, nlon)), updated in place sample-by-sample, so
this never needs to hold the whole time series in memory at once.
"""

import numpy as np


def Welford_mean_aggregate_init(array_in):
    r"""
    Welford_mean_aggregate_init: a fresh (count, mean) aggregate for a running mean.

    Welford's method computes statistics in one pass, one sample at a time, so a record too long for
    memory can be averaged. This starts it.

        Parameters:
            array_in (float): template array giving the aggregate's shape

        Returns:
            aggregate (tuple): (count, mean)

        Formula:
            (n, \mu) = (0, 0)

        See also:
            Welford_mean, Welford_mean_finalize
    """
    count = 0
    mean = np.zeros_like(array_in)
    return (count, mean)


def Welford_mean(previous_aggregate, new_sample):
    r"""
    Welford_mean: add one sample to a running mean.

        Parameters:
            previous_aggregate (tuple): (count, mean) from the previous step
            new_sample         (float): next sample field

        Returns:
            aggregate (tuple): updated (count, mean)

        Formula:
            n \leftarrow n + 1,\qquad \mu \leftarrow \mu + \frac{x - \mu}{n}

        See also:
            Welford_mean_aggregate_init, Welford_mean_finalize
    """
    (count, mean) = previous_aggregate
    # Update counter
    count = count + 1
    # Update mean
    delta_1 = new_sample - mean
    mean = mean + delta_1 / count
    return (count, mean)


def Welford_mean_finalize(previous_aggregate):
    r"""
    Welford_mean_finalize: the mean of a (count, mean) aggregate.

    NaN with fewer than 2 samples, like the variance.

        Parameters:
            previous_aggregate (tuple): (count, mean)

        Returns:
            mean (float): the mean field, or NaN if count < 2

        Formula:
            \mu

        See also:
            Welford_mean
    """
    (count, mean) = previous_aggregate
    if count < 2:
        return float("nan")
    else:
        return mean


def Welford_aggregate_init(array_in):
    r"""
    Welford_aggregate_init: a fresh (count, mean, M2) aggregate for a running variance.

    Welford's update of the variance is numerically stable, unlike summing x and x^2 separately, and
    needs one sample at a time.

        Parameters:
            array_in (float): template array giving the aggregate's shape

        Returns:
            aggregate (tuple): (count, mean, M2)

        Formula:
            (n, \mu, M_2) = (0, 0, 0)

        See also:
            Welford_variance, Welford_finalize
    """
    count = 0
    mean = np.zeros_like(array_in)
    M2 = np.zeros_like(array_in)
    return (count, mean, M2)


def Welford_variance(previous_aggregate, new_sample):
    r"""
    Welford_variance: add one sample to a running mean and variance.

        Parameters:
            previous_aggregate (tuple): (count, mean, M2) from the previous step
            new_sample         (float): next sample field

        Returns:
            aggregate (tuple): updated (count, mean, M2)

        Formula:
            n \leftarrow n + 1,\quad \delta = x - \mu,\quad \mu \leftarrow \mu + \frac{\delta}{n},\quad M_2 \leftarrow M_2 + \delta\,(x - \mu)

        See also:
            Welford_aggregate_init, Welford_finalize
    """
    (count, mean, M2) = previous_aggregate
    # Update counter
    count = count + 1
    # Update mean
    delta_1 = new_sample - mean
    mean = mean + delta_1 / count
    # Update M2
    delta_2 = new_sample - mean
    M2 = M2 + delta_1 * delta_2
    return (count, mean, M2)


def Welford_finalize(previous_aggregate):
    r"""
    Welford_finalize: mean, population variance and sample variance of a (count, mean, M2) aggregate.

    NaN when there are too few samples.

        Parameters:
            previous_aggregate (tuple): (count, mean, M2)

        Returns:
            result (tuple): (mean, variance, sample_variance), or NaN if
                             count < 2

        Formula:
            \mu,\qquad \sigma^2 = \frac{M_2}{n},\qquad s^2 = \frac{M_2}{n - 1}

        See also:
            Welford_variance
    """
    (count, mean, M2) = previous_aggregate
    if count < 2:
        return float("nan")
    else:
        (mean, variance, sample_variance) = (mean, M2 / count, M2 / (count - 1))
        return (mean, variance, sample_variance)


def Welford_aggregate_xy_init(array_in):
    r"""
    Welford_aggregate_xy_init: a fresh aggregate for the running means, variances and covariance of two fields.

        Parameters:
            array_in (float): template array giving the aggregate's shape

        Returns:
            aggregate (tuple): (count, mean_x, mean_y, M2x, M2y, C2xy)

        Formula:
            (n, \mu_x, \mu_y, M_{2x}, M_{2y}, C_{xy}) = 0

        See also:
            Welford_covariance, Welford_cov_finalize
    """
    count = 0
    mean_x = np.zeros_like(array_in)
    mean_y = np.zeros_like(array_in)
    M2x = np.zeros_like(array_in)
    M2y = np.zeros_like(array_in)
    C2xy = np.zeros_like(array_in)
    return (count, mean_x, mean_y, M2x, M2y, C2xy)


def Welford_covariance(previous_aggregate_xy, new_sample_x, new_sample_y):
    r"""
    Welford_covariance: add one pair of samples to a running mean, variance and covariance.

    The covariance update uses the mean of x before and of y after the update, which keeps it exact.

        Parameters:
            previous_aggregate_xy (tuple): (count, mean_x, mean_y, M2x, M2y, C2xy)
            new_sample_x          (float): next sample field for x
            new_sample_y          (float): next sample field for y

        Returns:
            aggregate (tuple): updated (count, mean_x, mean_y, M2x, M2y, C2xy)

        Formula:
            n \leftarrow n+1,\quad \delta_x = x - \mu_x,\quad \mu_x \leftarrow \mu_x + \tfrac{\delta_x}{n},\quad \mu_y \leftarrow \mu_y + \tfrac{y - \mu_y}{n},\quad C_{xy} \leftarrow C_{xy} + \delta_x\,(y - \mu_y)

        See also:
            Welford_aggregate_xy_init, Welford_cov_finalize
    """
    (count, mean_x, mean_y, M2x, M2y, C2xy) = previous_aggregate_xy
    # Update counter
    count = count + 1
    # Update means
    delta_x = new_sample_x - mean_x
    mean_x = mean_x + delta_x / count
    delta_y = new_sample_y - mean_y
    mean_y = mean_y + delta_y / count
    # Update M2x
    delta_2x = new_sample_x - mean_x
    M2x = M2x + delta_x * delta_2x
    # Update M2y
    delta_2y = new_sample_y - mean_y
    M2y = M2y + delta_y * delta_2y
    # Update C2 (mixed old/new delta, NOT delta_x * delta_y: using the
    # pre-update delta_y here would add a spurious pre-update product on
    # every step and break the online covariance identity, most visibly
    # at count == 1 where it would add new_sample_x * new_sample_y
    # instead of the correct contribution of 0)
    C2xy = C2xy + delta_x * delta_2y
    return (count, mean_x, mean_y, M2x, M2y, C2xy)


def Welford_cov_finalize(previous_aggregate_xy):
    r"""
    Welford_cov_finalize: means, population variances and covariance of a joint aggregate.

    NaN when there are too few samples.

        Parameters:
            previous_aggregate_xy (tuple): (count, mean_x, mean_y, M2x, M2y, C2xy)

        Returns:
            result (tuple): (mean_x, mean_y, var_x, var_y, cov_xy), or NaN
                             if count < 1

        Formula:
            \mu_x,\ \mu_y,\qquad \sigma_x^2 = \frac{M_{2x}}{n},\ \ \sigma_y^2 = \frac{M_{2y}}{n},\qquad c_{xy} = \frac{C_{xy}}{n}

        See also:
            Welford_covariance
    """
    (count, mean_x, mean_y, M2x, M2y, C2xy) = previous_aggregate_xy
    if count < 1:
        return float("nan")
    else:
        (mean_x, mean_y, var_x, var_y, cov_xy) = (mean_x, mean_y, M2x / count, M2y / count, C2xy / count)
        return (mean_x, mean_y, var_x, var_y, cov_xy)
