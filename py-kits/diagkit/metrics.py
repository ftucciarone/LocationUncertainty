#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module metrics

Weighted comparison metrics. Spatial axes are the trailing weights.ndim axes.
Ported from DiagnosticMetrics of the earlier analysis code.
"""
import numpy as np


def _wmean(q, weights):
    r"""
    _wmean: weighted mean over the trailing axes covered by the weights.

        Parameters:
            q       (float): Field, whose last axes match weights
            weights (float): Weights

        Returns:
            mean (float): Weighted mean, with the leading axes of q

        Formula:
            \langle q \rangle_w = \frac{\sum q\, w}{\sum w}

        See also:
            rms
    """
    axes = tuple(range(-weights.ndim, 0))
    return np.sum(q*weights, axis=axes) / np.sum(weights, axis=axes)


def rms(observable, weights):
    r"""
    rms: weighted root mean square.

        Parameters:
            observable (float): Field
            weights    (float): Weights over its trailing axes

        Returns:
            rms (float): sqrt of the weighted mean of observable^2

        Formula:
            \sqrt{\langle o^2 \rangle_w}

        See also:
            rmse
    """
    return np.sqrt(_wmean(observable**2, weights))


def rmse(estimator, observable, weights):
    r"""
    rmse: weighted root mean square error of an estimate against a reference.

    The usual measure of how far a model or a reduced model is from the reference, in the field's units.

        Parameters:
            estimator  (float): Estimate
            observable (float): Reference, same shape
            weights    (float): Weights over the trailing axes

        Returns:
            rmse (float): rms of estimator - observable

        Formula:
            \sqrt{\langle (e - o)^2 \rangle_w}

        See also:
            rms, pattern_correlation
    """
    return rms(estimator - observable, weights)


def pattern_correlation(estimator, observable, weights):
    r"""
    pattern_correlation: weighted, uncentred pattern correlation, 1 for identical patterns.

    Measures whether two fields have the same shape, regardless of amplitude. An earlier version had no
    square root in the denominator, so it was not bounded by 1.

        Parameters:
            estimator  (float): Estimate
            observable (float): Reference, same shape
            weights    (float): Weights over the trailing axes

        Returns:
            corr (float): <e o> / sqrt(<e^2> <o^2>), weighted means

        Formula:
            \frac{\langle e\, o \rangle_w}{\sqrt{\langle e^2 \rangle_w\, \langle o^2 \rangle_w}}

        See also:
            rmse
    """
    num = _wmean(estimator*observable, weights)
    return num / np.sqrt(_wmean(estimator**2, weights) * _wmean(observable**2, weights))


def gaussian_relative_entropy(est_avg, obs_avg, est_vnc, obs_vnc):
    r"""
    gaussian_relative_entropy: relative entropy of a Gaussian estimate against a Gaussian reference, point by point.

    It measures how much information is lost when the estimate replaces the reference, penalising both a
    wrong mean and a wrong spread, which suits comparing ensembles with their target.

        Parameters:
            est_avg (float): Mean of the estimate
            obs_avg (float): Mean of the reference
            est_vnc (float): Variance of the estimate
            obs_vnc (float): Variance of the reference

        Returns:
            gre (float): 0.5 [(est_avg - obs_avg)^2 / est_vnc + r - 1 - log r], r = obs_vnc / est_vnc

        Formula:
            \mathcal{R} = \tfrac12\Big[\frac{(\mu_e - \mu_o)^2}{\sigma_e^2} + \frac{\sigma_o^2}{\sigma_e^2} - 1 - \log\frac{\sigma_o^2}{\sigma_e^2}\Big]

        See also:
            averaged_gre
    """
    tmp1 = np.divide((est_avg - obs_avg)**2, est_vnc)
    tmp2 = np.divide(obs_vnc, est_vnc)
    return 0.5 * (tmp1 + tmp2 - (1 + np.log(tmp2)))


def averaged_gre(est_avg, obs_avg, est_vnc, obs_vnc, weights):
    r"""
    averaged_gre: weighted mean of the Gaussian relative entropy.

        Parameters:
            est_avg (float): Mean of the estimate
            obs_avg (float): Mean of the reference
            est_vnc (float): Variance of the estimate
            obs_vnc (float): Variance of the reference
            weights (float): Weights over the trailing axes

        Returns:
            gre (float): Weighted mean

        Formula:
            \langle \mathcal{R} \rangle_w

        See also:
            gaussian_relative_entropy
    """
    return _wmean(gaussian_relative_entropy(est_avg, obs_avg, est_vnc, obs_vnc), weights)
