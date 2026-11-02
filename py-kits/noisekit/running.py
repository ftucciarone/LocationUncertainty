#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module running

Statistics over a sliding window of snapshots.
"""
from collections import deque

import numpy as np


def rolling_stats(snapshots, window, stride=1):
    r"""
    rolling_stats: mean and standard deviation over a sliding window of snapshots.

    Statistics over a short window follow how variability itself changes in time, for example from season
    to season. Only `window` snapshots are held at once, so a long record can be streamed.

    The standard deviation is the population one, divided by the window length.

        Parameters:
            snapshots (iterable): Arrays of one shape, in time order
            window         (int): Snapshots per window
            stride         (int): Keep one window in every `stride`

        Returns:
            generator: (index, mean, std) per kept window, index being that of its newest
                       snapshot; std is the population standard deviation

        Formula:
            \mu_t = \frac{1}{W}\sum_{s=t-W+1}^{t} x_s,\qquad \sigma_t^2 = \frac{1}{W}\sum_{s=t-W+1}^{t} (x_s - \mu_t)^2,\qquad t = W-1,\ W-1+S,\ \dots

        See also:
            running_stat, running_cov
    """
    buf = deque(maxlen=window)
    for i, snap in enumerate(snapshots):
        buf.append(snap)
        if len(buf) == window and (i - window + 1) % stride == 0:
            stack = np.stack(buf)
            yield i, stack.mean(axis=0), stack.std(axis=0)


def running_stat(snapshots, window, stride=1, stat="var"):
    r"""
    running_stat: one statistic over a sliding window of snapshots.

    The mean, variance or standard deviation of rolling_stats, one at a time.

        Parameters:
            snapshots (iterable): Arrays of one shape, in time order
            window         (int): Snapshots per window
            stride         (int): Keep one window in every `stride`
            stat           (str): 'var', 'std' or 'mean'; population variance and std

        Returns:
            generator: (index, statistic) per kept window, index being that of its newest snapshot

        Formula:
            \mu_t,\ \sigma_t^2 \text{ or } \sigma_t \text{ over } x_{t-W+1}, \dots, x_t

        See also:
            rolling_stats, running_cov
    """
    for i, mean, std in rolling_stats(snapshots, window, stride):
        yield i, {"mean": mean,
                  "std": std,
                  "var": std**2}[stat]


def running_cov(pairs, window, stride=1):
    r"""
    running_cov: covariance of two fields over a sliding window.

    How two components vary together, for example u and v, sets the orientation of variability. The two
    fields must be at the same points: on the C grid, put U and V on T-points first.

    Population covariance, as running_stat's variance. NaN in either field gives NaN.

        Parameters:
            pairs (iterable): (x, y) arrays of one shape, in time order
            window     (int): Pairs per window
            stride     (int): Keep one window in every `stride`

        Returns:
            generator: (index, cov) per kept window, index being that of its newest pair; NaN
                       where x or y is NaN in the window

        Formula:
            c_t = \frac{1}{W}\sum_{s=t-W+1}^{t} (x_s - \bar x_t)(y_s - \bar y_t)

        See also:
            running_stat, orcakit.grid_c.interp_UV2T
    """
    buf = deque(maxlen=window)
    for i, (x, y) in enumerate(pairs):
        buf.append((x, y))
        if len(buf) == window and (i - window + 1) % stride == 0:
            xs, ys = np.stack([b[0] for b in buf]), np.stack([b[1] for b in buf])
            yield i, ((xs - xs.mean(axis=0))*(ys - ys.mean(axis=0))).mean(axis=0)
