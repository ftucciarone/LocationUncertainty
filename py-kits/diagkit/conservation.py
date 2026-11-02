#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module conservation

Conservation-checking routines for NEMO tracer fields: volume-weighted
averages over arbitrary depth ranges, and the 3D wet-point mask they need.

Promoted out of the earlier salinity conservation notebook (the volume-weighted
salinity check developed while debugging a surface-salinity drop), so the
same method can be reused for other tracers/runs without copy-pasting the
notebook cells.
"""

import numpy as np


def wet_mask_3d(bottom_level, jpk):
    r"""
    wet_mask_3d: the 3D wet-point mask from the number of wet levels of each column.

    Domain files store, for each column, how many levels are wet; volume integrals need the mask on every
    level. A level k, counted from 0, is wet when k + 1 <= bottom_level, so land columns, with 0, are dry
    everywhere.

        Parameters:
            bottom_level  (int): (y, x) array from domain_cfg.nc
            jpk           (int): number of vertical levels

        Returns:
            mask3d      (float): (k, y, x) wet-point mask, 1.0/0.0

        Formula:
            m_{k,j,i} = 1 \text{ if } k + 1 \le \text{bottom\_level}_{j,i},\ 0 \text{ otherwise}

        See also:
            volume_weighted_mean, orcakit.orca.create_tmask_from_top_and_bottom
    """
    k_idx = np.arange(1, jpk + 1)[:, None, None]
    return (k_idx <= bottom_level[None, :, :]).astype(np.float64)


def volume_weighted_mean(so, e3t, mask3d, area, k_lo=0, k_hi=None):
    r"""
    volume_weighted_mean: volume-weighted mean of a tracer over a range of levels, at every time.

    A conserved tracer's volume mean should not drift in a closed basin; tracking it, and the layer
    volume beside it, checks conservation. The layer volume should stay nearly constant.

        Parameters:
            so       (float): (time, k, y, x) tracer field
            e3t      (float): (time, k, y, x) cell thickness
            mask3d   (float): (k, y, x) wet-point mask, from wet_mask_3d
            area     (float): (y, x) cell area (e.g. e1t*e2t)
            k_lo       (int): first level of the range (inclusive)
            k_hi       (int): last level of the range (exclusive); defaults
                               to all remaining levels

        Returns:
            mean_s   (float): (time,) volume-weighted mean of `so`
            layer_volume (float): (time,) total volume of the layer - a
                               conservation sanity check, not just an
                               intermediate: it should stay ~constant.

        Formula:
            \bar s(t) = \frac{\sum_{k,j,i} s\, e_{3t}\, A\, m}{\sum_{k,j,i} e_{3t}\, A\, m},\qquad V(t) = \sum_{k,j,i} e_{3t}\, A\, m,\qquad k_{lo} \le k < k_{hi}

        See also:
            wet_mask_3d
    """
    if k_hi is None:
        k_hi = so.shape[1]
    m = mask3d[k_lo:k_hi]  # (k, y, x)
    vol = e3t[:, k_lo:k_hi] * area[None, :, :] * m[None]  # (time, k, y, x)
    salt = so[:, k_lo:k_hi] * vol
    vol_sum = vol.sum(axis=(1, 2, 3))
    mean_s = salt.sum(axis=(1, 2, 3)) / vol_sum
    return mean_s, vol_sum
