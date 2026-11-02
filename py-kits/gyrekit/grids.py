#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module grids

GYRE masks and high-to-low resolution helpers. Ported from
the mask and filter modules of the earlier GYRE analysis code.
coarse_grain_T was not ported: it needs scipy's removed interp2d.
"""
import numpy as np
import scipy.ndimage


def jet_mask(lat, lon, dx, dy, latc=29, lonc=-85, latw=5, lonw=20):
    r"""
    jet_mask: a mask of the box around the GYRE jet.

    1 inside a latitude-longitude box around the jet, 0 outside it and within 150 grid units of the north
    and west edges, where boundary currents would dominate.

        Parameters:
            lat  (float): (y, x) latitudes, in degrees
            lon  (float): (y, x) longitudes, in degrees
            dx   (float): Spacing in x
            dy   (float): Spacing in y
            latc (float): Centre latitude of the box
            lonc (float): Centre longitude of the box
            latw (float): Half-width of the box in latitude
            lonw (float): Half-width of the box in longitude

        Returns:
            mask (int): (y, x) mask

        Formula:
            m = 1 \text{ if } |\varphi - \varphi_c| < w_\varphi,\ |\lambda - \lambda_c| < w_\lambda, \text{ away from the north and west edges}

        See also:
            area_int
    """
    mask = ((lat > latc-latw) & (lat < latc+latw) & (lon > lonc-lonw) & (lon < lonc+lonw)).astype(int)
    mask[-int(np.floor(150/dy)):, :] = 0
    mask[:, :int(np.floor(150/dx))] = 0
    return mask


def find_overlap(HRlat, HRlon, LRlat, LRlon):
    r"""
    find_overlap: columns and rows of a high-resolution grid that match the lines of a low-resolution one.

    Original logic, not verified.

        Parameters:
            HRlat (float): (y, x) high-resolution latitudes
            HRlon (float): (y, x) high-resolution longitudes
            LRlat (float): (y, x) low-resolution latitudes
            LRlon (float): (y, x) low-resolution longitudes

        Returns:
            cols (int): Matching columns of the high-resolution grid
            rows (int): Matching rows of the high-resolution grid

        Formula:
            \text{columns and rows where repeated low-resolution coordinates match the high-resolution ones}

        See also:
            coarse_grain_P
    """
    HRny, HRnx = HRlat.shape
    LRny, LRnx = LRlat.shape
    xR = (HRnx - 2) // (LRnx - 2)
    yR = (HRny - 2) // (LRny - 2)
    repLon = np.repeat(np.repeat(LRlon[1:-1, 1:-1], yR, axis=0), xR, axis=1)
    repLat = np.repeat(np.repeat(LRlat[1:-1, 1:-1], yR, axis=0), xR, axis=1)
    mask = (repLon - HRlon[1:-1, 1:-1]) == (repLat - HRlat[1:-1, 1:-1])
    cols = np.arange(HRnx)[1:-1]
    rows = np.arange(HRny)[1:-1]
    return cols[np.sum(mask, axis=0) < 60], rows[np.sum(mask, axis=1) < 1]


def area_int(p, wefc, snfc):
    r"""
    area_int: area integral of a 2D array in grid units, with weights for the boundary rows and columns.

    Boundary points of a pressure grid represent half cells, or quarter cells at the corners.

        Parameters:
            p  (float): (y, x) array
            wefc (float): Weight of the west and east columns
            snfc (float): Weight of the south and north rows; corners get both

        Returns:
            total (float): Weighted sum

        Formula:
            \sum_{\text{interior}} p + s\sum_{\text{S,N rows}} p + w\sum_{\text{W,E columns}} p + w s \sum_{\text{corners}} p

        See also:
            coarse_grain_P
    """
    return (np.sum(p[1:-1, 1:-1]) + snfc*(np.sum(p[0, 1:-1]) + np.sum(p[-1, 1:-1]))
            + wefc*(np.sum(p[1:-1, 0]) + np.sum(p[1:-1, -1]))
            + wefc*snfc*(p[0, 0] + p[-1, 0] + p[0, -1] + p[-1, -1]))


def coarse_grain_P(ph, ratio):
    r"""
    coarse_grain_P: coarse-grain a high-resolution pressure-grid field, keeping its area integral.

    The field is smoothed with a Gaussian of sigma ratio/4 and subsampled every ratio points; a constant
    boundary value then restores the area integral of the original.

        Parameters:
            ph  (float): (y, x) high-resolution field
            ratio (int): Resolution ratio

        Returns:
            pl (float): Coarse field, with one boundary row and column on each side

        Formula:
            p_c = (G_{r/4} * p_h)_{\,r-1::r},\qquad p_{bc} = \frac{\int p_h / r^2 - \int p_c}{\text{boundary points} + 1}

        See also:
            area_int
    """
    ph_m = area_int(ph, 0.5, 0.5) / (ratio**2.)
    pf = scipy.ndimage.gaussian_filter(ph, sigma=0.25*ratio, mode='constant')[ratio-1:-ratio+1:ratio, ratio-1:-ratio+1:ratio]
    pf_m = area_int(pf, 1.0, 1.0)
    pbc = (ph_m - pf_m) / (np.sum(pf.shape) + 1.)
    pl = np.zeros((pf.shape[0]+2, pf.shape[1]+2), dtype=np.float64)
    pl[1:-1, 1:-1] = pf
    pl[0, :] = pl[-1, :] = pl[1:-1, 0] = pl[1:-1, -1] = pbc
    return pl
