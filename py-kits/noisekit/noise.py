#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module noise

Noise field construction. Mode handling is in modredkit.modes.
"""


def uv_from_ssh_dA(ssh, e1t, e2t, ff_t, dx, dy):
    r"""
    uv_from_ssh_dA: geostrophic velocity from sea surface height, on interior points only.

    The earlier form of orcakit.grid_c.geostrophic_uv, which supersedes it: that one keeps the full grid,
    leaves land and the equator out, and obeys the seam and fold.

    SSH is averaged to F-points with area weights, then differenced onto U and V points, with g = 9.806
    and f averaged from T-points. The result covers rows and columns 1 to n-2.

        Parameters:
            ssh  (float): (..., ny, nx) sea surface height, in metres
            e1t  (float): (ny, nx) T-cell widths along x, in metres
            e2t  (float): (ny, nx) T-cell widths along y, in metres
            ff_t (float): (ny, nx) Coriolis parameter at T-points, in s-1
            dx   (float): (ny, nx) e1v
            dy   (float): (ny, nx) e2u

        Returns:
            Ug (float): (..., ny-2, nx-2) zonal velocity on rows and columns 1 to n-2, in m/s
            Vg (float): (..., ny-2, nx-2) meridional velocity, same points, in m/s

        Formula:
            U_g = -\frac{g}{f_u}\,\frac{\delta_j \eta_f}{e_{2u}},\qquad V_g = \frac{g}{f_v}\,\frac{\delta_i \eta_f}{e_{1v}},\qquad \eta_f = \frac{\sum_{4} A\,\eta}{\sum_{4} A}

        See also:
            orcakit.grid_c.geostrophic_uv
    """
    # Compute area metric
    dA = e1t * e2t
    # Interpolate ssh on f
    ssh_f = ( ssh[...,0:-1,0:-1]*dA[0:-1,0:-1] +
              ssh[...,1:  ,0:-1]*dA[1:  ,0:-1] +
              ssh[...,1:  ,1:  ]*dA[1:  ,1:  ] +
              ssh[...,0:-1,1:  ]*dA[0:-1,1:  ] ) / ( dA[0:-1,0:-1] +
              dA[1:  ,0:-1] +
              dA[1:  ,1:  ] +
              dA[0:-1,1:  ] )
    # Interpolate f on u and v
    ff_u = ( ff_t[1:-1,1:-1] + ff_t[1:-1,2:] ) * 0.5
    ff_v = ( ff_t[1:-1,1:-1] + ff_t[2:,1:-1] ) * 0.5
    Ug = - 9.806 * ( ssh_f[...,1:,1:] - ssh_f[...,0:-1,1:] ) / (dy[1:-1,1:-1] * ff_u)
    Vg = + 9.806 * ( ssh_f[...,1:,1:] - ssh_f[...,1:,0:-1] ) / (dx[1:-1,1:-1] * ff_v)
    return Ug, Vg
