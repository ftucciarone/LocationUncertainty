#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module flow

Strain, vorticity and Okubo-Weiss on a C-grid with uniform spacing (dx, dy
scalars) and no boundary conditions; edge rows and columns stay zero. On ORCA
grids this is only a local approximation: use a box away from the north fold.
Ported from grad_py of the earlier analysis code.
"""
import numpy as np


def _t_from_f(F):
    r"""
    _t_from_f: average an F-point field onto interior T-points.

        Parameters:
            F (float): (..., ny, nx) F-point field

        Returns:
            T (float): (..., ny, nx) mean of the four F-points around each interior T-point,
                       0 on the outer rows and columns

        Formula:
            T_{j,i} = \tfrac14\big(F_{j-1,i-1} + F_{j,i-1} + F_{j,i} + F_{j-1,i}\big)

        See also:
            okubo_weiss
    """
    out = np.zeros_like(F)
    out[..., 1:-1, 1:-1] = (F[..., 0:-2, 0:-2] + F[..., 1:-1, 0:-2] +
                            F[..., 1:-1, 1:-1] + F[..., 0:-2, 1:-1]) / 4.0
    return out


def strain_normal(u, v, dx, dy):
    r"""
    strain_normal: normal strain on interior T-points, on a uniform grid.

    The normal strain measures stretching along x against y. With the shear strain and the vorticity it
    separates straining flow from rotating flow.

        Parameters:
            u  (float): (..., ny, nx) U-point velocity
            v  (float): (..., ny, nx) V-point velocity
            dx (float): Grid spacing in x
            dy (float): Grid spacing in y

        Returns:
            sn (float): (..., ny, nx), 0 on the outer rows and columns

        Formula:
            s_n = \frac{u_{j,i} - u_{j,i-1}}{\Delta x} - \frac{v_{j,i} - v_{j-1,i}}{\Delta y}

        See also:
            strain_shear, okubo_weiss
    """
    out = np.zeros_like(u)
    out[..., 1:-1, 1:-1] = ((u[..., 1:-1, 1:-1] - u[..., 1:-1, 0:-2]) / dx -
                            (v[..., 1:-1, 1:-1] - v[..., 0:-2, 1:-1]) / dy)
    return out


def strain_shear(u, v, dx, dy):
    r"""
    strain_shear: shear strain on F-points, on a uniform grid.

        Parameters:
            u  (float): (..., ny, nx) U-point velocity
            v  (float): (..., ny, nx) V-point velocity
            dx (float): Grid spacing in x
            dy (float): Grid spacing in y

        Returns:
            ss (float): (..., ny, nx), 0 on the last two rows and columns

        Formula:
            s_s = \frac{v_{j,i+1} - v_{j,i}}{\Delta x} + \frac{u_{j+1,i} - u_{j,i}}{\Delta y}

        See also:
            strain_normal, okubo_weiss
    """
    out = np.zeros_like(u)
    out[..., 0:-2, 0:-2] = ((v[..., 0:-2, 1:-1] - v[..., 0:-2, 0:-2]) / dx +
                            (u[..., 1:-1, 0:-2] - u[..., 0:-2, 0:-2]) / dy)
    return out


def vorticity(u, v, dx, dy):
    r"""
    vorticity: relative vorticity on F-points, on a uniform grid.

    The local rotation of the flow, positive anticlockwise seen from above. For the ORCA grids, with
    varying cell sizes, use orcakit.grid_c.vort_Vgrid.

        Parameters:
            u  (float): (..., ny, nx) U-point velocity
            v  (float): (..., ny, nx) V-point velocity
            dx (float): Grid spacing in x
            dy (float): Grid spacing in y

        Returns:
            zeta (float): (..., ny, nx), 0 on the last two rows and columns

        Formula:
            \zeta = \frac{v_{j,i+1} - v_{j,i}}{\Delta x} - \frac{u_{j+1,i} - u_{j,i}}{\Delta y}

        See also:
            orcakit.grid_c.vort_Vgrid, okubo_weiss
    """
    out = np.zeros_like(u)
    out[..., 0:-2, 0:-2] = ((v[..., 0:-2, 1:-1] - v[..., 0:-2, 0:-2]) / dx -
                            (u[..., 1:-1, 0:-2] - u[..., 0:-2, 0:-2]) / dy)
    return out


def okubo_weiss(u, v, dx, dy):
    r"""
    okubo_weiss: the Okubo-Weiss parameter on T-points.

    It separates eddies from the straining flow around them: negative where rotation dominates, inside
    eddy cores, and positive where strain dominates.

        Parameters:
            u  (float): (..., ny, nx) U-point velocity
            v  (float): (..., ny, nx) V-point velocity
            dx (float): Grid spacing in x
            dy (float): Grid spacing in y

        Returns:
            W (float): (..., ny, nx) on T-points

        Formula:
            W = s_n^2 + s_s^2 - \zeta^2

        See also:
            strain_normal, strain_shear, vorticity
    """
    return (strain_normal(u, v, dx, dy)**2 + _t_from_f(strain_shear(u, v, dx, dy))**2
            - _t_from_f(vorticity(u, v, dx, dy))**2)
