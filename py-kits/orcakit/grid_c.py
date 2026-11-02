#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module grid_c

Arakawa C-grid interpolation, differentiation, divergence and vorticity
operators. These are grid-generic: every operator that needs boundary
padding takes an `exchange` callable rather than importing one, so this
module has no dependency on any specific grid family (e.g. ORCA's
tripolar north-fold). Pass `orca.halo` (see orca.py) for ORCA-family grids.

    exchange signature: exchange(field, halo_array, isMeridional=False, point="T") -> haloed_field

The north fold flips the sign of both velocity components, and its mapping
depends on the grid point (T, U, V or F).
"""

import numpy as np


def ke_compute(u, v):
    """
    ke_compute: sum of the squares of the four velocities around each interior T-point.

    Kinetic energy is defined at T-points, but on the C grid the velocities are on the faces. Summing the
    squares of the two U-points and two V-points around each T-point puts it there.

    The result is 4 times the kinetic energy, kept for older callers; diagkit.stats.kinetic_energy
    returns the energy itself. The outer rows and columns are 0.

        Parameters:
            u (float): (..., ny, nx) U-point velocity
            v (float): (..., ny, nx) V-point velocity

        Returns:
            ke (float): (..., ny, nx) u(i-1)^2 + u(i)^2 + v(j-1)^2 + v(j)^2 on interior T-points,
                        0 on the outer rows and columns

        Formula:
            K_{i,j} = u_{i-1,j}^2 + u_{i,j}^2 + v_{i,j-1}^2 + v_{i,j}^2

        See also:
            uv2t, diagkit.stats.kinetic_energy
    """
    ke = np.zeros_like(u)
    up2 = u**2
    vp2 = v**2
    ke[..., 1:-1, 1:-1] = (up2[..., 1:-1, 0:-2] + up2[..., 1:-1, 1:-1] +
                           vp2[..., 0:-2, 1:-1] + vp2[..., 1:-1, 1:-1])
    return ke


def uv2t(u, v):
    """
    uv2t: sum of the four velocities around each interior T-point.

    The companion of ke_compute, with the velocities themselves in place of their squares.

    The outer rows and columns are 0. For a proper average of one component onto T-points, see
    interp_UV2T.

        Parameters:
            u (float): (..., ny, nx) U-point velocity
            v (float): (..., ny, nx) V-point velocity

        Returns:
            t (float): (..., ny, nx) u(i-1) + u(i) + v(j-1) + v(j) on interior T-points, 0 on the
                       outer rows and columns

        Formula:
            S_{i,j} = u_{i-1,j} + u_{i,j} + v_{i,j-1} + v_{i,j}

        See also:
            ke_compute, interp_UV2T
    """
    t = np.zeros_like(u)
    t[..., 1:-1, 1:-1] = (u[..., 1:-1, 0:-2] + u[..., 1:-1, 1:-1] +
                          v[..., 0:-2, 1:-1] + v[..., 1:-1, 1:-1])
    return t


def interp_UV2T(field_in, exchange, weights_in=None, dir=None, isMeridional=False):
    r"""
    interp_UV2T: average a U-point or V-point field onto T-points.

    On the C grid, u sits on the east face of each cell and v on its north face, so the two
    components, and anything computed from both, are never at the same place. Products such as
    u v, kinetic energy or the covariance of u and v need them co-located first: this averages
    the two faces on either side of each T-point, along x for U fields and along y for V fields.

    Along x the western neighbour of column 0 comes from the periodic seam, through exchange.
    Along y, row 0 has no V-point below it and keeps V(0); on the ORCA grids that row is land.
    With weights, the average is weighted, e.g. by face areas e2u e3u to average transports.

        Parameters:
            field_in    (float): (..., ny, nx) field on U-points (dir='lon') or V-points (dir='lat')
            exchange (callable): Boundary/halo exchange, e.g. orca.halo
            weights_in  (float): Weights at the same points, same shape, or None for a plain mean
            dir           (str): 'lon' for U to T, 'lat' for V to T
            isMeridional (bool): Unused: neither direction reads across the north fold

        Returns:
            field_out (float): (..., ny, nx) field on T-points

        Formula:
            \overline{u}^{\,i}_{i,j} = \tfrac12\,(u_{i-1,j} + u_{i,j}),\qquad \overline{v}^{\,j}_{i,j} = \tfrac12\,(v_{i,j-1} + v_{i,j})
            \text{with weights: } \overline{u}^{\,i}_{i,j} = \frac{w_{i-1,j}\,u_{i-1,j} + w_{i,j}\,u_{i,j}}{w_{i-1,j} + w_{i,j}}

        See also:
            interp_T2UV, interp_1D, noisekit.running.running_cov
    """
    if dir == "lon":
        # First, add a halo point at western boundary (since U is on the right of T)
        field_u = exchange(field_in, np.array([1, 0, 0]))
        if weights_in is None:
            field_out = (field_u[..., :, 1:] + field_u[..., :, :-1]) * 0.5
        else:
            weight_t = exchange(weights_in, np.array([1, 0, 0]))
            field_out = (weight_t[..., :, 1:] * field_u[..., :, 1:]
                         + weight_t[..., :, :-1] * field_u[..., :, :-1]) \
                        / (weight_t[..., :, 1:] + weight_t[..., :, :-1])
    elif dir == "lat":
        # V->T needs the BACKWARD neighbour V(p-1); no exchange() call needed
        # since it's already available directly in the un-haloed array for
        # every interior row (the north halo only ever extends forward/at
        # the high-index end, so it can't supply a backward neighbour here
        # -- see orca.halo's docstring). Row 0 (T(0) would need the
        # nonexistent V(-1)) uses a Neumann-style V(-1)~=V(0) fallback.
        if weights_in is None:
            field_out = np.zeros_like(field_in)
            field_out[..., 1:, :] = (field_in[..., 1:, :] + field_in[..., :-1, :]) * 0.5
            field_out[..., 0, :] = field_in[..., 0, :]
        else:
            field_out = np.zeros_like(field_in)
            field_out[..., 1:, :] = (weights_in[..., 1:, :] * field_in[..., 1:, :]
                                      + weights_in[..., :-1, :] * field_in[..., :-1, :]) \
                                     / (weights_in[..., 1:, :] + weights_in[..., :-1, :])
            field_out[..., 0, :] = field_in[..., 0, :]
    else:
        raise ValueError("dir must be 'lon' or 'lat'")
    return field_out


def interp_T2UV(field_in, exchange, weights_in=None, dir=None, isMeridional=False):
    r"""
    interp_T2UV: average a T-point field onto U-points or V-points.

    Tracers, sea level and cell thicknesses live at the centre of the cells, but terms of the momentum or
    transport budgets are needed on the faces, next to the velocities. This averages the two T-points on
    either side of each face, along x for U-points and along y for V-points.

    The neighbour beyond the last column or row comes from the halo: the periodic seam in x, the north
    fold in y, where isMeridional flips the sign for a velocity component. With weights, the average is
    weighted.

        Parameters:
            field_in    (float): Field to be interpolated
            exchange (callable): Boundary/halo exchange, see module docstring
            weights_in  (float): Interpolation weights
            dir           (str): Direction of interpolation, 'lon' or 'lat'
            isMeridional (bool): if True, change the sign of the north halo

        Returns:
            field_out   (float): Resulting field

        Formula:
            \overline{T}^{\,i}_{i,j} = \tfrac12\,(T_{i,j} + T_{i+1,j}),\qquad \overline{T}^{\,j}_{i,j} = \tfrac12\,(T_{i,j} + T_{i,j+1})

        See also:
            interp_UV2T, interp_1D
    """
    if dir == "lon":
        field_t = exchange(field_in, np.array([0, 0, 1]))
        if weights_in is None:
            field_out = (field_t[..., :, 1:] + field_t[..., :, :-1]) * 0.5
        else:
            weight_t = exchange(weights_in, np.array([0, 0, 1]))
            field_out = (weight_t[..., :, 1:] * field_t[..., :, 1:]
                         + weight_t[..., :, :-1] * field_t[..., :, :-1]) \
                        / (weight_t[..., :, 1:] + weight_t[..., :, :-1])
    elif dir == "lat":
        field_t = exchange(field_in, np.array([0, 1, 0]), isMeridional)
        if weights_in is None:
            field_out = (field_t[..., 1:, :] + field_t[..., :-1, :]) * 0.5
        else:
            weight_t = exchange(weights_in, np.array([0, 1, 0]))
            field_out = (weight_t[..., 1:, :] * field_t[..., 1:, :]
                         + weight_t[..., :-1, :] * field_t[..., :-1, :]) \
                        / (weight_t[..., 1:, :] + weight_t[..., :-1, :])
    else:
        raise ValueError("dir must be 'lon' or 'lat'")
    return field_out


def interp_1D(field_in, exchange, weights_in=None, dir=None, isMeridional=False, point="T"):
    r"""
    interp_1D: average a field with its next neighbour along x or along y.

    The general form of interp_T2UV for any C-grid point: each value is averaged with the one after it, so
    a field moves half a cell forward, for example from U-points to F-points along y.

    The neighbour after the last column or row comes from the halo, with the fold rule of `point` and the
    sign flip given by isMeridional. With weights, the average is weighted.

        Parameters:
            field_in    (float): Field to be interpolated
            exchange (callable): Boundary/halo exchange, see module docstring
            weights_in  (float): Interpolation weights
            dir           (str): Direction of interpolation, 'lon' or 'lat'
            isMeridional (bool): if True, change the sign of the north halo
            point         (str): Where field_in lives, 'T', 'U', 'V' or 'F'

        Returns:
            field_out   (float): Resulting field

        Formula:
            \overline{f}^{\,i}_{i,j} = \tfrac12\,(f_{i,j} + f_{i+1,j}),\qquad \overline{f}^{\,j}_{i,j} = \tfrac12\,(f_{i,j} + f_{i,j+1})

        See also:
            interp_T2UV, fdiff_1D
    """
    if dir == "lon":
        field_t = exchange(field_in, np.array([0, 0, 1]))
        if weights_in is None:
            field_out = (field_t[..., :, 1:] + field_t[..., :, :-1]) * 0.5
        else:
            weight_t = exchange(weights_in, np.array([0, 0, 1]))
            field_out = (weight_t[..., :, 1:] * field_t[..., :, 1:]
                         + weight_t[..., :, :-1] * field_t[..., :, :-1]) \
                        / (weight_t[..., :, 1:] + weight_t[..., :, :-1])
    elif dir == "lat":
        field_t = exchange(field_in, np.array([0, 1, 0]), isMeridional, point)
        if weights_in is None:
            field_out = (field_t[..., 1:, :] + field_t[..., :-1, :]) * 0.5
        else:
            weight_t = exchange(weights_in, np.array([0, 1, 0]), False, point)
            field_out = (weight_t[..., 1:, :] * field_t[..., 1:, :]
                         + weight_t[..., :-1, :] * field_t[..., :-1, :]) \
                        / (weight_t[..., 1:, :] + weight_t[..., :-1, :])
    else:
        raise ValueError("dir must be 'lon' or 'lat'")
    return field_out


def fdiff_1D(field_in, weights_in, exchange, dir=None, isMeridional=False, point="T"):
    r"""
    fdiff_1D: forward difference along x or along y, divided by a scale factor.

    Gradients on the C grid are differences between neighbours divided by the distance between them, and
    they land half a cell forward: the x-gradient of a T-point field sits on U-points. The weights are
    that distance, e.g. e1u for d/dx of a T field.

    The neighbour after the last column or row comes from the halo, with the fold rule of `point`.

        Parameters:
            field_in    (float): Field to be differentiated
            weights_in  (float): Metric (e.g. e1/e2 scale factor)
            exchange (callable): Boundary/halo exchange, see module docstring
            dir           (str): Direction of differentiation, 'lon' or 'lat'
            isMeridional (bool): if True, change the sign of the north halo
            point         (str): Where field_in lives, 'T', 'U', 'V' or 'F'

        Returns:
            field_out   (float): Resulting field

        Formula:
            (\delta_i f)_{i,j} = \frac{f_{i+1,j} - f_{i,j}}{e_{i,j}},\qquad (\delta_j f)_{i,j} = \frac{f_{i,j+1} - f_{i,j}}{e_{i,j}}

        See also:
            bdiff_1D, interp_1D
    """
    if dir == "lon":
        field_t = exchange(field_in, np.array([0, 0, 1]))
        weight = exchange(weights_in, np.array([0, 0, 1]))
        field_out = (field_t[..., :, 1:] - field_t[..., :, :-1]) / weight[..., :, :-1]
    elif dir == "lat":
        field_t = exchange(field_in, np.array([0, 1, 0]), isMeridional, point)
        weight = exchange(weights_in, np.array([0, 1, 0]), False, point)
        field_out = (field_t[..., 1:, :] - field_t[..., :-1, :]) / weight[..., :-1, :]
    else:
        raise ValueError("dir must be 'lon' or 'lat'")
    return np.array(field_out)


def bdiff_1D(field_in, weights_in, exchange, dir=None, isMeridional=False):
    r"""
    bdiff_1D: backward difference along x or along y, divided by a scale factor.

    The backward counterpart of fdiff_1D, for gradients that land half a cell back, for example from
    U-points to T-points along x.

    Along x the neighbour before column 0 comes from the periodic seam. Along y row 0 has no neighbour
    below and is set to 0; on the ORCA grids that row is land.

        Parameters:
            field_in    (float): Field to be differentiated
            weights_in  (float): Metric (e.g. e1/e2 scale factor)
            exchange (callable): Boundary/halo exchange, see module docstring
            dir           (str): Direction of differentiation, 'lon' or 'lat'
            isMeridional (bool): if True, change the sign of the north halo

        Returns:
            field_out   (float): Resulting field

        Formula:
            (\delta_i f)_{i,j} = \frac{f_{i,j} - f_{i-1,j}}{e_{i,j}},\qquad (\delta_j f)_{i,j} = \frac{f_{i,j} - f_{i,j-1}}{e_{i,j}}

        See also:
            fdiff_1D
    """
    if dir == "lon":
        field_t = exchange(field_in, np.array([1, 0, 0]))
        weight = exchange(weights_in, np.array([1, 0, 0]))
        field_out = (field_t[..., :, 1:] - field_t[..., :, :-1]) / weight[..., :, 1:]
    elif dir == "lat":
        field_out = np.zeros_like(field_in)
        field_out[..., 1:, :] = (field_in[..., 1:, :] - field_in[..., :-1, :]) / weights_in[..., 1:, :]
    else:
        raise ValueError("dir must be 'lon' or 'lat'")
    return np.array(field_out)


def dxFu_1D(u, e1t, e2u, exchange):
    r"""
    dxFu_1D: the zonal part of a flux divergence on T-points.

    A divergence on the C grid is the difference of the fluxes through the faces of a cell. This is the x
    part alone: the flux e2u u through each east face, differenced with the west face and divided by e1t.

    The west face of column 0 comes from the periodic seam.

        Parameters:
            u           (float): Zonal, X-directed velocity (U-point)
            e1t         (float): Metric term on T-point (output grid)
            e2u         (float): Metric term on U-point (building the flux)
            exchange (callable): Boundary/halo exchange, see module docstring

        Returns:
            field_out   (float): d(Fu)/dx, on T-points

        Note: the parameter names used to be swapped (e1u, e2t) and Fu was
        built as u*e2t -- a T-point metric applied to a U-point velocity
        with no interpolation. Also used an east-halo/forward-difference,
        which is the wrong direction for a T-point divergence (needs the
        backward neighbour, see hDiv_Vgrid). Both fixed here, matching
        hDiv_Vgrid's (verified) Dlon_Fu convention.

        Formula:
            \frac{(e_{2u}\,u)_{i,j} - (e_{2u}\,u)_{i-1,j}}{e_{1t}}

        See also:
            hDiv_Vgrid
    """
    Fu = u * e2u
    Fu = exchange(Fu, np.array([1, 0, 0]))
    Dlon_Fu = (Fu[..., :, 1:] - Fu[..., :, :-1]) / e1t
    return np.array(Dlon_Fu)


def hDiv_Vgrid(u, v, dx_t, dy_t, dx_v, dy_u, exchange, isMeridional=False):
    r"""
    hDiv_Vgrid: horizontal divergence on T-points of velocities on U-points and V-points.

    The divergence measures how much flow leaves a cell: the net flux through its four faces, each the
    velocity times the face width, divided by the cell area. It is the C-grid form NEMO uses, and with
    transports e3u u and e3v v in place of u and v it gives the transport divergence that sets the
    vertical velocity.

    The west face of column 0 comes from the periodic seam. Row 0 has no south face and gets no
    meridional term; on the ORCA grids that row is land.

        Parameters:
            u           (float): Zonal, X-directed velocity
            v           (float): Meridional, Y-directed velocity
            dx_t,dy_t   (float): Metric terms on T-point
            dx_v        (float): Metric terms on V-point
            dy_u        (float): Metric terms on U-point
            exchange (callable): Boundary/halo exchange, see module docstring
            isMeridional (bool): if True, change the sign of the north halo

        Returns:
            field_out   (float): Resulting field

        Usage:
            divH = hDiv_Vgrid(u, v, e1t, e2t, e1v, e2u, orca.halo)

        Formula:
            \nabla_h\cdot\mathbf{u} = \frac{1}{e_{1t}e_{2t}}\Big[(e_{2u}u)_{i,j} - (e_{2u}u)_{i-1,j} + (e_{1v}v)_{i,j} - (e_{1v}v)_{i,j-1}\Big]

        NEMO:
            src/OCE/DYN/divhor.F90: hdiv(ji,jj,jk) = (  (  e2u(ji  ,jj) * e3u(ji  ,jj,jk,Kmm) * pu(ji  ,jj,jk)

        See also:
            vertical_velocity, hDiv_Tgrid
    """
    # T-point divergence needs the BACKWARD neighbour in both directions
    # (U(p)-U(p-1), V(p)-V(p-1)). West halo genuinely gives a backward
    # neighbour (see orca.halo docstring), so Fu uses that; Fv has no such
    # halo available (north only ever extends forward) so it's diffenced
    # directly against the un-haloed array, with row 0 left at 0 (no V(-1)).
    Fu = u * dy_u
    Fv = v * dx_v
    Fu = exchange(Fu, np.array([1, 0, 0]))
    dxdy = dx_t * dy_t
    Dlon_Fu = Fu[..., :, 1:] - Fu[..., :, :-1]
    Dlat_Fv = np.zeros_like(Fv)
    Dlat_Fv[..., 1:, :] = Fv[..., 1:, :] - Fv[..., :-1, :]
    hdiff = (Dlon_Fu + Dlat_Fv) / dxdy
    return np.array(hdiff)


def w_from_divergence(div, e3t, tmask):
    r"""
    w_from_divergence: vertical velocity by integrating the divergence up from the bottom.

    Incompressible flow that converges in a layer must leave it vertically. Starting from w = 0 below the
    deepest level, each level adds the divergence times its thickness, which gives w at the top of every
    T-cell, as NEMO computes it.

    Levels are axis -3, surface first. The divergence is counted only on wet points.

        Parameters:
            div   (float): (..., nz, ny, nx) horizontal divergence at T-points, in s-1
            e3t   (float): (..., nz, ny, nx) T-cell thicknesses, in metres, or broadcastable
            tmask (float): (..., nz, ny, nx) T-point wet mask

        Returns:
            w (float): (..., nz, ny, nx) vertical velocity at the top of T-cells, in m/s

        Formula:
            w_k = w_{k+1} - e_{3t,k}\,(\nabla_h\cdot\mathbf{u})_k\,t_k,\qquad w = 0 \text{ below the deepest level}

        NEMO:
            src/OCE/DYN/sshwzv.F90: pww(ji,jj,jk) = pww(ji,jj,jk+1) - (  e3t(ji,jj,jk,Kmm) * hdiv(ji,jj,jk)  ) * tmask(ji,jj,jk)

        See also:
            vertical_velocity, hDiv_Vgrid
    """
    w = np.zeros_like(div)
    for k in range(div.shape[-3]-2, -1, -1):
        w[..., k, :, :] = w[..., k+1, :, :] - e3t[..., k, :, :] * div[..., k, :, :] * tmask[..., k, :, :]
    return w


def vertical_velocity(u, v, e1t, e2t, e1v, e2u, e3u, e3v, umask, vmask, tmask, exchange):
    r"""
    vertical_velocity: vertical velocity from horizontal velocity, by continuity, as NEMO computes it.

    Many products store u and v but not w, though w follows from them: what converges horizontally in a
    layer must leave it vertically. Here the fluxes through each face carry the face's thickness, so the
    partial cells at the sea floor are counted with their true size and no T-cell thickness is needed.
    On a NEMO grid this reproduces the model's own w.

    w sits at the top of each T-cell, zero below the deepest level. Velocities on land count as 0, and
    the result is NaN on land.

        Parameters:
            u      (float): (nz, ny, nx) U-point velocity, may hold NaN on land
            v      (float): (nz, ny, nx) V-point velocity, may hold NaN on land
            e1t    (float): (ny, nx) T-cell widths along x, in metres
            e2t    (float): (ny, nx) T-cell widths along y, in metres
            e1v    (float): (ny, nx) V-cell widths along x, in metres
            e2u    (float): (ny, nx) U-cell widths along y, in metres
            e3u    (float): (nz, ny, nx) U-cell thicknesses, in metres
            e3v    (float): (nz, ny, nx) V-cell thicknesses, in metres
            umask  (array): (nz, ny, nx) U-point wet mask
            vmask  (array): (nz, ny, nx) V-point wet mask
            tmask  (array): (nz, ny, nx) T-point wet mask
            exchange (callable): Boundary/halo exchange, e.g. orca.halo

        Returns:
            w (float): (nz, ny, nx) vertical velocity at the top of T-cells, in m/s; NaN on land

        Formula:
            w_k = w_{k+1} - \frac{(e_{2u}e_{3u}u)_{i,j} - (e_{2u}e_{3u}u)_{i-1,j} + (e_{1v}e_{3v}v)_{i,j} - (e_{1v}e_{3v}v)_{i,j-1}}{e_{1t}e_{2t}}\bigg|_k,\qquad w = 0 \text{ at the bottom}

        NEMO:
            src/OCE/DYN/divhor.F90: hdiv(ji,jj,jk) = (  (  e2u(ji  ,jj) * e3u(ji  ,jj,jk,Kmm) * pu(ji  ,jj,jk)
            src/OCE/DYN/sshwzv.F90: pww(ji,jj,jk) = pww(ji,jj,jk+1) - (  e3t(ji,jj,jk,Kmm) * hdiv(ji,jj,jk)  ) * tmask(ji,jj,jk)

        See also:
            w_from_divergence, hDiv_Vgrid, orca.e3uv_from_partial_steps
    """
    u, v = np.nan_to_num(u)*umask, np.nan_to_num(v)*vmask
    transport = hDiv_Vgrid(u*e3u, v*e3v, e1t, e2t, e1v, e2u, exchange)
    return np.where(tmask, w_from_divergence(transport, np.ones((u.shape[0], 1, 1)), tmask), np.nan)


def geostrophic_uv(ssh, tmask, umask, vmask, e1t, e2t, e1v, e2u, gphiu, gphiv, exchange, min_lat=5.0, gravity=9.80665, omega=7.292116e-5):
    r"""
    geostrophic_uv: surface geostrophic velocity from sea surface height.

    Away from the equator and the coasts, the large-scale surface flow is close to geostrophic balance:
    it follows the lines of constant sea level, faster where they are closer. This gives that flow on the
    full grid, so it can be compared with the model's velocity or used where only sea level is known.

    SSH is averaged to F-points over wet T-points only, then differenced onto U-points and V-points, with
    the Coriolis parameter at each of them. Within about 5 degrees of the equator geostrophy does not hold
    and the result is NaN; within about two cells of land the wet-only average biases it.

        Parameters:
            ssh    (float): (ny, nx) sea surface height, in metres; NaN counts as land
            tmask   (bool): (ny, nx) T-point wet mask
            umask   (bool): (ny, nx) U-point wet mask
            vmask   (bool): (ny, nx) V-point wet mask
            e1t    (float): (ny, nx) T-cell widths along x, in metres
            e2t    (float): (ny, nx) T-cell widths along y, in metres
            e1v    (float): (ny, nx) V-cell widths along x, in metres
            e2u    (float): (ny, nx) U-cell widths along y, in metres
            gphiu  (float): (ny, nx) latitude of U-points, in degrees
            gphiv  (float): (ny, nx) latitude of V-points, in degrees
            exchange (callable): Boundary/halo exchange, e.g. orca.halo
            min_lat (float): Latitude, in degrees, below which the result is NaN
            gravity (float): Gravity, in m s-2; NEMO's value by default
            omega   (float): Earth's rotation rate, in s-1

        Returns:
            Ug (float): (ny, nx) zonal velocity on U-points, in m/s; NaN on land, for |lat| < min_lat
                        and on row 0
            Vg (float): (ny, nx) meridional velocity on V-points, in m/s; NaN on land and for
                        |lat| < min_lat

        Formula:
            \eta_f = \frac{\sum_{4\,T} A\,\eta}{\sum_{4\,T} A}\ \text{(wet T-points only)},\qquad A = e_{1t}e_{2t}
            U_g = -\frac{g}{f_u}\,\frac{\eta_f(i,j) - \eta_f(i,j-1)}{e_{2u}},\qquad V_g = \frac{g}{f_v}\,\frac{\eta_f(i,j) - \eta_f(i-1,j)}{e_{1v}},\qquad f = 2\Omega\sin\varphi

        See also:
            noisekit.noise.uv_from_ssh_dA
    """
    wet = tmask & np.isfinite(ssh)
    corners = np.array([1, 1, 1])
    area = exchange(np.where(wet, e1t*e2t, 0.0), corners)
    weighted = exchange(np.where(wet, ssh*e1t*e2t, 0.0), corners)

    def corner_sum(a):
        """
        corner_sum: sum of each 2 x 2 block, the four T-points around an F-point.

        A helper of geostrophic_uv, for the area-weighted average of SSH onto F-points.

            Parameters:
                a (float): (ny+1, nx+2) padded array

            Returns:
                s (float): (ny, nx+1) sums

            Formula:
                s_{j,i} = a_{j,i} + a_{j+1,i} + a_{j+1,i+1} + a_{j,i+1}

            See also:
                geostrophic_uv
        """
        return a[:-1, :-1] + a[1:, :-1] + a[1:, 1:] + a[:-1, 1:]

    with np.errstate(invalid="ignore", divide="ignore"):
        # Column c holds the F-points of column c - 1, from -1 to nx - 1
        ssh_f = corner_sum(weighted)/corner_sum(area)
        Ug = np.full(ssh.shape, np.nan)
        Ug[1:] = -gravity*(ssh_f[1:, 1:] - ssh_f[:-1, 1:])/(e2u[1:]*2*omega*np.sin(np.radians(gphiu[1:])))
        Vg = gravity*(ssh_f[:, 1:] - ssh_f[:, :-1])/(e1v*2*omega*np.sin(np.radians(gphiv)))
    return (np.where(umask & (np.abs(gphiu) >= min_lat), Ug, np.nan),
            np.where(vmask & (np.abs(gphiv) >= min_lat), Vg, np.nan))


def hDiv_Tgrid(u, v, e1t, e2t, e1v, e2u, exchange, isMeridional=True):
    r"""
    hDiv_Tgrid: horizontal divergence on T-points of velocities that are themselves on T-points.

    For velocities already averaged or interpolated to the cell centres, as some datasets store them.
    They are first averaged onto the faces, then their fluxes are differenced as in hDiv_Vgrid.

    The halo supplies the neighbours beyond the last column and row, with the sign flip at the fold for
    v. Row 0 gets no meridional term.

        Parameters:
            u           (float): Zonal, X-directed velocity
            v           (float): Meridional, Y-directed velocity
            e1t, e2t    (float): Metric terms on T-point
            e1v         (float): Metric terms on V-point
            e2u         (float): Metric terms on U-point
            exchange (callable): Boundary/halo exchange, see module docstring
            isMeridional (bool): if True, change the sign of v's north halo

        Returns:
            field_out   (float): Resulting field. Row 0 has no meridional term.

        Usage:
            divH = hDiv_Tgrid(u, v, e1t, e2t, e1v, e2u, orca.halo)

        Formula:
            \nabla_h\cdot\mathbf{u} = \frac{1}{e_{1t}e_{2t}}\Big[(e_{2u}\overline{u}^{\,i})_{i,j} - (e_{2u}\overline{u}^{\,i})_{i-1,j} + (e_{1v}\overline{v}^{\,j})_{i,j} - (e_{1v}\overline{v}^{\,j})_{i,j-1}\Big]

        See also:
            hDiv_Vgrid, interp_T2UV
    """
    # Stage 1: interpolate T-point u,v to U,V-points (forward pairing).
    u = exchange(u, np.array([0, 0, 1]))
    v = exchange(v, np.array([0, 1, 0]), isMeridional)

    Fu = (u[..., :, 1:] + u[..., :, :-1]) * e2u / 2
    Fv = (v[..., :-1, :] + v[..., 1:, :]) * e1v / 2

    # Stage 2: T-point divergence needs the BACKWARD neighbour of the U/V-point
    # fluxes. The west halo gives Fu's; Fv has none, so row 0 is left at 0.
    Fu = exchange(Fu, np.array([1, 0, 0]))
    Dlon_Fu = Fu[..., :, 1:] - Fu[..., :, :-1]
    Dlat_Fv = np.zeros_like(Fv)
    Dlat_Fv[..., 1:, :] = Fv[..., 1:, :] - Fv[..., :-1, :]

    hdiff = (Dlon_Fu + Dlat_Fv) / (e1t * e2t)
    return np.array(hdiff)


def vort_Vgrid(u, v, dx_f, dy_f, dx_u, dy_v, exchange):
    r"""
    vort_Vgrid: relative vorticity on F-points of velocities on U-points and V-points.

    Vorticity measures the local rotation of the flow, positive for anticlockwise rotation seen
    from above. On the C grid it sits at the corners of the cells, from the circulation of u and v
    around each corner divided by its area, as in NEMO's vorticity term.

    The neighbours beyond the last column and row come from the halo, with the sign flip of both
    velocity components at the fold.

        Parameters:
            u           (float): Zonal, X-directed velocity
            v           (float): Meridional, Y-directed velocity
            dx_f,dy_f   (float): Metric terms on F-point, where the vorticity lives
            dx_u        (float): Metric terms on U-point
            dy_v        (float): Metric terms on V-point
            exchange (callable): Boundary/halo exchange, see module docstring

        Returns:
            field_out   (float): Resulting field, on F-points

        Usage:
            vort = vort_Vgrid(u, v, e1f, e2f, e1u, e2v, orca.halo)

        Formula:
            \zeta_{i,j} = \frac{1}{e_{1f}e_{2f}}\Big[(e_{2v}v)_{i+1,j} - (e_{2v}v)_{i,j} - (e_{1u}u)_{i,j+1} + (e_{1u}u)_{i,j}\Big]

        NEMO:
            src/OCE/DYN/dynvor.F90: zwz(ji,jj) = (  ( e2v(ji+1,jj) * pv(ji+1,jj,jk) - e2v(ji,jj) * pv(ji,jj,jk) )

        See also:
            hDiv_Vgrid, diagkit.flow.vorticity
    """
    Fu = u * dx_u
    Fv = v * dy_v
    Fu = exchange(Fu, np.array([0, 1, 0]), True, "U")
    Fv = exchange(Fv, np.array([0, 0, 1]), True, "V")
    Dlon_Fv = Fv[..., :, 1:] - Fv[..., :, :-1]
    # Forward difference in j, to match Dlon_Fv's forward difference in i --
    # verified against synthetic shear flows with known analytic vorticity
    # (pure du/dy and solid-body rotation cases); the previous backward
    # difference here made vort = dv/dx + du/dy instead of dv/dx - du/dy.
    Dlat_Fu = Fu[..., 1:, :] - Fu[..., :-1, :]
    vort = (Dlon_Fv - Dlat_Fu) / (dx_f * dy_f)
    return np.array(vort)
