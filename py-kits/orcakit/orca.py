#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module orca

ORCA/eORCA-family grid specifics: halo (boundary) exchange including the
tripolar north-fold sign flip, 2D smoothing built on that exchange, nearest
-neighbour grid search, downsampling, and NEMO mask construction.

`halo`/`dehalo` are the `exchange` callable expected by grid_c.py's
operators (interp_UV2T, hDiv_Vgrid, ...).
"""

import os

import numpy as np

whalo = np.array([1, 0, 0])
nhalo = np.array([0, 1, 0])
ehalo = np.array([0, 0, 1])


# North-fold rule on a T-pivot grid (eORCA025 family), checked on the GLORYS04 and
# CMEMS meshes: x(ny-1, i) = x(ny-1-(row-1), nx+col-i), for each grid point.
#   point: (row, col), ghost row ny+p is copied from row ny-row-p
FOLD = {"T": (4, 0),
        "U": (4, -1),
        "V": (5, 0),
        "F": (5, -1)}


def wrap_columns(i, nx):
    r"""
    wrap_columns: bring column indices outside the grid back onto it.

    The last two columns of the ORCA grids repeat the first two, so the grid is periodic with period
    nx - 2, not nx.

    Indices already on the grid are unchanged.

        Parameters:
            i  (int): Column indices, any shape
            nx (int): Number of columns

        Returns:
            i_out (int): Indices in [0, nx); those already inside are unchanged

        Formula:
            i \mapsto i \bmod (n_x - 2)\quad\text{for } i < 0 \text{ or } i \ge n_x

        See also:
            fold_columns, halo
    """
    i = np.asarray(i)
    return np.where((i < 0) | (i >= nx), i % (nx - 2), i)


def fold_columns(i, nx, point="T"):
    r"""
    fold_columns: the column that the north fold maps each column onto.

    At the north fold the grid is mirrored about a pivot, so a column near one end continues as a column
    near the other. The mirror column depends on the point type: T and V points, U and F points.

        Parameters:
            i     (int): Column indices, any shape
            nx    (int): Number of columns
            point (str): 'T', 'U', 'V' or 'F'; the mapping differs by point, see FOLD

        Returns:
            i_out (int): Mirror columns, in [0, nx - 2)

        Formula:
            i \mapsto (n_x + c_p - i) \bmod (n_x - 2),\qquad c_p \text{ from FOLD: } 0 \text{ for T and V, } -1 \text{ for U and F}

        NEMO:
            src/OCE/LBC/lbc_nfd_generic.h90: CASE ( 'T' , 'W' )

        See also:
            halo, wrap_columns
    """
    return (nx + FOLD[point][1] - np.asarray(i)) % (nx - 2)


def halo(eOrca_field, halo_array=np.array([0, 0, 0]), isMeridional=False, point="T"):
    r"""
    halo: add a halo to a field on the ORCA grid, from its own seam and fold.

    Operations that look at neighbours, such as filters and differences, need values beyond the edges of
    the array. On the ORCA grids those are known: east and west, the grid wraps with period nx - 2; north,
    it folds back on itself about a T-point pivot, as in NEMO, with the sign flipped for a velocity
    component. This returns the field with that many extra columns and rows filled.

    The halo is added at the west, the east and the north; there is none to the south, which is closed.
    It holds for the fine grids with overlap columns and a T pivot (GLORYS04, GLORYS12, CMEMS), not for
    eORCA1's domain file, which has neither. The north halo cannot be wider than the rows it mirrors.

        Parameters:
            eOrca_field (float): Field to be haloed, real or complex
            halo_array    (int): np.array([west, north, east])
            isMeridional (bool): if True, change the sign of the north halo. Needed
                                 by BOTH velocity components, not by scalars or variances
            point         (str): Where the field lives: 'T', 'U', 'V' or 'F'

        Returns:
            eOrca_haloed  (float): Resulting field, at least double precision

        Usage:
            halo(t, halo_array=np.array([0,0,1]))
            halo(u, halo_array=np.array([1,2,1]), isMeridional=True, point="U")
            halo(v, halo_array=np.array([0,2,0]), isMeridional=True, point="V")

        Formula:
            x_{i,\,n_y+p} = s\;x_{\mathrm{fold}(i),\;n_y-r-p},\qquad s = -1 \text{ for velocity components, } +1 \text{ otherwise},\qquad r \text{ from FOLD}
            x_{i,j} = x_{i \bmod (n_x-2),\,j} \text{ east and west}

        NEMO:
            src/OCE/LBC/lbc_nfd_generic.h90: ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)

        See also:
            dehalo, fold_columns, mpi_orca.add_halo
    """
    hwest, hnorth, heast = (int(h) for h in halo_array)
    field = np.asarray(eOrca_field)
    ny, nx = field.shape[-2:]
    row, _ = FOLD[point]
    if hnorth > ny - row + 1:
        raise ValueError(f"halo: north halo {hnorth} is too wide for {ny} rows")

    columns = np.arange(-hwest, nx + heast)
    eOrca_haloed = np.zeros(field.shape[:-2] + (ny + hnorth, columns.size), dtype=np.result_type(field, np.float64))
    eOrca_haloed[..., :ny, :] = field[..., wrap_columns(columns, nx)]
    if hnorth:
        rows = ny - row - np.arange(hnorth)
        sign = -1 if isMeridional else 1
        eOrca_haloed[..., ny:, :] = sign * field[..., rows, :][..., fold_columns(columns, nx, point)]
    return eOrca_haloed


def dehalo(eOrca_haloed, halo_array=np.array([0, 0, 0])):
    r"""
    dehalo: remove the halo that halo added.

    Returns the interior after an operation on the padded field.

        Parameters:
            eOrca_haloed (float): Field to be dehaloed
            halo_array     (int): np.array([west, north, east])

        Returns:
            eOrca_dehaloed (float): Resulting field

        Usage:
            dehalo(u, halo_array=np.array([1,2,1]))
            dehalo(v, halo_array=np.array([0,2,0]))
            dehalo(t, halo_array=np.array([0,0,1]))

        Formula:
            y_{j,i} = x_{j,\ i + h_w},\qquad 0 \le j < n_y - h_n,\quad 0 \le i < n_x - h_w - h_e

        See also:
            halo
    """
    # Map indexes
    iwest = 0
    ieast = 2
    inorth = 1

    # Map values
    hwest = halo_array[iwest]
    heast = halo_array[ieast]
    hnorth = halo_array[inorth]
    ny = np.shape(eOrca_haloed)[-2] - hnorth

    # Safety measure in case of zero halo
    if heast == 0:
        heast = -np.shape(eOrca_haloed)[-1] - hwest

    # Return the interior domain
    return eOrca_haloed[..., :ny, hwest:-heast]


def smooth_2D(field_in, sigma, mask=None, isMeridional=False, point="T"):
    r"""
    smooth_2D: Gaussian smoothing in the horizontal, across the ORCA seam and fold.

    Smoothing an ORCA field as a plain array would treat the seam and the north edge as walls. Here the
    field is first padded with its own continuation (halo), smoothed, then cut back, so the filter sees
    the grid as it really connects.

    It does not know about land: land values are averaged in. For a land-aware low-pass use lowpass or
    masked_smooth.

        Parameters:
            field_in    (float): Field to be smoothed
            sigma       (float): Gaussian filter sigma, in gridpoints
            mask        (float): Optional mask applied to the result
            isMeridional (bool): if True, change the sign of the north halo
            point         (str): Where the field lives, see halo

        Returns:
            field_out   (float): Resulting field

        Formula:
            (G_\sigma * f)_{i,j} = \sum_{|m|,|n| \le 4\sigma} \frac{e^{-(m^2+n^2)/2\sigma^2}}{\sum e^{-(m^2+n^2)/2\sigma^2}}\, f_{i+m,\,j+n}

        See also:
            filters.lowpass, filters.masked_smooth, halo
    """
    from scipy.ndimage import gaussian_filter
    # Halo width must cover gaussian_filter's actual kernel reach
    # (truncate=4.0 by default, so ~4*sigma) -- a fixed width of 5 silently
    # falls back to scipy's own generic edge-reflection past that point for
    # any sigma bigger than ~1.25, which produced wildly wrong values near
    # the north fold at low ocean-fraction rows for sigma=6 (smoothed-mask
    # denominator ~1e-5 in filters.masked_smooth, amplifying noise ~300x).
    hwidth = int(np.ceil(4 * sigma))
    halo_array = np.array([hwidth, hwidth, hwidth])
    sig_array = np.zeros(np.ndim(field_in))
    sig_array[-1] = sigma
    sig_array[-2] = sigma
    field = halo(field_in, halo_array=halo_array, isMeridional=isMeridional, point=point)
    if np.iscomplexobj(field):
        field = gaussian_filter(field.real, sigma=sig_array) + 1j * gaussian_filter(field.imag, sigma=sig_array)
    else:
        field = gaussian_filter(field, sigma=sig_array)
    field = dehalo(field, halo_array=halo_array)
    if mask is not None:
        field = field * mask
    return field


def _unit_sphere(lat, lon):
    r"""
    _unit_sphere: points on the unit sphere.

    Distances between Cartesian points need no special case for the longitude wrap or the poles.

        Parameters:
            lat (float): Latitudes, in degrees, any shape
            lon (float): Longitudes, in degrees, same shape

        Returns:
            xyz (float): (..., 3) Cartesian coordinates

        Formula:
            \mathbf{x} = (\cos\varphi\cos\lambda,\ \cos\varphi\sin\lambda,\ \sin\varphi)

        See also:
            find_nearest_grid
    """
    la, lo = np.deg2rad(np.asarray(lat, dtype=np.float64)), np.deg2rad(np.asarray(lon, dtype=np.float64))
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=-1)


def find_nearest_LatLon(lat_array, lon_array, lat_value, lon_value):
    r"""
    find_nearest_LatLon: the grid point nearest to one position, on the sphere.

    For locating a station, a mooring or a box corner on a curvilinear grid, where latitude and longitude
    are not along the axes.

        Parameters:
           lat_array   (float): 2D latitude array in which to search
           lon_array   (float): 2D longitude array in which to search
           lat_value   (float): Latitude to look for in the array
           lon_value   (float): Longitude to look for in the array

        Returns:
          lat_idx_out   (int): Latitude index of the closest gridpoint
          lon_idx_out   (int): Longitude index of the closest gridpoint

        Formula:
            (j, i) = \arg\min_{j,i} \|\mathbf{x}(\varphi_{j,i}, \lambda_{j,i}) - \mathbf{x}(\varphi_0, \lambda_0)\|

        See also:
            find_nearest_grid, _unit_sphere
    """
    distance = np.sum((_unit_sphere(lat_array, lon_array) - _unit_sphere(lat_value, lon_value))**2, axis=-1)
    return np.unravel_index(np.argmin(distance), distance.shape)


def find_nearest_grid(lat_coarse, lon_coarse, lat_fine_in, lon_fine_in, mask=None):
    r"""
    find_nearest_grid: for each point of a coarse grid, the nearest point of a fine grid.

    Sampling a high-resolution field onto a coarser grid needs, for each coarse point, the fine point
    closest to it. Distances are taken on the sphere, so the seam and the fold need no special case.

    The last two fine rows are left out, since they only repeat rows south of the fold. Coarse points
    outside the mask get -999.

        Parameters:
            lat_coarse   (float): Latitude array of the coarse grid
            lon_coarse   (float): Longitude array of the coarse grid
           lat_fine_in   (float): Latitude array of the fine grid
           lon_fine_in   (float): Longitude array of the fine grid
                  mask    (bool): Restrict the search to these coarse-grid points

        Returns:
          lat_fine_out   (int): Latitude index of the closest fine-gridpoint to each coarse gridpoint, -999 off the mask
          lon_fine_out   (int): Longitude index of the closest fine-gridpoint to each coarse gridpoint, -999 off the mask

        Callees:
           scipy.spatial.cKDTree

        Formula:
            \mathbf{x} = (\cos\varphi\cos\lambda,\ \cos\varphi\sin\lambda,\ \sin\varphi),\qquad \text{nearest} = \arg\min_f \|\mathbf{x}_c - \mathbf{x}_f\|

        See also:
            load_or_make_map, downsample_field, assign_to_coarse_grid
    """
    from scipy.spatial import cKDTree
    lat_fine = np.asarray(lat_fine_in)[:-2]
    lat_fine_out = np.zeros_like(lat_coarse, dtype=np.int16) - 999
    lon_fine_out = np.zeros_like(lat_coarse, dtype=np.int16) - 999
    wanted = np.ones(np.shape(lat_coarse), dtype=bool) if mask is None else np.asarray(mask, dtype=bool)
    tree = cKDTree(_unit_sphere(lat_fine, np.asarray(lon_fine_in)[:-2]).reshape(-1, 3))
    _, nearest = tree.query(_unit_sphere(lat_coarse, lon_coarse)[wanted])
    lat_fine_out[wanted], lon_fine_out[wanted] = np.unravel_index(nearest, lat_fine.shape)
    return lat_fine_out, lon_fine_out


def downsample_field(field_HR, ilat_HR, ilon_HR, mask_LR, nz=None):
    r"""
    downsample_field: sample a fine field onto a coarse grid at the nearest fine points.

    Each coarse point takes the value of its nearest fine point, from find_nearest_grid. It keeps the
    fine values unaveraged, so a noisy field stays noisy; coarse_grain averages instead.

    Points outside the coarse mask are 0. Leading axes, such as levels, are kept.

        Parameters:
            field_HR   (float): High-resolution field to sample from
            ilat_HR      (int): Latitude index array, from find_nearest_grid
            ilon_HR      (int): Longitude index array, from find_nearest_grid
            mask_LR     (bool): Low-resolution mask; only True points are filled
            nz           (int): Unused, kept for older callers; leading axes follow field_HR

        Returns:
            field_ds   (float): Downsampled field, field_HR's leading axes + mask_LR's shape

        Formula:
            f^{c}_{j,i} = f^{f}_{\,\mathrm{ilat}_{j,i},\ \mathrm{ilon}_{j,i}} \ \text{where } m^{c}_{j,i},\ 0 \text{ elsewhere}

        See also:
            find_nearest_grid, load_or_make_map, coarse_grain
    """
    mask_LR = np.asarray(mask_LR, dtype=bool)
    field_HR = np.asarray(field_HR)
    field_ds = np.zeros(field_HR.shape[:-2] + mask_LR.shape, dtype=field_HR.dtype)
    field_ds[..., mask_LR] = field_HR[..., ilat_HR[mask_LR], ilon_HR[mask_LR]]
    return field_ds


def create_tmask_from_top_and_bottom(top, bot, nz, nlat, nlon):
    r"""
    create_tmask_from_top_and_bottom: the 3D T-point mask from the first and last wet levels of each column.

    Domain files store only those two levels per column; the mask is 1 between them and 0 elsewhere, as
    NEMO builds it.

        Parameters:
                top   (float): "top_level" array from domain_cfg.nc
                bot   (float): "bottom_level" array from domain_cfg.nc
                 nz   (float): Vertical dimension
               nlat   (float): Latitudinal dimension
               nlon    (bool): Longitudinal dimension

        Returns:
          tmask_out   (float): mask at t-point

        Last check date:  17/12/2025

        Formula:
            t_{i,j,k} = 1 \text{ if } \mathrm{top}_{i,j} \le k \le \mathrm{bottom}_{i,j},\ 0 \text{ otherwise}

        See also:
            create_umask_from_tmask, create_vmask_from_tmask
    """
    tmask_out = np.zeros([nz, nlat, nlon])
    for k in range(1, np.shape(tmask_out)[0] + 1):
        for j in range(0, np.shape(tmask_out)[1]):
            for i in range(0, np.shape(tmask_out)[2]):
                if k < top[j, i]:
                    tmask_out[k - 1, j, i] = 0
                if top[j, i] <= k <= bot[j, i]:
                    tmask_out[k - 1, j, i] = 1
                if k > bot[j, i]:
                    tmask_out[k - 1, j, i] = 0
    return tmask_out


def create_umask_from_tmask(tmask_in):
    r"""
    create_umask_from_tmask: the U-point mask from the T-point mask.

    A U-point is wet only if the T-cells on both sides of it are, as in NEMO.

    The east neighbour of the last column comes from the periodic seam.

        Parameters:
              tmask   (float): mask at t-point

        Returns:
          umask_out   (float): mask at u-point

        Formula:
            u_{i,j,k} = t_{i,j,k}\;t_{i+1,j,k}

        NEMO:
            src/OCE/DOM/dommsk.F90: umask(ji,jj,jk) = tmask(ji,jj  ,jk) * tmask(ji+1,jj  ,jk)

        See also:
            create_vmask_from_tmask, create_tmask_from_top_and_bottom
    """
    tmask_e = halo(tmask_in, halo_array=np.array([0, 0, 1]))
    return (tmask_e[..., :, :-1] * tmask_e[..., :, 1:]).astype(np.asarray(tmask_in).dtype)


def create_vmask_from_tmask(tmask_in):
    r"""
    create_vmask_from_tmask: the V-point mask from the T-point mask.

    A V-point is wet only if the T-cells on both sides of it are, as in NEMO.

    The north neighbour of the top row comes from the fold.

        Parameters:
              tmask   (float): mask at t-point

        Returns:
          vmask_out   (float): mask at v-point

        Formula:
            v_{i,j,k} = t_{i,j,k}\;t_{i,j+1,k}

        NEMO:
            src/OCE/DOM/dommsk.F90: vmask(ji,jj,jk) = tmask(ji,jj  ,jk) * tmask(ji  ,jj+1,jk)

        See also:
            create_umask_from_tmask
    """
    tmask_n = halo(tmask_in, halo_array=np.array([0, 1, 0]))
    return (tmask_n[..., :-1, :] * tmask_n[..., 1:, :]).astype(np.asarray(tmask_in).dtype)


def e3uv_from_partial_steps(e3t_0, e3t_ps, tmask, umask, vmask):
    r"""
    e3uv_from_partial_steps: U-point and V-point cell thicknesses on a grid with partial cells at the sea floor.

    With partial steps, the deepest cell of each column is thinner than its level, to follow the
    bathymetry. Older mesh files store only the T-point thicknesses; NEMO takes the thickness at a U or V
    point as the smaller of its two T-cells, which this rebuilds, so transports through the faces near
    the bottom have the right size.

    The deepest wet level of each column comes from the T mask, not from the U or V masks, which already
    require both neighbours: using them left some wet points with no thickness. The result is 0 on land.

        Parameters:
            e3t_0  (float): (nz,) 1D reference thickness per z-level
            e3t_ps (float): (y, x) partial-step bottom-cell thickness
            tmask  (bool):  (nz, y, x) T-point wet mask
            umask  (bool):  (nz, y, x) U-point wet mask
            vmask  (bool):  (nz, y, x) V-point wet mask

        Returns:
            e3u, e3v (float): (nz, y, x) reconstructed cell thicknesses,
                               zero outside umask/vmask respectively

        Note: bottom_idx (each T-column's own deepest wet level) MUST come
        from tmask, not from umask/vmask -- umask(k,i,j) already requires
        BOTH tmask(k,i,j) and tmask(k,i+1,j), so deriving it from umask
        instead double-applies that condition and can leave e3u spuriously
        zero at some genuinely-wet U-points. Hit exactly this failure
        (NaNs from a few zero-weight-but-wet columns) building the first
        version of this function -- see the GLORYS04 variance notebook's
        mesh_mask cell for the full failure mode, before it was folded in
        here.

        Callees:
            halo

        Last check date: 2026-09-05 (earlier noise-generation code)

        Formula:
            e_{3t,k} = e_{3t,k}^{0} \text{ above the deepest wet level},\quad e_{3t}^{ps} \text{ at it},\quad 0 \text{ below}
            e_{3u} = \min(e_{3t,i},\,e_{3t,i+1})\,u,\qquad e_{3v} = \min(e_{3t,j},\,e_{3t,j+1})\,v

        NEMO:
            tools/DOMAINcfg/src/domzgr.F90: e3u_0 (ji,jj,jk) = MIN( e3t_0(ji,jj,jk), e3t_0(ji+1,jj,jk) )

        See also:
            split_barotropic, grid_c.vertical_velocity
    """
    nz = e3t_0.shape[0]
    bottom_idx = tmask.sum(axis=0).astype(int) - 1     # T-point bottom level index, -1 if dry
    k_idx = np.arange(nz)[:, None, None]

    e3t_3d = np.where(k_idx < bottom_idx[None, :, :], e3t_0[:, None, None], 0.0)
    is_bottom = (k_idx == bottom_idx[None, :, :]) & (bottom_idx[None, :, :] >= 0)
    e3t_3d = np.where(is_bottom, e3t_ps[None, :, :], e3t_3d)
    e3t_3d = e3t_3d * tmask   # belt-and-braces: exactly zero on land

    # e3t is a cell-thickness (a positive scalar geometry field, not a
    # vector component), so isMeridional=False for both shifts.
    e3t_ip1 = halo(e3t_3d, halo_array=np.array([0, 0, 1]), isMeridional=False)[..., 1:]
    e3u = np.minimum(e3t_3d, e3t_ip1) * umask

    e3t_jp1 = halo(e3t_3d, halo_array=np.array([0, 1, 0]), isMeridional=False)[..., 1:, :]
    e3v = np.minimum(e3t_3d, e3t_jp1) * vmask

    return e3u, e3v


def split_barotropic(field_in, mask3d, e3=None):
    r"""
    split_barotropic: split a field into its depth mean and the departure from it.

    The barotropic part is the same at all depths, the baroclinic part varies with depth and averages to
    zero over the column. With the cell thicknesses as weights the mean is the true depth mean, as in
    NEMO's barotropic forcing; without them every level counts equally.

    Masked points are left out of the mean, and both parts are 0 on land.

        Parameters:
            field_in (float): (..., k, y, x) field, k = vertical levels
            mask3d   (float): (k, y, x) wet-point mask (1 valid, 0 for
                               land/below-bottom), broadcastable against
                               field_in
            e3       (float): Optional cell thickness, broadcastable
                               against field_in, for a volume-weighted
                               mean; if None, an unweighted mean over wet
                               levels is used instead (an approximation,
                               for when no vertical metric is available)

        Returns:
            barotropic (float): (..., y, x) vertical-mean field, masked
            baroclinic (float): (..., k, y, x) field_in - barotropic, masked

        Last check date: 2026-09-05 (earlier noise-generation code)

        Formula:
            \bar u = \frac{\sum_k e_{3,k}\,u_k\,m_k}{\sum_k e_{3,k}\,m_k},\qquad u'_k = (u_k - \bar u)\,m_k

        NEMO:
            src/OCE/DYN/dynspg_ts.F90: zu_frc(ji,jj) = SUM( e3u_0(ji,jj,:) * puu(ji,jj,:,Krhs) * umask(ji,jj,:) ) * r1_hu_0(ji,jj)

        See also:
            e3uv_from_partial_steps
    """
    k_axis = -3
    weight = mask3d if e3 is None else mask3d * e3
    num = np.sum(field_in * weight, axis=k_axis)
    den = np.sum(weight, axis=k_axis)
    with np.errstate(invalid="ignore", divide="ignore"):
        barotropic = num / den
    barotropic = np.where(den > 0, barotropic, 0.0)
    baroclinic = (field_in - barotropic[..., None, :, :]) * mask3d
    return barotropic, baroclinic


def assign_to_coarse_grid(lat_fine, lon_fine, lat_coarse, lon_coarse, mask_coarse):
    r"""
    assign_to_coarse_grid: for each fine point, the nearest wet cell of a coarse grid.

    The direction opposite to find_nearest_grid: many fine points go to each coarse cell, so a coarse
    value can be the average of all the fine points inside it rather than one sample. Distances are taken
    on the sphere, so the seam and the fold need no special case.

    Only wet coarse cells receive points. The index refers to the flattened coarse grid.

        Parameters:
            lat_fine, lon_fine     (float): fine-grid coords (deg), 2D
            lat_coarse, lon_coarse (float): coarse-grid coords (deg), 2D
            mask_coarse            (bool):  coarse wet mask, 2D

        Returns:
            idx (int): shape of lat_fine; values are flat indices into the
                       FULL flattened coarse grid (0 .. lat_coarse.size-1),
                       restricted to wet cells. Pass n_coarse_cells =
                       lat_coarse.size to bin_to_coarse, and reshape a
                       per-coarse-cell array straight back to lat_coarse's
                       2D shape (land cells just carry count 0).

        Callees:
            scipy.spatial.cKDTree

        Last check date: 2026-09-05 (earlier noise-generation code)

        Formula:
            c(p) = \arg\min_{c\ \text{wet}} \|\mathbf{x}_p - \mathbf{x}_c\|,\qquad \mathbf{x} = (\cos\varphi\cos\lambda,\ \cos\varphi\sin\lambda,\ \sin\varphi)

        See also:
            coarse_grain, bin_to_coarse
    """
    from scipy.spatial import cKDTree

    def _xyz(lat, lon):
        r"""
        _xyz: points on the unit sphere, flattened.

        A helper of assign_to_coarse_grid.

            Parameters:
                lat (float): Latitudes, in degrees
                lon (float): Longitudes, in degrees

            Returns:
                xyz (float): (n, 3) Cartesian coordinates

            Formula:
                \mathbf{x} = (\cos\varphi\cos\lambda,\ \cos\varphi\sin\lambda,\ \sin\varphi)

            See also:
                _unit_sphere
        """
        la, lo = np.deg2rad(lat), np.deg2rad(lon)
        return np.stack([(np.cos(la) * np.cos(lo)).ravel(),
                         (np.cos(la) * np.sin(lo)).ravel(),
                         np.sin(la).ravel()], axis=-1)

    valid_idx = np.where(mask_coarse.ravel())[0]
    tree = cKDTree(_xyz(lat_coarse, lon_coarse)[mask_coarse.ravel()])
    _, nearest = tree.query(_xyz(lat_fine, lon_fine))
    return valid_idx[nearest].reshape(np.shape(lat_fine))


def bin_to_coarse(field, weight, coarse_flat_idx, n_coarse_cells):
    r"""
    bin_to_coarse: weighted average of a field over the fine points of each coarse cell.

    The averaging half of coarse-graining, once assign_to_coarse_grid has grouped the points.

    Cells with no weight get 0, with their weight count 0 beside them. coarse_grain does the same with
    NaN, any number of leading axes, and the result on the coarse grid.

        Parameters:
            field, weight   (float): (y, x) or (k, y, x); `weight` is the
                                      fine wet mask (0/1) or any per-point
                                      weight
            coarse_flat_idx   (int): from assign_to_coarse_grid
            n_coarse_cells    (int): number of wet coarse cells

        Returns:
            avg    (float): (n_coarse_cells,) or (k, n_coarse_cells);
                             0 where a coarse cell got no weight
            counts (float): summed weight per coarse cell, same shape

        Reconstruct on the fine grid with
        ``avg[..., coarse_flat_idx.ravel()].reshape(field.shape)``.

        Last check date: 2026-09-05 (earlier noise-generation code)

        Formula:
            \bar f_c = \frac{\sum_{p:\,c(p)=c} w_p\, f_p}{\sum_{p:\,c(p)=c} w_p}

        See also:
            coarse_grain, assign_to_coarse_grid
    """
    flat = coarse_flat_idx.ravel()
    if field.ndim == 2:
        num = np.bincount(flat, weights=(field * weight).reshape(-1), minlength=n_coarse_cells)
        den = np.bincount(flat, weights=np.asarray(weight, dtype=float).reshape(-1), minlength=n_coarse_cells)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(den > 0, num / den, 0.0), den
    nk = field.shape[0]
    avg = np.zeros((nk, n_coarse_cells))
    cnt = np.zeros((nk, n_coarse_cells))
    for k in range(nk):
        avg[k], cnt[k] = bin_to_coarse(field[k], weight[k], coarse_flat_idx, n_coarse_cells)
    return avg, cnt


def coarse_grain(field, weight, coarse_flat_idx, shape):
    r"""
    coarse_grain: weighted average of a field onto a coarser target grid.

    A coarse model cell represents the average of what lies inside it, not one sample of it. Each fine
    point goes to the coarse cell found by assign_to_coarse_grid, and each cell gets the weighted mean of
    its points; with cell areas as weights the area integral is kept.

    NaN counts as no data, and cells with no weight are NaN. Levels and other leading axes are averaged
    separately.

        Parameters:
            field          (float): (..., ny, nx) source field; NaN counts as no data
            weight         (float): Weights, e.g. cell areas, broadcastable to field; 0 leaves a
                                    point out
            coarse_flat_idx  (int): (ny, nx) flat index of each point's target cell
            shape          (tuple): (ny, nx) of the target grid

        Returns:
            field_out (float): (..., *shape), NaN on target cells that got no weight

        Formula:
            \bar f_c = \frac{\sum_{p:\,c(p)=c} w_p\, f_p}{\sum_{p:\,c(p)=c} w_p}\quad\text{over points with data}

        See also:
            assign_to_coarse_grid, bin_to_coarse, downsample_field
    """
    field = np.asarray(field, dtype=float)
    weight = np.broadcast_to(weight, field.shape)
    valid = np.isfinite(field)
    n = int(np.prod(shape))
    flat = np.broadcast_to(np.ravel(coarse_flat_idx), (int(np.prod(field.shape[:-2])), np.size(coarse_flat_idx)))
    offset = np.arange(flat.shape[0])[:, None]*n
    index = (flat + offset).ravel()
    num = np.bincount(index, weights=(np.where(valid, field, 0.0)*weight).ravel(), minlength=flat.shape[0]*n)
    den = np.bincount(index, weights=(valid*weight).ravel(), minlength=flat.shape[0]*n)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, num/den, np.nan).reshape(field.shape[:-2] + tuple(shape))


def load_or_make_map(path, lat_coarse, lon_coarse, lat_fine, lon_fine, mask=None):
    r"""
    load_or_make_map: find_nearest_grid, saved to a file and reloaded on later calls.

    The search is slow on a fine grid and its result never changes for a pair of grids, so it is computed
    once.

    The file is checked only for its shape: delete it if the grids or the mask change.

        Parameters:
            path         (str): File of the map
            lat_coarse (float): Latitudes of the coarse grid
            lon_coarse (float): Longitudes of the coarse grid
            lat_fine   (float): Latitudes of the fine grid
            lon_fine   (float): Longitudes of the fine grid
            mask        (bool): Coarse points to map, all by default

        Returns:
            ilat (int): Fine-grid row of each coarse point, -999 off the mask
            ilon (int): Fine-grid column of each coarse point, -999 off the mask

        Formula:
            (\mathrm{ilat}, \mathrm{ilon}) = \text{find\_nearest\_grid}(\varphi^c, \lambda^c, \varphi^f, \lambda^f),\ \text{cached in a file}

        See also:
            find_nearest_grid, downsample_field
    """
    import netCDF4
    if os.path.exists(path):
        with netCDF4.Dataset(path) as nc:
            ilat, ilon = nc["ilat"][:].data, nc["ilon"][:].data
        if ilat.shape != np.shape(lat_coarse):
            raise ValueError(f"{path} has shape {ilat.shape}, expected {np.shape(lat_coarse)}")
        return ilat, ilon
    ilat, ilon = find_nearest_grid(lat_coarse, lon_coarse, lat_fine, lon_fine, mask=mask)
    with netCDF4.Dataset(path, "w") as nc:
        nc.createDimension("y", ilat.shape[0])
        nc.createDimension("x", ilat.shape[1])
        for name, values in (("ilat", ilat), ("ilon", ilon)):
            nc.createVariable(name, "i4", ("y", "x"))[...] = values
    return ilat, ilon


def coastal_taper(mask, sigma, npass):
    r"""
    coastal_taper: a weight that fades fields towards the coasts.

    Fields sampled near land are often unreliable, or would excite spurious responses in a model. This
    weight is 1 far from land and falls smoothly to 0 at the coast: the mask is smoothed and multiplied
    back by itself, npass times, each pass pushing the fade further offshore.

    It filters all axes of the mask, with no seam or fold handling.

        Parameters:
            mask  (array): Wet mask, any shape
            sigma (float): Gaussian lengthscale, in grid points
            npass   (int): Number of passes; more passes push the fade further offshore

        Returns:
            taper (float): float32 weight, 1 far from land, falling towards 0 at the coast, 0 on land

        Formula:
            \tau_0 = m,\qquad \tau_{n+1} = (G_\sigma * \tau_n)\,m,\qquad n = 0, \dots, n_{\mathrm{pass}}-1

        See also:
            smooth_2D
    """
    from scipy.ndimage import gaussian_filter
    mask = np.asarray(mask, dtype=np.float32)
    out = mask
    for _ in range(npass):
        out = gaussian_filter(out, sigma=sigma) * mask
    return out
