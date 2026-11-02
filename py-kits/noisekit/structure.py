#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module structure

Spatial structure of the noise by an elliptic (SPDE, Matern) operator: assemble
A = kappa2*M - L on a masked grid, solve it with white noise, normalise to unit
variance. Prototype status: from the earlier noise-generation notebooks, not yet used by a run.
assemble_A does not wrap periodically, so it suits a regional box.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

NU = 1.0


def tensor_field(w_maj, w_min, theta, shape, nu=NU):
    r"""
    tensor_field: diffusion tensor and screening coefficient of an anisotropic correlation.

    Random fields with a given correlation length and orientation can be drawn by solving a diffusion-type
    equation. The tensor stretches correlations along the major axis by (w_maj/w_min)^2 relative to the
    minor one, rotated by theta, and kappa^2 sets the minor length.

        Parameters:
            w_maj (float): Major length scale, in metres; scalar or (ny, nx)
            w_min (float): Minor length scale, in metres; scalar or (ny, nx)
            theta (float): Angle of the major axis from x, in radians; scalar or (ny, nx)
            shape (tuple): (ny, nx) of the result
            nu    (float): Smoothness parameter of the correlation

        Returns:
            Dxx    (float): (ny, nx), dimensionless
            Dyy    (float): (ny, nx), dimensionless
            Dxy    (float): (ny, nx), dimensionless
            kappa2 (float): (ny, nx) screening coefficient 8 nu / w_min^2, in m-2

        Formula:
            D = R(\theta)\,\mathrm{diag}\big((w_{maj}/w_{min})^2,\ 1\big)\,R(\theta)^{T},\qquad \kappa^2 = \frac{8\nu}{w_{min}^2}

        See also:
            assemble_A, build_filter
    """
    Dmaj, Dmin = (w_maj / w_min) ** 2, 1.0
    c, s = np.cos(theta), np.sin(theta)
    Dxx = Dmaj * c ** 2 + Dmin * s ** 2
    Dyy = Dmaj * s ** 2 + Dmin * c ** 2
    Dxy = (Dmaj - Dmin) * s * c
    kappa2 = 8.0 * nu / w_min ** 2
    return (np.broadcast_to(Dxx, shape).copy(), np.broadcast_to(Dyy, shape).copy(),
            np.broadcast_to(Dxy, shape).copy(), np.broadcast_to(kappa2, shape).copy())


def assemble_A(Dxx, Dyy, Dxy, kappa2, e1, e2, wet):
    r"""
    assemble_A: the sparse system of the correlation operator on a masked grid.

    A = kappa^2 M - L, where L is the anisotropic diffusion operator with the tensor of tensor_field and M
    the cell areas. Solving A x = b for white noise b gives a correlated field. It is discretised with the
    metrics e1, e2 on any C-grid point type.

    Faces towards land or off the box carry no flux; there is no periodic wrap. A is symmetrised, so it
    is symmetric positive definite.

        Parameters:
            Dxx    (float): (ny, nx) tensor component
            Dyy    (float): (ny, nx) tensor component
            Dxy    (float): (ny, nx) cross component, centred on the cells
            kappa2 (float): (ny, nx) screening coefficient, in m-2
            e1     (float): (ny, nx) cell widths along x, in metres
            e2     (float): (ny, nx) cell widths along y, in metres
            wet    (array): (ny, nx) 1 on water, 0 on land

        Returns:
            A (csr_matrix): (ny*nx, ny*nx) system matrix
            M (dia_matrix): (ny*nx, ny*nx) lumped mass matrix, e1*e2 on the diagonal

        Formula:
            A = \kappa^2 M - \nabla\cdot(D\nabla),\qquad M = \mathrm{diag}(e_1 e_2),\qquad \text{flux through a face} = D\,\frac{e_{across}}{e_{along}}\,\delta

        See also:
            tensor_field, build_filter, pcg_jacobi
    """
    ny, nx = Dxx.shape
    idx = np.arange(ny * nx).reshape(ny, nx)
    idxp = np.pad(idx, 1, constant_values=-1)             # -1 = off the edge of the box: no periodic wrap
    nb = dict(E=idxp[1:-1, 2:], W=idxp[1:-1, :-2], N=idxp[2:, 1:-1], S=idxp[:-2, 1:-1],
              NE=idxp[2:, 2:], NW=idxp[2:, :-2], SE=idxp[:-2, 2:], SW=idxp[:-2, :-2])
    wetf = wet.astype(float)

    def gather(field, nbidx):
        r"""
        gather: values of a field at one neighbour of every cell.

            Parameters:
                field  (float): (ny, nx) field
                nbidx    (int): (ny, nx) flat index of the neighbour, -1 off the box

            Returns:
                valid   (bool): (ny, nx) True where the neighbour exists
                values (float): (ny, nx) neighbour values; the cell's own value off the box

            Formula:
                g_{j,i} = f_{\text{nb}(j,i)},\ \text{or } f_{j,i} \text{ off the box}

            See also:
                assemble_A
        """
        valid = nbidx >= 0
        return valid, np.where(valid, field.ravel()[np.where(valid, nbidx, 0)].reshape(ny, nx), field)

    def face(dcoef, nbidx, e_along, e_across):
        r"""
        face: coefficient of the flux through one face of every cell.

        The face carries no flux if the neighbour is off the box or either cell is land.

            Parameters:
                dcoef  (float): (ny, nx) diffusion coefficient at the face
                nbidx    (int): (ny, nx) flat index of the neighbour across the face
                e_along (float): (ny, nx) cell width across the face, averaged with the neighbour's
                e_across (float): (ny, nx) length of the face, averaged with the neighbour's

            Returns:
                cols   (int): (ny*nx,) column of each coefficient; the cell itself where there is no flux
                coeff (float): (ny*nx,) coefficient dcoef e_across / e_along, 0 where there is no flux

            Formula:
                c = D\,\frac{\overline{e}_{across}}{\overline{e}_{along}}\ \text{if both cells are wet, else } 0

            See also:
                assemble_A
        """
        valid, wet_nb = gather(wetf, nbidx)
        ok = valid & (wet_nb > 0) & (wet > 0)
        _, ea = gather(e_along, nbidx); _, ec = gather(e_across, nbidx)
        coeff = np.where(ok, dcoef * (0.5 * (e_across + ec)) / (0.5 * (e_along + ea)), 0.0)
        return np.where(ok, nbidx, idx).ravel(), coeff.ravel()

    _, Dxx_e = gather(Dxx, nb["E"]); Dxx_e = 0.5 * (Dxx + Dxx_e)
    _, Dxx_w = gather(Dxx, nb["W"]); Dxx_w = 0.5 * (Dxx + Dxx_w)
    _, Dyy_n = gather(Dyy, nb["N"]); Dyy_n = 0.5 * (Dyy + Dyy_n)
    _, Dyy_s = gather(Dyy, nb["S"]); Dyy_s = 0.5 * (Dyy + Dyy_s)

    rows, cols, data = [], [], []
    diag = kappa2 * e1 * e2
    for dcoef, key, ea, ec in ((Dxx_e, "E", e1, e2), (Dxx_w, "W", e1, e2),
                               (Dyy_n, "N", e2, e1), (Dyy_s, "S", e2, e1)):
        col, coeff = face(dcoef, nb[key], ea, ec)
        rows.append(idx.ravel()); cols.append(col); data.append(-coeff)
        diag += coeff.reshape(ny, nx)
    for key, sgn in (("NE", 1), ("SW", 1), ("NW", -1), ("SE", -1)):        # mixed (cross) term, cell-centred Dxy
        valid, wet_nb = gather(wetf, nb[key])
        ok = valid & (wet_nb > 0) & (wet > 0)
        coeff = np.where(ok, sgn * Dxy / 2.0, 0.0)
        rows.append(idx.ravel()); cols.append(np.where(ok, nb[key], idx).ravel()); data.append(-coeff.ravel())
    rows.append(idx.ravel()); cols.append(idx.ravel()); data.append(diag.ravel())

    A = sp.coo_matrix((np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
                      shape=(ny * nx, ny * nx)).tocsr()
    A = 0.5 * (A + A.T)
    M = sp.diags((e1 * e2).ravel())
    return A, M


def pcg_jacobi(A, b, Pinv, tol=1e-10, maxit=5000):
    r"""
    pcg_jacobi: conjugate gradients with a Jacobi (diagonal) preconditioner.

    An iterative solver for symmetric positive definite systems, such as the correlation operator, that
    needs only products with the matrix.

        Parameters:
            A      (matrix): Symmetric positive definite system
            b       (float): Right-hand side
            Pinv    (float): Inverse of the diagonal preconditioner, same shape as b
            tol     (float): Stop when |r| / |b| falls below this
            maxit     (int): Largest number of iterations

        Returns:
            x     (float): Solution
            hist   (list): Relative residual |r| / |b| after each iteration

        Formula:
            x_{k+1} = x_k + \alpha_k p_k,\quad r_{k+1} = r_k - \alpha_k A p_k,\quad z = P^{-1} r,\quad p_{k+1} = z_{k+1} + \beta_k p_k

        See also:
            cheby_jacobi, assemble_A
    """
    x = np.zeros_like(b)
    r = b - A @ x
    z = Pinv * r
    p = z.copy()
    rz_old = r @ z
    bnorm = np.linalg.norm(b)
    hist = []
    for k in range(maxit):
        Ap = A @ p
        alpha = rz_old / (p @ Ap)
        x = x + alpha * p
        r = r - alpha * Ap
        rn = np.linalg.norm(r) / bnorm
        hist.append(rn)
        if rn < tol:
            break
        z = Pinv * r
        rz_new = r @ z
        beta = rz_new / rz_old
        p = z + beta * p
        rz_old = rz_new
    return x, hist


def cheby_jacobi(A, b, Pinv, lam_min, lam_max, tol=1e-10, maxit=5000):
    r"""
    cheby_jacobi: Chebyshev iteration with a Jacobi preconditioner.

    Unlike conjugate gradients, it needs no inner products, only one matrix-vector product per step, so it
    parallelises without global reductions. It needs bounds of the preconditioned spectrum instead.

        Parameters:
            A       (matrix): Symmetric positive definite system
            b        (float): Right-hand side
            Pinv     (float): Inverse of the diagonal preconditioner, same shape as b
            lam_min  (float): Lower bound of the eigenvalues of Pinv A
            lam_max  (float): Upper bound of the eigenvalues of Pinv A
            tol      (float): Stop when |r| / |b| falls below this
            maxit      (int): Largest number of iterations

        Returns:
            x     (float): Solution
            hist   (list): Relative residual |r| / |b| after each iteration

        Formula:
            x_k = x_{k-1} + d_k,\quad d_{k+1} = \rho_k\rho_{k+1}\, d_k + \frac{2\rho_{k+1}}{\delta}\, P^{-1} r_k,\quad \theta = \tfrac{\lambda_{max}+\lambda_{min}}{2},\ \delta = \tfrac{\lambda_{max}-\lambda_{min}}{2}

        See also:
            pcg_jacobi
    """
    theta = 0.5 * (lam_max + lam_min)
    delta = 0.5 * (lam_max - lam_min)
    sigma = theta / delta
    x = np.zeros_like(b)
    r = b - A @ x                     # r_0
    d = (Pinv * r) / theta            # d_1
    rho = 1.0 / sigma                 # rho_1 = sigma_0/sigma_1
    bnorm = np.linalg.norm(b)
    hist = []
    for k in range(maxit):
        x = x + d                     # x_k = x_{k-1} + d_k
        r = r - A @ d                 # r_k = r_{k-1} - A d_k  (the ONE matvec per iteration)
        rn = np.linalg.norm(r) / bnorm
        hist.append(rn)
        if rn < tol:
            break
        rho_new = 1.0 / (2 * sigma - rho)     # rho_{k+1} = sigma_k/sigma_{k+1}
        gamma = 2 * rho_new / delta           # gamma_k
        eta = rho * rho_new                   # eta_k = sigma_{k-1}/sigma_{k+1}
        d = eta * d + gamma * (Pinv * r)      # d_{k+1}  (boxed recursion, D3)
        rho = rho_new
    return x, hist


def build_filter(w_maj, w_min, theta, e1, e2, wet):
    r"""
    build_filter: factorise the correlation operator once, for many draws.

    Drawing many correlated fields with the same structure only needs the factorisation once; each draw
    is then a cheap solve.

        Parameters:
            w_maj (float): Major length scale, in metres; scalar or (ny, nx)
            w_min (float): Minor length scale, in metres; scalar or (ny, nx)
            theta (float): Angle of the major axis from x, in radians; scalar or (ny, nx)
            e1    (float): (ny, nx) cell widths along x, in metres
            e2    (float): (ny, nx) cell widths along y, in metres
            wet   (array): (ny, nx) 1 on water, 0 on land

        Returns:
            lu (SuperLU): LU factorisation of A, see assemble_A

        Formula:
            LU = A(w_{maj}, w_{min}, \theta)

        See also:
            draw, assemble_A
    """
    Dxx, Dyy, Dxy, kappa2 = tensor_field(w_maj, w_min, theta, wet.shape)
    A, _ = assemble_A(Dxx, Dyy, Dxy, kappa2, e1, e2, wet)
    return spla.splu(A.tocsc())


def draw(lu, e1, e2, wet, rng):
    r"""
    draw: one spatially correlated random field.

    White noise is scaled by sqrt(e1 e2), so the result does not depend on the cell size, and passed
    through the factorised operator.

    The variance of the result varies in space; unit_variance_scale gives the factor to make it 1.

        Parameters:
            lu (SuperLU): From build_filter
            e1   (float): (ny, nx) cell widths along x, in metres
            e2   (float): (ny, nx) cell widths along y, in metres
            wet  (array): (ny, nx) 1 on water, 0 on land
            rng (Generator): numpy random generator

        Returns:
            field (float): (ny, nx) correlated field, 0 on land; its variance varies in space, see
                           unit_variance_scale

        Formula:
            x = A^{-1}\big(\sqrt{e_1 e_2}\ \xi\big),\qquad \xi \sim \mathcal{N}(0, 1)

        See also:
            build_filter, unit_variance_scale
    """
    b = np.sqrt(e1 * e2) * rng.standard_normal(wet.shape) * wet
    return lu.solve(b.ravel()).reshape(wet.shape) * wet


def unit_variance_scale(lu, e1, e2, wet, rng, samples=128):
    r"""
    unit_variance_scale: pointwise standard deviation of the draws, by Monte Carlo.

    Near coasts and where the grid changes, draws have more or less variance than elsewhere. Dividing
    draws by this scale makes their variance 1 everywhere, so their amplitude can be set by another field.

    Points with almost no variance, such as land, get 1.

        Parameters:
            lu      (SuperLU): From build_filter
            e1        (float): (ny, nx) cell widths along x, in metres
            e2        (float): (ny, nx) cell widths along y, in metres
            wet       (array): (ny, nx) 1 on water, 0 on land
            rng   (Generator): numpy random generator
            samples     (int): Number of draws

        Returns:
            scale (float): (ny, nx) standard deviation; 1 where it is below 1e-8, e.g. on land

        Formula:
            s = \sqrt{\frac{1}{N}\sum_{n=1}^{N} x_n^2}

        See also:
            draw
    """
    s2 = np.zeros(wet.shape)
    for _ in range(samples):
        s2 += draw(lu, e1, e2, wet, rng) ** 2
    s = np.sqrt(s2 / samples)
    return np.where(s > 1e-8, s, 1.0)
