#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module modes

Working with snapshot matrices and modes. Snapshot matrices are (n, m) and
modes (n, r), as in dmd.py. Selectors return a boolean mask over modes, applied
with apply_mask, so that they combine with & and |.
"""
import numpy as np


def pack(fields, wets):
    r"""
    pack: stack the wet points of several fields into one snapshot matrix.

    Decompositions want one column per snapshot and no land. Several fields, such as u and v, can share
    the matrix, so their joint variability is decomposed.

        Parameters:
            fields (list): (m, ...) arrays, m snapshots each
            wets   (list): Wet mask of each field, its shape without the first axis

        Returns:
            X (float): (n, m) matrix, the wet points of each field in turn, one snapshot per column

        Formula:
            X = \begin{pmatrix} f^{(1)}[w^{(1)}] \\ f^{(2)}[w^{(2)}] \\ \vdots \end{pmatrix},\ \text{one column per snapshot}

        See also:
            unpack
    """
    return np.vstack([np.asarray(f)[:, w].T for f, w in zip(fields, wets)])


def unpack(X, wets, fill=np.nan):
    r"""
    unpack: the inverse of pack, for snapshots or modes.

        Parameters:
            X  (float): (n, m) matrix built by pack
            wets (list): The same wet masks
            fill (float): Value put on land

        Returns:
            fields (list): (m, ...) arrays, one per mask

        Formula:
            f^{(k)}[w^{(k)}] = X_{\text{rows of field } k},\qquad \text{land} = \text{fill}

        See also:
            pack
    """
    out, start = [], 0
    for w in wets:
        n = int(np.count_nonzero(w))
        f = np.full((X.shape[1],) + w.shape, fill, dtype=np.result_type(X.dtype, np.float32))
        f[:, w] = X[start:start+n].T
        out.append(f)
        start += n
    return out


def gram(X, kind="linear", degree=2, scale=None):
    r"""
    gram: the snapshot kernel matrix for kernel DMD.

    A linear kernel gives plain DMD in kernel form; a polynomial one captures some nonlinear interactions
    between snapshots.

        Parameters:
            X     (float): (n, m) snapshots, one per column
            kind    (str): 'linear', X^T X, or 'poly', (1 + X^T X / scale)^degree
            degree  (int): Degree of the polynomial kernel
            scale (float): Scale of the polynomial kernel; the mean of the diagonal of X^T X by default

        Returns:
            S2 (float): (m, m) kernel matrix. Raises ValueError for another kind

        Formula:
            K = X^{T} X \ \text{(linear)},\qquad K = \big(1 + X^{T} X / c\big)^{d} \ \text{(poly)},\qquad c = \overline{\mathrm{diag}(X^{T} X)}

        See also:
            modredkit.dmd.kernelDMD
    """
    S2 = X.T @ X
    if kind == "linear":
        return S2
    if kind == "poly":
        return (1 + S2 / (np.mean(np.diag(S2)) if scale is None else scale))**degree
    raise ValueError(f"unknown kernel {kind}")


def mask_unit_circle(dt_eig, tol=1e-3):
    r"""
    mask_unit_circle: the modes near the unit circle, neither growing nor decaying much.

        Parameters:
            dt_eig (complex): Discrete-time eigenvalues
            tol      (float): Largest | |dt_eig| - 1 | kept

        Returns:
            mask (bool): True for the modes kept

        Formula:
            \big|\,|\lambda_k| - 1\,\big| \le \text{tol}

        See also:
            mask_not_growing, apply_mask
    """
    return np.abs(np.abs(dt_eig) - 1) <= tol


def mask_not_growing(dt_eig, tol=1e-3):
    r"""
    mask_not_growing: the modes that do not grow.

        Parameters:
            dt_eig (complex): Discrete-time eigenvalues
            tol      (float): Margin above 1

        Returns:
            mask (bool): True for the modes kept

        Formula:
            |\lambda_k| \le 1 + \text{tol}

        See also:
            mask_unit_circle, apply_mask
    """
    return np.abs(dt_eig) <= 1 + tol


def mask_period(ct_eig, min_period, max_period):
    r"""
    mask_period: the modes whose period lies in a range.

        Parameters:
            ct_eig   (complex): Continuous-time eigenvalues
            min_period (float): Shortest period kept, in the time unit of 1/ct_eig
            max_period (float): Longest period kept, same unit

        Returns:
            mask (bool): True for the modes kept

        Formula:
            \frac{2\pi}{T_{max}} \le |\mathrm{Im}\,\lambda^{c}_k| \le \frac{2\pi}{T_{min}}

        See also:
            apply_mask
    """
    omega = np.abs(ct_eig.imag)
    return (omega >= 2*np.pi/max_period) & (omega <= 2*np.pi/min_period)


def mask_one_per_pair(ct_eig):
    r"""
    mask_one_per_pair: one mode of each complex-conjugate pair; real modes are all kept.

        Parameters:
            ct_eig (complex): Continuous-time eigenvalues

        Returns:
            mask (bool): True for the modes kept

        Formula:
            \mathrm{Im}\,\lambda^{c}_k \ge 0

        See also:
            modredkit.dmd.getConjugatePairsMask, apply_mask
    """
    return ct_eig.imag >= 0


def mask_top_energy(mds, amp, n):
    r"""
    mask_top_energy: the n modes that contribute most to the record.

        Parameters:
            mds (complex): (n_points, r) modes, one per column
            amp (complex): (r,) amplitudes
            n       (int): Number of modes kept

        Returns:
            mask (bool): True for the modes kept

        Formula:
            \text{the } n \text{ largest } |b_k|\,\|\phi_k\|

        See also:
            apply_mask, modredkit.dmd.DMD_amplitudes_optimal
    """
    energy = np.abs(amp) * np.linalg.norm(mds, axis=0)
    mask = np.zeros(energy.shape, dtype=bool)
    mask[np.argsort(-energy)[:n]] = True
    return mask


def apply_mask(mask, mds, *rest):
    r"""
    apply_mask: keep the selected modes, in the modes and in any per-mode arrays.

        Parameters:
            mask  (bool): (r,) modes to keep
            mds (complex): (..., r) modes, one per last-axis entry
            *rest (array): (r,) arrays, e.g. eigenvalues and amplitudes

        Returns:
            out (tuple): mds and each of rest, restricted to the modes kept

        Formula:
            \phi_{:,\text{mask}},\ a_{\text{mask}},\ \dots

        See also:
            mask_unit_circle, mask_top_energy
    """
    return (mds[..., mask], *(a[mask] for a in rest))


def time_dynamics(dt_eig, nt, kind="exact"):
    r"""
    time_dynamics: the time behaviour of each mode over a number of snapshots.

    'undamped' keeps only the phase of each eigenvalue, so modes oscillate without growing or decaying,
    for example to extend a noise forcing beyond the record.

        Parameters:
            dt_eig (complex): (r,) discrete-time eigenvalues
            nt         (int): Number of snapshots
            kind       (str): 'exact', dt_eig^t, or 'undamped', keeping only the phase of dt_eig

        Returns:
            dyn (complex): (r, nt). Raises ValueError for another kind

        Formula:
            T_{k,t} = \lambda_k^{\,t}\ \text{(exact)},\qquad T_{k,t} = e^{i t \arg \lambda_k}\ \text{(undamped)}

        See also:
            reconstruct
    """
    if kind == "undamped":
        dt_eig = np.exp(1j*np.angle(dt_eig))
    elif kind != "exact":
        raise ValueError(f"unknown kind {kind}")
    return dt_eig[:, None] ** np.arange(nt)[None, :]


def reconstruct(mds, amp, dt_eig, nt, kind="exact"):
    r"""
    reconstruct: the real snapshot matrix from modes, amplitudes and eigenvalues.

        Parameters:
            mds    (complex): (n, r) modes
            amp    (complex): (r,) amplitudes
            dt_eig (complex): (r,) discrete-time eigenvalues
            nt         (int): Number of snapshots
            kind       (str): See time_dynamics

        Returns:
            X (float): (n, nt) real part of (mds amp) dt_eig^t

        Formula:
            X = \mathrm{Re}\big[(\Phi\,\mathrm{diag}(b))\, T\big]

        See also:
            time_dynamics, modredkit.dmd.DMD_dynamics_complex
    """
    return ((mds*amp) @ time_dynamics(dt_eig, nt, kind)).real
