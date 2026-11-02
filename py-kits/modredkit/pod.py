#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module pod

Snapshot POD, written for py-kits from the theory in
the POD notes of the earlier model-reduction work. An older eigenproblem-based
version existed in the earlier DMD code; neither is in this repository.

Snapshot matrices are (n, m): n spatial points, m snapshots, as in dmd.py.
Build them with dmd.st_unravel.
"""
import numpy as np


def rank_from_energy(S, energy=0.99):
    r"""
    rank_from_energy: the smallest rank whose modes hold a share of the total variance.

        Parameters:
            S      (float): Singular values, in decreasing order
            energy (float): Share of the variance, sum of S^2, to keep

        Returns:
            r (int): Rank, at most len(S)

        Formula:
            r = \min\Big\{ r : \frac{\sum_{k \le r} S_k^2}{\sum_k S_k^2} \ge \tau \Big\}

        See also:
            pod, modredkit.dmd.RelativeInformationContent
    """
    frac = np.cumsum(S**2) / np.sum(S**2)
    return int(min(np.searchsorted(frac, energy) + 1, S.size))


def pod(X, r=None, weights=None, remove_mean=True):
    r"""
    pod: proper orthogonal decomposition of a snapshot matrix.

    POD finds the spatial patterns that capture the most variance, each with a time series: the basis of
    reduced models. With weights, such as cell volumes, the modes are orthonormal in that weighted inner
    product, so cells count by their size.

    Weights must be strictly positive: remove land points first.

        Parameters:
            X         (float): (n, m) snapshots, one per column
            r           (int): Modes kept, None for all
            weights   (float): (n,) strictly positive weights, or None; remove land points first
            remove_mean (bool): Subtract the time mean first

        Returns:
            modes  (float): (n, r) spatial modes
            S      (float): (r,) singular values, of the weighted matrix
            coeffs (float): (r, m) time coefficients, S times the right singular vectors
            mean   (float): (n,) time mean, zeros without remove_mean

        Formula:
            W^{1/2}(X - \bar X) = U S V^{T},\qquad \text{modes} = W^{-1/2} U_r,\qquad \text{coeffs} = S_r V_r^{T}

        See also:
            reconstruct, pod_multi, rank_from_energy
    """
    X = np.asarray(X, dtype=np.float64)
    mean = X.mean(axis=1) if remove_mean else np.zeros(X.shape[0])
    A = X - mean[:, None]
    if weights is None:
        U, S, Vh = np.linalg.svd(A, full_matrices=False)
        modes = U
    else:
        w = np.sqrt(np.asarray(weights, dtype=np.float64))
        if np.any(w <= 0):
            raise ValueError("weights must be strictly positive; remove land points first")
        U, S, Vh = np.linalg.svd(w[:, None] * A, full_matrices=False)
        modes = U / w[:, None]
    r = S.size if r is None else min(int(r), S.size)
    return modes[:, :r], S[:r], S[:r, None] * Vh[:r], mean


def reconstruct(modes, coeffs, mean=None):
    r"""
    reconstruct: rebuild the snapshot matrix from the output of pod.

        Parameters:
            modes  (float): (n, r) spatial modes
            coeffs (float): (r, m) time coefficients
            mean   (float): (n,) time mean to add back, or None

        Returns:
            X (float): (n, m) snapshots

        Formula:
            X = \Phi\, A + \bar X

        See also:
            pod
    """
    X = modes @ coeffs
    return X if mean is None else X + mean[:, None]


def pod_multi(fields, r=None, weights=None, remove_mean=True, tol=1e-12):
    r"""
    pod_multi: joint POD of several fields that share time, such as u, v and w.

    The modes of the fields are found together from one time correlation matrix summed over the fields, so
    each mode is a coherent pattern across all of them. This is the method of snapshots, cheap when there
    are fewer snapshots than points.

    Modes with an eigenvalue below tol times the largest are dropped.

        Parameters:
            fields   (list): (n_i, m) snapshot matrix of each field
            r         (int): Modes kept, None for all
            weights  (list): Per field: None, (n_i,) diagonal weights or an (n_i, n_i) matrix
            remove_mean (bool): Subtract each field's time mean
            tol     (float): Modes with eigenvalue below tol times the largest are dropped

        Returns:
            modes  (list): (n_i, r) spatial modes of each field, jointly orthonormal in the
                           weighted product
            S     (float): (r,) singular values
            coeffs (float): (r, m) time coefficients
            means  (list): (n_i,) time mean of each field

        Formula:
            C = \sum_f A_f^{T} W_f A_f = V S^2 V^{T},\qquad \text{modes}_f = A_f V S^{-1},\qquad \text{coeffs} = S V^{T}

        See also:
            pod
    """
    weights = [None]*len(fields) if weights is None else weights
    fields = [np.asarray(X, dtype=np.float64) for X in fields]
    means = [X.mean(axis=1) if remove_mean else np.zeros(X.shape[0]) for X in fields]
    anom = [X - mu[:, None] for X, mu in zip(fields, means)]

    def apply_w(W, A):
        r"""
        apply_w: apply one field's weights to its anomalies.

            Parameters:
                W (float): None, (n,) diagonal weights or an (n, n) matrix
                A (float): (n, m) anomalies

            Returns:
                WA (float): (n, m) weighted anomalies; A itself when W is None

            Formula:
                W A,\ \text{with } W \text{ diagonal or full}

            See also:
                pod_multi
        """
        if W is None:
            return A
        W = np.asarray(W, dtype=np.float64)
        return W[:, None]*A if W.ndim == 1 else W @ A

    C = sum(A.T @ apply_w(W, A) for A, W in zip(anom, weights))
    eigval, eigvec = np.linalg.eigh(C)
    ids = np.argsort(eigval)[::-1]
    eigval, eigvec = eigval[ids], eigvec[:, ids]
    keep = eigval > tol*eigval[0]
    if r is not None:
        keep[int(r):] = False
    S = np.sqrt(eigval[keep])
    V = eigvec[:, keep]
    modes = [A @ V / S for A in anom]
    return modes, S, S[:, None]*V.T, means
