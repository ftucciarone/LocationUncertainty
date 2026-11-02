#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Dynamic mode decomposition: exact, kernel and windowed multi-level DMD, with amplitudes, eigenvalue
filters, reshaping helpers and plots.

Written by F. Tucciarone, 2023. Taken from the classDMD module of the earlier model-reduction code
(2026-09-06) and renamed dmd.py -- it is a module of functions, not a class.

Working / used by the kDMD notebooks:
    st_unravel, st_ravel, s_unravel, s_ravel,
    SVD, exactDMD, exactDMD_Light, kernelDMD,
    DMD_continuous_eig, DMD_amplitudes_pinv, DMD_dynamics_complex,
    getConjugatePairsMask, RemGtT, CleanUnitCircle, PlotCircle, PlotAmplitudes.

The classes ``DMD_vecmodes`` and ``DMD_phyModes`` were broken and unused, and are
gone; they are not in this repository.
"""
import copy
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.colors as colors

def RelativeInformationContent(S, r=None, verbose=False):
    r"""
    RelativeInformationContent: the rank that keeps a given share of the singular values.

    Truncating a decomposition keeps the dominant structures and drops noise; the question is where to cut.
    This picks the rank whose singular values add up to a given share of their total, or takes a rank
    given directly.

    The share uses the singular values S, not S^2, so it is not the share of variance; rank_from_energy
    in pod uses S^2.

        Parameters:
            S     (float): Singular values, in decreasing order
            r  (int/float): An int is a rank, capped at len(S); a float is the share to keep, as a
                            fraction, or as a percentage if above 1; None keeps all
            verbose (bool): Plot the cumulative share and the rank chosen

        Returns:
            r (int): Rank

        Formula:
            r = \min\Big\{ r : \frac{\sum_{k \le r} S_k}{\sum_k S_k} \gtrsim \tau \Big\}

        See also:
            SVD, modredkit.pod.rank_from_energy
    """
    # Compute dimension
    ns = np.size(S)
    # Dimensionality reduction parameter
    if r is not None:
        if isinstance(r, int):
            r = min(r, ns)
        elif isinstance(r, float):
            ric = np.cumsum(S, dtype=np.float64)/np.sum(S, dtype=np.float64)
            if r > 1.0:
                r = int(np.abs(ric - r/100).argmin()) + 1
            else:
                r = int(np.abs(ric - r).argmin()) + 1
        else:
            r = ns

    else:
        r = ns

    # If verbose then plot
    if verbose:
        func = np.cumsum(S, dtype=np.float64)/np.sum(S, dtype=np.float64)
        t = r-1
        plt.title("Relative Information Content for SVD")
        plt.xlim([0,ns])
        plt.ylim([0,1.1])
        plt.plot(func)
        # Vertical line
        plt.plot([t, t], [0, func[t]], color='gray', linewidth=1.5, linestyle="--") 
        # Horizontal line
        plt.plot([0, t], [func[t], func[t]], color='gray', linewidth=1.5, linestyle="--")
        # Point on the curve
        plt.scatter([t, ], [func[t], ], 10, color='gray')
        # Point on x-axis
        plt.scatter([t, ], [0, ], 10, color='gray')
        # Point on y-axis
        plt.scatter([0, ], [func[t], ], 10, color='gray')
        plt.show()
    # Return
    return r

def SVD( X, r, verbose=False):
    r"""
    SVD: thin singular value decomposition of a snapshot matrix, truncated to a rank.

    The SVD gives the orthogonal structures that best represent the snapshots, ordered by importance; it is
    the first step of DMD and of POD.

        Parameters:
            X     (float): (n, m) matrix of snapshots, one per column
            r  (int/float): Rank, or share of information to keep (see RelativeInformationContent);
                            None keeps all
            verbose (bool): Print the shapes

        Returns:
            U  (float): (n, r) left singular vectors
            S  (float): (r,) singular values
            Vh (float): (r, m) right singular vectors, conjugate-transposed

        Formula:
            X \approx U_r\, \mathrm{diag}(S_r)\, V_r^{*}

        See also:
            RelativeInformationContent, exactDMD
    """
    # Compute dimension
    ns = np.size(X,1) # snapshots number
    # Compute SVD
    U, S, Vh = np.linalg.svd(np.array(X, dtype=np.float64), full_matrices=False)
    # Print out information
    if verbose:
        print(" ")
        print(" Arrays dimensions within SVD")
        print("                        np.shape(U)",np.shape(U))
        print("     np.shape(U.conj().transpose())",np.shape(U.conj().transpose()))
        print("                        np.shape(X)",np.shape(X))
        print("                       np.shape(Vh)",np.shape(Vh))
        print(" ")
    # Dimensionality reduction parameter: None keeps all, an int is the rank, a
    # float is an information threshold, as a fraction or a percentage (100.0 keeps
    # all). py-kits fix: r used to be ignored.
    if r is None:
        r = S.size
    elif isinstance(r, (float, np.floating)):
        r = RelativeInformationContent( S, float(r), verbose=verbose)
    r = min(int(r), S.size)
    U = U[:,:r]
    S = S[:r]
    Vh = Vh[:r,:]
    # Return
    return U, S, Vh

def exactDMD( Xb, Xa, r, SVD_out=False, verbose=False):
    r"""
    exactDMD: exact dynamic mode decomposition of pairs of snapshots.

    DMD finds spatial modes that each evolve with one complex frequency and growth rate, from snapshots
    one time step apart: the best linear map from each snapshot to the next, reduced to the rank of the
    truncated SVD.

        Parameters:
            Xb     (float): (n, m) snapshots before, one per column
            Xa     (float): (n, m) the same snapshots one step later
            r   (int/float): Rank of the SVD, or share of information to keep
            SVD_out  (bool): Also return the SVD and the eigenvectors of A~
            verbose  (bool): Print the shapes

        Returns:
            dmd_mds (complex): (n, r) DMD modes
            dt_eig  (complex): (r,) discrete-time eigenvalues
            ns          (int): Number of modes
            U, S, Vh  (float): The SVD of Xb, with SVD_out only
            W       (complex): (r, r) eigenvectors of A~, with SVD_out only

        Formula:
            X_b = U S V^{*},\qquad \tilde A = U^{*} X_a V S^{-1},\qquad \tilde A W = W \Lambda,\qquad \Phi = X_a V S^{-1} W

        See also:
            exactDMD_Light, DMD_continuous_eig, DMD_amplitudes_optimal
    """
    # Compute SVD with reduced order 'r'
    U, S, Vh = SVD(np.array(Xb, dtype=np.float64), r)
    # Compute V
    V = Vh.conj().transpose()
    if verbose:
        print(" ")
        print(" Arrays dimensions within SVD")
        print("                        np.shape(U)",np.shape(U))
        print("     np.shape(U.conj().transpose())",np.shape(U.conj().transpose()))
        print("                       np.shape(Xa)",np.shape(Xa))
        print("                        np.shape(V)",np.shape(V))
        print(" ")
    # Low rank dynamics
    #   \tilde{A} = U^{\dagger} D' V S^{-1}
    Atilde = np.matmul(np.matmul(np.matmul(U.conj().transpose(), np.array(Xa, dtype=np.float64)), V), np.diag(1./S))
    # Discrete-time eigenvalues 
    #   \tilde{A} W = W \Lambda
    dt_eig, W = np.linalg.eig(Atilde) 
    # DMD modes
    #   \Phi = D' V S^{-1} W
    dmd_mds = np.matmul(np.matmul(np.matmul(np.array(Xa, dtype=np.float64), V), np.diag(1./S)), W)
    # number of modes
    ns = np.shape(dmd_mds)[1]
    # Return
    if SVD_out:
        return dmd_mds, dt_eig, ns, U, S, Vh, W
    else:
        return dmd_mds, dt_eig, ns

def exactDMD_Light( data, r, verbose=False):
    r"""
    exactDMD_Light: exactDMD of a single snapshot matrix, keeping as few copies in memory as possible.

    For large fields: the pairs are the snapshot matrix without its last and first columns, read as views,
    and the intermediate arrays are reused and freed as soon as possible.

        Parameters:
            data   (float): (n, nt) snapshots, one per column
            r   (int/float): Rank of the SVD, or share of information to keep
            verbose  (bool): Print sizes and shapes

        Returns:
            Phi    (complex): (n, r) DMD modes
            dt_eig (complex): (r,) discrete-time eigenvalues
            ns         (int): Number of modes

        Formula:
            X_b = X_{:,:-1},\quad X_a = X_{:,1:},\quad \tilde A = U^{*} X_a V S^{-1},\quad \Phi = X_a V S^{-1} W

        See also:
            exactDMD
    """
    # Print shapes
    if verbose:
        print("Performing SVD with full data matrix (classDMD.DMD_LightRAM)")
        print(" ")
        print("  Snapshot data input ")
        print("                size: "+str(np.size(data)))
        print("               shape: "+str(np.shape(data)))
        print("           in memory: "+str(round(data.nbytes/1024/1024,2))+" MB")

    # Compute SVD (with reduced order 'r' to be implemented)
    U, S, work = SVD(np.array(data[:,:-1], dtype=np.float64), r)
    # Compute Complex-Conjugate of V and scale it by the eigenvalues
    #   thus define work = V S^{-1}
    work = np.matmul(work.conj().transpose(), np.diag(1./S))

    # Print SVD shapes
    if verbose:
        print(" ")
        print(" Arrays dimensions within SVD (classDMD.DMD_LightRAM)")
        print("                        np.shape(U)",np.shape(U))
        print("     np.shape(U.conj().transpose())",np.shape(U.conj().transpose()))
        print("                        np.shape(V)",np.shape(work))

    # Low rank dynamics
    #   \tilde{A} = U^{\dagger} D' ( V S^{-1} )
    #             = U^{\dagger} * ( D' * work )
    #
    # 1. Define \Phi <- D^{\prime} V S^{-1} = data[:,1:] * work, that has the same dimensions as the 
    #    modes that result from the procedure. The low rank dynamics is thus \tilde{A} = U^{\dagger} * Phi
    Phi = np.matmul( np.array(data[:,1:], dtype=np.float64), work)
    # 2. Notice that \tilde{A} has the same dimension as V, and that work has finished its scope, hence
    #    one can define work = \tilde{A}
    work = np.matmul(U.conj().transpose(), Phi)
    # 3. Notice that U is now out of scope, one can safely delete it (along with S)
    del U
    # Discrete-time eigenvalues 
    #   \tilde{A} W = W \Lambda
    # That with our substitutions has become
    #   work * W = W \Lambda
    dt_eig, W = np.linalg.eig(work) 

    if verbose:
        print(" ")
        print(" Constructing DMD modes (classDMD.DMD_LightRAM)")
        print(" ")

    # DMD modes
    #   \Phi = D' V S^{-1} W
    #        = \Phi * W
    Phi = np.matmul(Phi, W)
    # number of modes
    ns = np.shape(Phi)[1]
    # Return
    return Phi, dt_eig, ns

def kernelDMD( S2, X, dt, r,Y=None, verbose=False, noModes=False):
    r"""
    kernelDMD: kernel dynamic mode decomposition from a Gram matrix (Williams, Rowley and Kevrekidis 2015).

    Kernel DMD works with the inner products of the snapshots rather than the snapshots themselves, and a
    nonlinear kernel lets it approximate the Koopman operator of nonlinear dynamics with modest cost.

    The Gram matrix of all snapshots gives G~ (pairs before) and A~ (pairs after); the eigendecomposition
    of G~, truncated to r, yields the Koopman approximation, whose eigenvalues and modes are returned.

        Parameters:
            S2     (float): (nt, nt) Gram matrix k(x_i, x_j) of all snapshots
            X      (float): (nt-1, n) snapshots before, one per row, to project the modes on
            dt     (float): Time between snapshots
            r        (int): Rank kept from the eigendecomposition of G~
            Y      (float): Unused
            verbose (bool): Unused
            noModes (bool): Return only the eigenvalues

        Returns:
            kDMD_mds (complex): (n, r) Koopman modes, None with noModes
            ct_eig   (complex): (r,) continuous-time eigenvalues, log(mu)/dt
            mu       (complex): (r,) discrete-time eigenvalues
            A        (complex): (nt-1, r) Q Sigma^+ V, None with noModes
            Vr       (complex): (r, r) right eigenvectors of K, None with noModes
            nx           (int): Rows of Phi = Q Sigma V, None with noModes
            nm           (int): Number of modes, None with noModes

        Formula:
            \tilde G = Q \Sigma^2 Q^{T},\qquad \hat K = (\Sigma^{+} Q^{T})\, \tilde A\, (Q \Sigma^{+}),\qquad \hat K V = V M,\qquad \lambda^{c} = \frac{\log \mu}{\Delta t}

        See also:
            exactDMD, modredkit.modes.gram
    """
    import scipy.linalg as sLA
    #
    # Algorithm 3, step 1: Compute the Gramiam and \tilde{A}
    # \tilde{G} = k(x_i, x_j)
    Gtilde = S2[:-1,:-1]
    # \tilde{A} = k(y_i, x_j)
    Atilde = S2[1:,:-1]
    #
    #  Algorithm 3, step 2: Compute the eigendecomposition of the Gramiam to
    #  obtain Q and \Sigma
    # (Right) Eigenvalue of \tilde{G}
    #   \tilde{G} Q = Q \Lambda
    # print("Gtilde",np.shape(Gtilde))
    G_eival, Q = sLA.eigh(Gtilde)
    idx = np.argsort(G_eival)[::-1][:r]
    G_eival = G_eival[idx]
    Q = Q[:, idx]
    # print("Q",np.shape(Q))

    G_sigma = np.diag(np.sqrt(np.abs(G_eival))) # Prevent non-positive issues
    S_psinv = np.diag(1.0 / np.sqrt(G_eival))

    #
    # Proposition 1: Definition of the approximation of the Koopman operator \hat{K}
    #   \hat{K} = ( \Sigma^{+} * Q^{T} ) \hat{A} ( Q * \Sigma^{+} ) = C * A * C^T (11)
    #   with ancilliary matrix C = \Sigma^{+} * Q^T
    C = S_psinv @ Q.T
    Khat = C @ Atilde @ C.T
    #
    # Left and Right Eigenvalue of \hat{K}
    #   \hat{K} \hat{V} = \hat{V} \Mu (Right eigenvalue, Vr[:,i])
    #   \hat{W} \hat{K} = \Mu \hat{W} (Left  eigenvalue, Wl[:,i])
    mu_Keival, Wl, Vr = sLA.eig(Khat, left=True, right=True)
    # Continuous-time eigenvalues
    ct_eig = np.log(mu_Keival, dtype=complex)/dt
    
    if noModes == False:
        weights = np.diag(1.0 / np.diag((Wl.conj().T @ Vr))) # not just real part or we loose info !
        #weights = np.diag(1.0/np.diag((Wl.conj().T @ Vr).real)) 
        #
        # Equation 20, compute Xi
        #   \Xi = ( \hat{W}^{T} \Sigma^{+} Q^{T} X )^{T}
        Xi = weights @ Wl.conj().T @ S_psinv @ Q.T @ X
        Xi = Xi.T
        #
        # Equation 17, compute Phi
        #   \Phi = Q \Sigma \hat{V}
        Phi = Q @ G_sigma @ Vr 
        # print("Phi",np.shape(Phi))
        # print("Xi",np.shape(Xi))
        #
        # Compute DMD modes
        kDMD_mds = Xi#np.matmul(Xi,Phi)
        # print("kDMD_mds",np.shape(kDMD_mds))
        # print("dt",dt)
        
        A = Q @ S_psinv @ Vr
        # Return "modes", continuous eigenvalue, discrete eigenvalue, spatial dimension, number of modes
        return kDMD_mds, ct_eig, mu_Keival, A, Vr, np.shape(Phi)[0], np.shape(Phi)[1]
    else:
        return None, ct_eig, mu_Keival, None, None, None, None
    
def DMD_amplitudes_pinv( dmd_mds, X0, apply=False ):

    r"""
    DMD_amplitudes_pinv: amplitudes of the modes that best fit one snapshot, by pseudoinverse.

    The simplest choice, usually the first snapshot. DMD_amplitudes_optimal fits the whole record instead.

        Parameters:
            dmd_mds (complex): (n, r) DMD modes
            X0        (float): (n,) snapshot to fit, usually the first
            apply      (bool): Also return the modes scaled by their amplitudes

        Returns:
            amp    (complex): (r,) amplitudes, pinv(dmd_mds) X0
            scaled (complex): (n, r) dmd_mds * amp, with apply only

        Formula:
            b = \Phi^{+} x_0

        See also:
            DMD_amplitudes_optimal
    """
    # DMD amplitudes (through Moore-Penrose inverse)
    amp = np.matmul(np.linalg.pinv(dmd_mds), X0) 
    # Return
    if not apply:
        return amp
    else:
        return amp, dmd_mds*amp

def DMD_continuous_eig( dt_eig, dt ):
    r"""
    DMD_continuous_eig: continuous-time eigenvalues from discrete-time ones.

    Discrete eigenvalues depend on the time step; continuous ones give the growth rate (real part) and
    angular frequency (imaginary part) in physical units.

        Parameters:
            dt_eig (complex): Discrete-time eigenvalues
            dt       (float): Time between snapshots

        Returns:
            ct_eig (complex): Continuous-time eigenvalues; real part the growth rate, imaginary
                              part the angular frequency

        Formula:
            \lambda^{c} = \frac{\log \lambda}{\Delta t}

        See also:
            DMD_dynamics_complex
    """
    # Continuous-time eigenvalues
    ct_eig = np.log(dt_eig, dtype=complex)/dt
    # Return 
    return ct_eig

def DMD_dynamics_complex( dmd_mds, dt_eig, dt, nt ):
    
    r"""
    DMD_dynamics_complex: rebuild snapshots from all the modes, each with its complex exponential.

        Parameters:
            dmd_mds (complex): (n, r) modes, already scaled by their amplitudes
            dt_eig  (complex): (r,) discrete-time eigenvalues
            dt        (float): Time between snapshots
            nt          (int): Number of snapshots to rebuild

        Returns:
            recX   (float): (n, nt) real part of the reconstruction
            ct_eig (complex): (r,) continuous-time eigenvalues

        Formula:
            X(t) = \mathrm{Re}\sum_k \phi_k\, e^{\lambda^{c}_k t},\qquad t = 0, \Delta t, \dots, (n_t - 1)\Delta t

        See also:
            DMD_dynamics_noPairs, modredkit.modes.reconstruct
    """
    # Compute dimension
    nm = np.size(dmd_mds,1) # modes number
    # Continuous-time eigenvalues
    ct_eig = DMD_continuous_eig(dt_eig,dt)
    # Time vector
    t = np.linspace(0, nt-1, nt, endpoint=True) * dt
    # Time dynamics
    time_dynamics = np.zeros((nm, nt), dtype=complex)
    for iter in range(nt):
        time_dynamics[:,iter] = np.exp( ct_eig * t[iter] )
    # Reconstruction
    recX = np.matmul(dmd_mds, time_dynamics) 
    # Return 
    return recX.real, ct_eig

def getConjugatePairsMask( ct_eig ):
    
    r"""
    getConjugatePairsMask: keep one eigenvalue of each complex-conjugate pair.

    For a real field, modes come in conjugate pairs that carry the same information. Keeping the one with
    Im >= 0, and every real mode, halves the storage, which is how NEMO's noise input stores them.

        Parameters:
            ct_eig (complex): Eigenvalues

        Returns:
            mask   (bool): True for the eigenvalues kept
            nm_out  (int): Number kept

        Formula:
            \text{keep } k \iff \mathrm{Im}\,\lambda_k \ge 0

        See also:
            DMD_dynamics_noPairs, modredkit.modes.mask_one_per_pair
    """
    mask = np.asarray(ct_eig).imag >= 0
    # Compute the length
    nm_out = np.sum(mask)
    # Return 
    return mask, nm_out # = removeConjugatePairs( dmd_mds, dt_eig, ct_eig, amp )

def DMD_dynamics_noPairs( dmd_mds, dt_eig, dt, nt ):
    
    r"""
    DMD_dynamics_noPairs: rebuild snapshots from one mode per complex-conjugate pair.

    A complex mode stands for itself and its conjugate, so it counts twice; a real mode counts once.

        Parameters:
            dmd_mds (complex): (n, r) modes kept, already scaled by their amplitudes
            dt_eig  (complex): (r,) their discrete-time eigenvalues
            dt        (float): Time between snapshots
            nt          (int): Number of snapshots to rebuild

        Returns:
            recX   (float): (n, nt) reconstruction
            ct_eig (complex): (r,) continuous-time eigenvalues

        Formula:
            X(t) = \mathrm{Re}\sum_k w_k\, \phi_k\, e^{\lambda^{c}_k t},\qquad w_k = 2 \text{ if } \mathrm{Im}\,\lambda_k \ne 0,\ 1 \text{ otherwise}

        See also:
            getConjugatePairsMask, DMD_dynamics_complex
    """
    # Continuous-time eigenvalues
    ct_eig = DMD_continuous_eig(dt_eig,dt)
    # Time vector
    t = np.linspace(0, nt-1, nt, endpoint=True) * dt
    # Reconstruction
    weight = np.where(np.abs(np.asarray(dt_eig).imag) > 0, 2.0, 1.0)
    recX = np.matmul(dmd_mds * weight, np.exp(np.outer(ct_eig, t))).real
    # Return 
    return recX, ct_eig

def DMD_amplitudes_2014Jovanovich( S, V, W, dt_eig):

    r"""
    DMD_amplitudes_2014Jovanovich: optimal amplitudes of Jovanovic et al. (2014), computed from the SVD.

    The amplitudes do not change under the rotation of the SVD, so the global least-squares problem of
    DMD_amplitudes_optimal is solved in SVD coordinates, which is cheaper. An earlier version built a
    square Vandermonde matrix and was right only without truncation.

        Parameters:
            S       (float): (r,) singular values of the snapshots before
            V       (float): (m, r) right singular vectors, Vh.conj().T; (r, m) is also accepted
            W     (complex): (r, r) eigenvectors of the reduced operator A~
            dt_eig (complex): (r,) discrete-time eigenvalues

        Returns:
            amp (complex): (r,) amplitudes

        Formula:
            b = \arg\min_b \big\| S V^{*} - W\, \mathrm{diag}(b)\, \mathrm{Vand}(\lambda) \big\|_F

        See also:
            DMD_amplitudes_optimal
    """
    S = np.asarray(S).ravel()
    V = np.asarray(V)
    W = np.asarray(W)
    dt_eig = np.asarray(dt_eig)
    r = dt_eig.size
    if V.ndim == 2 and V.shape[0] == r and V.shape[1] != r:   # tolerate V given as (r, m)
        V = V.conj().T
    G = S[:, None] * V.conj().T                               # Sigma V^H = X1 in POD coords, (r, m)
    # amplitudes are invariant under the orthogonal POD rotation U, so solve the
    # identical global-LSQ problem in POD coordinates: modes -> W (= Y), data -> G
    return DMD_amplitudes_optimal(W, dt_eig, G)


def DMD_amplitudes_optimal(dmd_mds, dt_eig, X, return_energy=False):
    r"""
    DMD_amplitudes_optimal: amplitudes that fit the whole record, not only its first snapshot.

    Fitting one snapshot leaves the amplitudes at the mercy of that snapshot's noise. This minimises the
    misfit over all snapshots at once (Jovanovic et al. 2014), from the modes and eigenvalues alone, so it
    suits both exactDMD and kernelDMD. Fit the same field the modes describe.

        Parameters:
            dmd_mds    (complex): (n, r) modes, one per column
            dt_eig     (complex): (r,) discrete-time eigenvalues
            X            (float): (n, m) snapshots to fit, X[:, j] ~ Phi diag(b) dt_eig^j
            return_energy (bool): Also return each mode's contribution

        Returns:
            b    (complex): (r,) amplitudes
            energy (float): (r,) |b_k| ||Phi_k|| ||dt_eig_k^j||, to rank modes; with return_energy only

        Formula:
            b = \arg\min_b \| X - \Phi\, \mathrm{diag}(b)\, \mathrm{Vand} \|_F,\qquad P b = q,\quad P = (\Phi^{*}\Phi) \circ \overline{\mathrm{Vand}\,\mathrm{Vand}^{*}},\quad q = \mathrm{diag}(\Phi^{*} X\, \mathrm{Vand}^{*})

        See also:
            DMD_amplitudes_pinv, modredkit.modes.mask_top_energy
    """
    Phi = np.asarray(dmd_mds)
    lam = np.asarray(dt_eig, dtype=complex)
    X = np.asarray(X)
    m = X.shape[1]
    j = np.arange(m)
    Vand = lam[:, None] ** j[None, :]                              # (r, m)
    P = (Phi.conj().T @ Phi) * np.conj(Vand @ Vand.conj().T)       # (r, r) Hermitian PD
    q = np.diag((Phi.conj().T @ X) @ Vand.conj().T)               # (r,)
    try:
        import scipy.linalg as _sla
        b = _sla.cho_solve(_sla.cho_factor(P), q)
    except Exception:
        b = np.linalg.lstsq(P, q, rcond=None)[0]
    if not return_energy:
        return b
    tenv = np.sqrt(np.sum(np.abs(Vand) ** 2, axis=1))             # ||lam_k**j||_2
    energy = np.abs(b) * np.linalg.norm(Phi, axis=0) * tenv
    return b, energy


def DefineXandXp( S, verbose=False ):
    
    r"""
    DefineXandXp: split a snapshot matrix into the pairs before and after one step.

        Parameters:
            S     (float): (n, nt) snapshots, one per column
            verbose (bool): Print sizes and shapes

        Returns:
            X  (float): (n, nt-1) S[:, :-1]
            Xp (float): (n, nt-1) S[:, 1:]

        Formula:
            X = S_{:,:-1},\qquad X' = S_{:,1:}

        See also:
            exactDMD
    """
    X = S[:,:-1]
    Xp = S[:,1:] 
    # Print shapes
    if verbose:
        print(r"Defining X and X^{\prime} (dmd.DefineXandXp)")
        print(" ")
        print("  Snapshot data input ")
        print("                size: "+str(np.size(S)))
        print("               shape: "+str(np.shape(S)))
        print(" ")
        print("  X data output ")
        print("                size: "+str(np.size(X)))
        print("               shape: "+str(np.shape(X)))
        print(" ")
        print(r"  X^{\prime} data output ")
        print("                size: "+str(np.size(Xp)))
        print("               shape: "+str(np.shape(Xp)))
        print(" ")
    # Return 
    return  X,  Xp
    

def ShiftStack( S, m, s ):
    
    r"""
    ShiftStack: stack time-shifted copies of the snapshots, a delay embedding.

    Stacking delayed copies lets DMD represent dynamics that one snapshot alone cannot, such as standing
    oscillations or a few dominant frequencies in a small state.

        Parameters:
            S (float): (n, m) snapshots, one per column
            m   (int): Number of snapshots
            s   (int): Number of shifted copies added

        Returns:
            ssS (float): ((s+1) n, m-s) S[:, :m-s] stacked over S[:, i+1:m-s+i+1], i = 0 .. s-1

        Formula:
            \begin{pmatrix} S_{:,0:m-s} \\ S_{:,1:m-s+1} \\ \vdots \\ S_{:,s:m} \end{pmatrix}

        See also:
            exactDMD
    """

    ssS = S[...,:m-s]
    for i in range(s):
        ssS = np.concatenate((ssS, S[...,i+1:m-s+i+1]), axis=0) 
    # Return 
    return ssS


def RemGtT( dmd_mds, dt_eig, ct_eig, amp, treshold = 1.0 ):
    
    r"""
    RemGtT: drop the modes whose discrete-time eigenvalue exceeds a threshold in modulus.

    With the default threshold, the growing modes, which would blow up a reconstruction or forecast.

        Parameters:
            dmd_mds (complex): (n, r) modes
            dt_eig  (complex): (r,) discrete-time eigenvalues
            ct_eig  (complex): (r,) continuous-time eigenvalues
            amp     (complex): (r,) amplitudes
            treshold  (float): Largest |dt_eig| kept, plus 1e-3

        Returns:
            dmd_mds, dt_eig, ct_eig, amp (complex): The same, for the modes kept

        Formula:
            \text{keep } k \iff |\lambda_k| \le \tau + 10^{-3}

        See also:
            CleanUnitCircle, modredkit.modes.mask_not_growing
    """
    
    # Identify modes with 
    lam_mod = np.abs(dt_eig)
    mask = lam_mod <= treshold + 1.0e-3
    
    dmd_mds = dmd_mds[...,mask]
    ct_eig = ct_eig[mask]
    dt_eig = dt_eig[mask]
    amp = amp[mask] 
    # Return 
    return dmd_mds, dt_eig, ct_eig, amp


def CleanUnitCircle( dmd_mds, ct_eig, dt_eig, amp, tol ):    
    r"""
    CleanUnitCircle: keep only the modes near the unit circle, neither growing nor decaying much.

        Parameters:
            dmd_mds (complex): (n, r) modes
            ct_eig  (complex): (r,) continuous-time eigenvalues
            dt_eig  (complex): (r,) discrete-time eigenvalues
            amp     (complex): (r,) amplitudes
            tol       (float): Largest | |dt_eig| - 1 | kept

        Returns:
            dmd_mds, ct_eig, dt_eig, amp (complex): The same, for the modes kept
            nm (int): Number of modes kept

        Formula:
            \text{keep } k \iff \big|\,|\lambda_k| - 1\,\big| \le \text{tol}

        See also:
            RemGtT, modredkit.modes.mask_unit_circle
    """
    
    non_osc_idx = np.where(np.abs(np.abs(dt_eig) - 1) > tol )

    circ_amp = np.delete(amp, non_osc_idx)
    circ_mds = np.delete(dmd_mds, non_osc_idx, axis=1)
    circ_ct_eig = np.delete(ct_eig, non_osc_idx)
    circ_dt_eig = np.delete(dt_eig, non_osc_idx)

    return circ_mds, circ_ct_eig, circ_dt_eig, circ_amp, np.shape(circ_mds)[-1]


def SortByFrequency( dmd_mds, ct_eig, dt_eig, amp ):
    r"""
    SortByFrequency: order modes from low to high frequency.

        Parameters:
            dmd_mds (complex): (n, r) modes
            ct_eig  (complex): (r,) continuous-time eigenvalues
            dt_eig  (complex): (r,) discrete-time eigenvalues
            amp     (complex): (r,) amplitudes

        Returns:
            dmd_mds, ct_eig, dt_eig, amp (complex): The same, sorted

        Formula:
            \text{order by } |\mathrm{Im}\,\lambda^{c}_k|

        See also:
            PlotAmplitudes
    """
    
    sort_index = np.argsort(np.abs(np.asarray(ct_eig).imag))

    # Sorting the arrays by omega
    sorted_amp = amp[sort_index]
    sorted_mds = dmd_mds[:, sort_index]
    sorted_ct_eig = ct_eig[sort_index]
    sorted_dt_eig = dt_eig[sort_index]

    return sorted_mds, sorted_ct_eig, sorted_dt_eig, sorted_amp

def st_unravel(phy_in, shape=False, verbose=False):
    r"""
    st_unravel: reshape a space-time array into a matrix with space in rows and time in columns.

    Decompositions work on matrices of snapshots, one per column.

        Parameters:
            phy_in (float): (nt, nz, ny, nx) or (nt, ...) array
            shape   (bool): Also return the sizes needed to undo it
            verbose (bool): Print shapes

        Returns:
            TaS_out (float): (nz*ny*nx, nt) matrix
            nt        (int): Number of snapshots, with shape only
            ns      (tuple): Spatial shape, with shape only

        Formula:
            X_{p,t} = f_{t, k, j, i},\qquad p = (k, j, i) \text{ flattened}

        See also:
            st_ravel, s_unravel
    """
    import math
    import sys
    nt = np.shape(phy_in)[0]
    ns = np.shape(phy_in)[1:]
    TaS_out = phy_in.reshape(nt, math.prod(ns)).T 
    # Print shapes
    if verbose:
        print("Unraveling data (classDMD.unravel)")
        print(" ")
        print("  Physical data input ")
        print("                size: "+str(np.size(phy_in)))
        print("               shape: "+str(np.shape(phy_in)))
        print(" ")
        print("  Computational data output ")
        print("                size: "+str(np.size(TaS_out)))
        print("               shape: "+str(np.shape(TaS_out)))
        print(" ")
        print("  Data output dimensions")
        print("                 nt: "+str(nt))
        print("                 ns: "+str(ns))
        print("          in memory: "+str(round(TaS_out.nbytes/1024/1024,2))+" MB")
        print(" ")
    # Return 
    if shape:
        return TaS_out, nt, ns
    else:
        return TaS_out

def s_unravel(phy_in, shape=False, verbose=False):
    r"""
    s_unravel: flatten a spatial array into one column.

        Parameters:
            phy_in (float): (nz, ny, nx) or other spatial array
            shape   (bool): Also return its shape
            verbose (bool): Print shapes

        Returns:
            TaS_out (float): (nz*ny*nx,) array
            ns      (tuple): Spatial shape, with shape only

        Formula:
            x_p = f_{k, j, i},\qquad p = (k, j, i) \text{ flattened}

        See also:
            s_ravel, st_unravel
    """
    import math
    import sys
    ns = np.shape(phy_in)
    TaS_out = phy_in.reshape(math.prod(ns)).T 
    # Print shapes
    if verbose:
        print("Unraveling data (classDMD.unravel)")
        print(" ")
        print("  Physical data input ")
        print("                size: "+str(np.size(phy_in)))
        print("               shape: "+str(np.shape(phy_in)))
        print(" ")
        print("  Computational data output ")
        print("                size: "+str(np.size(TaS_out)))
        print("               shape: "+str(np.shape(TaS_out)))
        print(" ")
        print("  Data output dimensions")
        print("                 ns: "+str(ns))
        print("          in memory: "+str(round(TaS_out.nbytes/1024/1024,2))+" MB")
        print(" ")
    # Return
    if shape:
        return TaS_out, ns
    else:
        return TaS_out

def st_ravel(TaS_in, nt, ns, verbose=False):
    """
    st_ravel: reshape a matrix with space in rows and time in columns back into a space-time array.

        Parameters:
            TaS_in (float): (nz*ny*nx, nt) matrix
            nt       (int): Number of columns
            ns     (tuple): Spatial shape
            verbose (bool): Print shapes

        Returns:
            phy_out (float): (nt, *ns) array

        Formula:
            f_{t, k, j, i} = X_{p,t}

        See also:
            st_unravel
    """
    phy_out = (TaS_in.T).reshape(tuple(np.append(nt, ns))) 
    # Print shapes
    if verbose:
        print("Raveling data (classDMD.ravel)")
        print(" ")
        print("  Computational data input ")
        print("                size: "+str(np.size(TaS_in)))
        print("               shape: "+str(np.shape(TaS_in)))
        print(" ")
        print("  Physical data output ")
        print("                size: "+str(np.size(phy_out)))
        print("               shape: "+str(np.shape(phy_out)))
        print(" ")
    # Return 
    return phy_out

def s_ravel(TaS_in, ns, verbose=False):
    """
    s_ravel: reshape one column back into a spatial array.

        Parameters:
            TaS_in (float): (nz*ny*nx,) array
            ns     (tuple): Spatial shape
            verbose (bool): Print shapes

        Returns:
            phy_out (float): Array of shape ns

        Formula:
            f_{k, j, i} = x_p

        See also:
            s_unravel
    """
    phy_out = (TaS_in.T).reshape(tuple(ns)) 
    # Print shapes
    if verbose:
        print("Raveling data (classDMD.ravel)")
        print(" ")
        print("  Computational data input ")
        print("                size: "+str(np.size(TaS_in)))
        print("               shape: "+str(np.shape(TaS_in)))
        print(" ")
        print("  Physical data output ")
        print("                size: "+str(np.size(phy_out)))
        print("               shape: "+str(np.shape(phy_out)))
        print(" ")
    # Return 
    return phy_out

def PlotAmplitudes(ct_eig, amp, title=None, split=None, freq_time=None, saveout=None, filename=None):
    r"""
    PlotAmplitudes: plot mode amplitudes against frequency.

    Shows which frequencies carry the variability, with an optional period marked to separate slow and
    fast modes.

        Parameters:
            ct_eig  (complex): Continuous-time eigenvalues, in s-1
            amp     (complex): Amplitudes; their modulus is plotted
            title       (str): Added to the title
            split     (float): Period, in days, marked by two vertical lines
            freq_time (float): Periods, in days, used as x ticks in place of frequencies
            saveout    (bool): Save to filename instead of showing
            filename    (str): File to save to

        Returns:
            None

        Formula:
            \big(\mathrm{Im}\,\lambda^{c}_k,\ |b_k|\big)

        See also:
            PlotCircle, SortByFrequency
    """
    plt.figure(figsize=(5,5))
    plt.rcParams['figure.dpi'] = 50 

    if split is not None:
        tresh = 1.0e5*np.pi/(split*86400.)
        plt.axvline(x=tresh, color="orange" )
        plt.axvline(x=-tresh, color="orange" )
    plt.plot(1.0e5*ct_eig.imag, np.abs(amp), '+')
     
    if title is not None:
        plt.title(r'Amplidute distribution '+ title)
    else:
        plt.title(r'Amplidute distribution')
    plt.grid(which='both', axis='both')
    plt.grid(which='both', axis='both')
 
    if freq_time is None:
        plt.xlabel(r'Im ($\omega$) $\times 10^{-5}$')
    else:   
        freqs =  np.append(freq_time, -freq_time[::-1])
        ticks = 1.0e5*np.pi / (freqs*86400.)
        labels = [ str(abs(x)) for x in freqs ]
        plt.xticks(ticks, labels, rotation=90)
        plt.xlabel(r'Im ($\omega$) (days)')
    plt.ylabel(r'Amplitude')
    plt.ylim([0, 60])
    if saveout is not None: 
        plt.savefig(filename, dpi=200, bbox_inches='tight', pad_inches=0)
    else:
        plt.show()


def PlotCircle(ct_eig, dt_eig, title=None, split=None, freq_time=None, saveout=None, filename=None):
    r"""
    PlotCircle: plot the discrete-time eigenvalues against the unit circle.

    Eigenvalues inside the circle decay, outside it grow, and on it oscillate steadily.

        Parameters:
            ct_eig  (complex): Continuous-time eigenvalues
            dt_eig  (complex): Discrete-time eigenvalues
            title       (str): Added to the title
            split     (float): Period, in days, separating slow and fast modes
            freq_time (float): Periods, in days, for scale separation
            saveout    (bool): Save to filename instead of showing
            filename    (str): File to save to

        Returns:
            None

        Formula:
            \big(\mathrm{Re}\,\lambda_k,\ \mathrm{Im}\,\lambda_k\big) \text{ with } |z| = 1

        See also:
            PlotAmplitudes, CleanUnitCircle
    """
    plt.figure(figsize=(5,5))
    plt.rcParams['figure.dpi'] = 50 

    theta = np.linspace(0,2*np.pi,100)
    plt.plot(np.sin(theta),np.cos(theta))
    if split is not None:
        id_split = np.argmin(np.abs(np.abs(ct_eig.imag) - np.pi/(split*86400.)))
        plt.plot([0, dt_eig[id_split].real], [0, dt_eig[id_split].imag], color="orange" )
        plt.plot([0, dt_eig[id_split].real], [0, -dt_eig[id_split].imag], color="orange")

    if freq_time is None: 
        plt.plot(dt_eig.real, dt_eig.imag, ".")
    else:
        for i in range(len(freq_time)):
            idm_temp = np.where(np.abs(ct_eig.imag) <= np.pi/(freq_time[i] * 86400.))
            mask_temp = np.full(len(ct_eig.imag), False, dtype=bool)
            mask_temp[idm_temp] = True
            plt.plot(dt_eig.real[mask_temp], dt_eig.imag[mask_temp], '.',label=str(freq_time[i]) + " days")
            del mask_temp, idm_temp
            plt.legend(loc='best')
    plt.grid(which='both', axis='both')

    if title is not None:
        plt.title(r'Eigenvalue $\lambda$ disposition '+ title)
    else:
        plt.title(r'Eigenvalue $\lambda$ disposition')

    plt.xlim([-1.05, 1.05])
    plt.ylim([-1.05, 1.05])
    plt.xlabel(r'Re ($\lambda$)')
    plt.ylabel(r'Im ($\lambda$)')
    if saveout is not None: 
        plt.savefig(filename, dpi=200, bbox_inches='tight', pad_inches=0)
    else:
        plt.show()
        




def SaveNetCDF(dmd_corr, dmd_uncorr, utrend, vtrend, filepath ):
    # Save outputs
    r"""
    SaveNetCDF: write corrected and uncorrected physical modes, their frequencies and the trends to one file.

    It expects objects of mode classes that were removed from this module, so it cannot run as it stands.

        Parameters:
            dmd_corr   (object): Modes written as mode_c (umode_real_c, ...)
            dmd_uncorr (object): Modes written as mode_r (umode_real_r, ...)
            utrend      (float): (z, y, x) u trend
            vtrend      (float): (z, y, x) v trend
            filepath      (str): File to create

        Returns:
            None

        Formula:
            \{u, v, w\}\text{mode\_\{real,imag\}\_\{c,r\}},\ \omega_{c},\ \omega_{r},\ \text{trends}

        See also:
            orcakit.netcdf.write_modes
    """
    fout = Dataset(filepath, 'w', format='NETCDF4')
    
    fout.createDimension('mode_c', dmd_corr.nm)
    fout.createDimension('mode_r', dmd_uncorr.nm)
    
    fout.createDimension('y', dmd_corr.ny)
    fout.createDimension('x', dmd_corr.nx)
    fout.createDimension('z', dmd_corr.nz)
    
    yid = fout.createVariable('y', 'f8', ('y',))
    xid = fout.createVariable('x', 'f8', ('x',))
    zid = fout.createVariable('z', 'f8', ('z',))
    
    mcid = fout.createVariable('mode_c', 'i4', ('mode_c',))
    mrid = fout.createVariable('mode_r', 'i4', ('mode_r',))
    
    umid = fout.createVariable('uco', 'f8', ('z','y','x',))
    vmid = fout.createVariable('vco', 'f8', ('z','y','x',))
    wmid = fout.createVariable('wco', 'f8', ('z','y','x',))
    
    wcid = fout.createVariable('omega_c', 'f8', ('mode_c',))
    wrid = fout.createVariable('omega_r', 'f8', ('mode_r',))
    
    ucreid = fout.createVariable('umode_real_c', 'f8', ('mode_c','z','y','x',))
    ucimid = fout.createVariable('umode_imag_c', 'f8', ('mode_c','z','y','x',))
    vcreid = fout.createVariable('vmode_real_c', 'f8', ('mode_c','z','y','x',))
    vcimid = fout.createVariable('vmode_imag_c', 'f8', ('mode_c','z','y','x',))
    wcreid = fout.createVariable('wmode_real_c', 'f8', ('mode_c','z','y','x',))
    wcimid = fout.createVariable('wmode_imag_c', 'f8', ('mode_c','z','y','x',))
    
    urreid = fout.createVariable('umode_real_r', 'f8', ('mode_r','z','y','x',))
    urimid = fout.createVariable('umode_imag_r', 'f8', ('mode_r','z','y','x',))
    vrreid = fout.createVariable('vmode_real_r', 'f8', ('mode_r','z','y','x',))
    vrimid = fout.createVariable('vmode_imag_r', 'f8', ('mode_r','z','y','x',))
    wrreid = fout.createVariable('wmode_real_r', 'f8', ('mode_r','z','y','x',))
    wrimid = fout.createVariable('wmode_imag_r', 'f8', ('mode_r','z','y','x',))
    
    # Add attributes
    mcid.long_name = 'Mode index (for correction drift)'
    mrid.long_name = 'Mode index (for random noise)'
    
    # xid.long_name = 'Ocean X axis'
    # yid.long_name = 'Ocean Y axis'
    # zid.long_name = 'Ocean Z axis'
    # xid.units = 'km'
    # yid.units = 'km'
    # zid.units = 'km'
    
    wcid.long_name = 'Continuous eigenvalues for correction drift'
    wrid.long_name = 'Continuous eigenvalues for random noise'
    wcid.units = 's^-1'
    wrid.units = 's^-1'
    
    umid.long_name = 'Mean of zonal correction drift'
    vmid.long_name = 'Mean of meridional correction drift'
    wmid.long_name = 'Mean of vertical correction drift'
    umid.units = 'm/s'
    vmid.units = 'm/s'
    wmid.units = 'm/s'
    
    ucreid.long_name = 'Zonal DMD modes for correction drift (real part)'
    vcreid.long_name = 'Meridional DMD modes for correction drift (real part)'
    wcreid.long_name = 'Vertical DMD modes for correction drift (real part)'
    ucreid.units = 'm/s'
    vcreid.units = 'm/s'
    wcreid.units = 'm/s'
    
    ucimid.long_name = 'Zonal DMD modes for correction drift (imag part)'
    vcimid.long_name = 'Meridional DMD modes for correction drift (imag part)'
    wcimid.long_name = 'Vertical DMD modes for correction drift (imag part)'
    ucimid.units = 'm/s'
    vcimid.units = 'm/s'
    wcimid.units = 'm/s'
    
    urreid.long_name = 'Zonal DMD modes for random noise (real part)'
    vrreid.long_name = 'Meridional DMD modes for random noise (real part)'
    wrreid.long_name = 'Vertical DMD modes for random noise (real part)'
    urreid.units = 'm/s'
    vrreid.units = 'm/s'
    wrreid.units = 'm/s'
    
    urimid.long_name = 'Zonal DMD modes for random noise (imag part)'
    vrimid.long_name = 'Meridional DMD modes for random noise (imag part)'
    wrimid.long_name = 'Vertical DMD modes for random noise (imag part)'
    vrimid.units = 'm/s'
    urimid.units = 'm/s'
    wrimid.units = 'm/s'
    
    # Write data
    mcid[:] = range(dmd_corr.nm)
    mrid[:] = range(dmd_uncorr.nm)
    
    # yid[:] = dmd_uncorr.y
    # xid[:] = dmd_uncorr.x
    # zid[:] = dmd_uncorr.z
    
    wcid[:] = dmd_corr.ct_eig.imag
    wrid[:] = dmd_uncorr.ct_eig.imag
    
    umid[...] = utrend 
    vmid[...] = vtrend
    #wmid[...] = wc
    
    ucreid[...] = dmd_corr.Umod.real
    vcreid[...] = dmd_corr.Vmod.real
    wcreid[...] = dmd_corr.Wmod.real
    
    ucimid[...] = dmd_corr.Umod.imag 
    vcimid[...] = dmd_corr.Vmod.imag
    wcimid[...] = dmd_corr.Wmod.imag
    
    urreid[...] = dmd_uncorr.Umod.real
    vrreid[...] = dmd_uncorr.Vmod.real
    wrreid[...] = dmd_uncorr.Wmod.real
    
    urimid[...] = dmd_uncorr.Umod.imag
    vrimid[...] = dmd_uncorr.Vmod.imag
    wrimid[...] = dmd_uncorr.Wmod.imag
    
    # Close file
    fout.close()
    
    print("")
    print('Program terminates')


def mrDMD( Xraw, r, dt, L ):
    r"""
    mrDMD: windowed multi-level DMD.

    Level l splits the record into 2^l windows and runs exactDMD on each, so slow modes come from the long
    windows and fast ones from the short. Growing modes are dropped and amplitudes fit each window's first
    snapshot.

    Unlike the full mrDMD of Kutz et al. (2016), slow modes are not subtracted between levels and the data
    are not subsampled. A remainder after the last whole window is dropped.

        Parameters:
            Xraw (float): (n, m) snapshots, one per column
            r      (int): Rank of each window's SVD
            dt   (float): Time between snapshots
            L      (int): Number of levels

        Returns:
            out_mds (complex): (n, K) modes of all windows, one per column
            out_dt  (complex): (K,) discrete-time eigenvalues
            out_ct  (complex): (K,) continuous-time eigenvalues
            out_amp (complex): (K,) amplitudes
            out_idx      (int): (K, 3) level, window and mode of each column

        Formula:
            \text{for } l = 0..L-1,\ j = 0..2^l-1:\quad \text{exactDMD}\big(X_{:,\ j s:(j+1) s}\big),\qquad s = \lfloor m / 2^l \rfloor

        See also:
            exactDMD, RemGtT
    """
    ns = np.size(Xraw,1)
    out_mds, out_dt, out_ct, out_amp, out_idx = [], [], [], [], []
    for l in range(L):
        for j in range(2**l):
            s = ns // 2**l
            X = Xraw[:,j*s:(j+1)*s]
            dmd_mds, dt_eig, nm = exactDMD( X[:,:-1], X[:,1:], r )
            ct_eig = DMD_continuous_eig( dt_eig, dt )
            amp = np.linalg.lstsq( dmd_mds, X[:,0], rcond=None )[0]
            # Drop growing modes only
            dmd_mds, dt_eig, ct_eig, amp = RemGtT( dmd_mds, dt_eig, ct_eig, amp )
            nm = np.size(dt_eig)
            out_mds.append(dmd_mds); out_dt.append(dt_eig); out_ct.append(ct_eig); out_amp.append(amp)
            out_idx.append(np.column_stack((np.full(nm,l), np.full(nm,j), np.arange(nm))))
    return ( np.concatenate(out_mds,axis=1), np.concatenate(out_dt), np.concatenate(out_ct),
             np.concatenate(out_amp), np.concatenate(out_idx).astype(int) )
