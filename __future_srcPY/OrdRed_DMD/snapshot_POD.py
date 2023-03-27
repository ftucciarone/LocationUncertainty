"""" Module to perform the Snapshot Proper Orthogonal Decomposition
     POD) procedure """

# Load modules
# ------------
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt

# Subroutines
# -----------


def spot_POD3F(U, V, W, wU, wV, wW, nm, check_POD):
    "Snapshot POD of 2D velocity data"

    [nt, nu] = np.shape(U)
    # nv = np.shape(V[0, :])

    # Detrending
    Um = np.mean(U, axis=0)
    Vm = np.mean(V, axis=0)
    Wm = np.mean(W, axis=0)

    Uf = U - Um
    Vf = V - Vm
    Wf = W - Wm

    # Build temporal correlation matrix
    C = inner_prod(Uf, Uf, wU)
    C += inner_prod(Vf, Vf, wV)
    C += inner_prod(Wf, Wf, wW)
    C /= nt

    # Solve eigen problem
    eigval, eigvec = linalg.eigh(C)
    eigval = eigval.real
    eigvec = eigvec.real
    ids = np.argsort(eigval)            # sort eigvals
    ids = ids[::-1]                     # reverse order
    PC = eigval[ids]
    PC = PC[:nm]
    Tmod = eigvec[:, ids]   # sort temporal modes
    Tmod = Tmod[:,:nm].transpose()
    del C, eigval, eigvec

    # Scaling temporal modes
    for i in range(nm):
        Tmod[i, :] *= np.sqrt(nt*np.maximum(PC[i], 0.))

    # Build spatial modes
    Umod = np.matmul(Tmod, Uf)
    Vmod = np.matmul(Tmod, Vf)
    Wmod = np.matmul(Tmod, Wf)
    for i in range(nm):
        if PC[i] >= 1e-15:
            Umod[i, :] /= (nt*PC[i])
            Vmod[i, :] /= (nt*PC[i])
            Wmod[i, :] /= (nt*PC[i])
        else:
            Umod[i, :] *= 0.
            Vmod[i, :] *= 0.
            Wmod[i, :] *= 0.

    if check_POD:

        print('Check if the spatial modes are orthonormal')
        Cx = inner_prod(Umod, Umod, wU)
        Cx += inner_prod(Vmod, Vmod, wV)
        Cx += inner_prod(Wmod, Wmod, wW)
        plt.figure()
        plt.imshow(Cx)
        plt.colorbar()
        plt.xlabel('modes')
        plt.ylabel('modes')
        plt.title('Orthogonality of spatial modes')
        plt.show()
        print(np.sum(np.diag(Cx))/nt)
        print('Check if the temporal modes are orthogonal')
        Ct = np.matmul(Tmod, Tmod.transpose())
        for m in range(nm):
            Ct[m, m] /= (nt*PC[m])
        plt.figure()
        plt.imshow(Ct)
        plt.colorbar()
        plt.xlabel('modes')
        plt.ylabel('modes')
        plt.title('Orthogonality of temporal modes (rescaled by PCs)')
        plt.show()

        print('Check POD reconstruction')
        Uerr = np.matmul(Tmod.transpose(), Umod)
        Verr = np.matmul(Tmod.transpose(), Vmod)
        Werr = np.matmul(Tmod.transpose(), Wmod)
        Uerr -= Uf
        Verr -= Vf
        Werr -= Wf
        print('     Max. zonal error = ', abs(Uerr).flatten().max())
        print('Max. meridional error = ', abs(Verr).flatten().max())
        print('  Max. vertical error = ', abs(Werr).flatten().max())

    return PC, Um, Vm, Wm, Umod, Vmod, Wmod, Tmod


def spot_POD2F(U, V, wU, wV, nm, check_POD):
    "Snapshot POD of 2D velocity data"

    [nt, nu] = np.shape(U)
    # nv = np.shape(V[0, :])
    # Detrending
    Um = np.mean(U, axis=0)
    Vm = np.mean(V, axis=0)

    Uf = U - Um
    Vf = V - Vm

    # Build temporal correlation matrix
    C = inner_prod(Uf, Uf, wU)
    C += inner_prod(Vf, Vf, wV)
    C /= nt

    # Solve eigen problem
    eigval, eigvec = linalg.eigh(C)
    eigval = eigval.real
    eigvec = eigvec.real
    ids = np.argsort(eigval)            # sort eigvals
    ids = ids[::-1]                     # reverse order
    PC = eigval[ids]
    PC = PC[:nm]
    Tmod = eigvec[:, ids]   # sort temporal modes
    Tmod = Tmod[:,:nm].transpose()
    del C, eigval, eigvec

    # Scaling temporal modes
    for i in range(nm):
        Tmod[i, :] *= np.sqrt(nt*np.maximum(PC[i], 0.))

    # Build spatial modes
    Umod = np.matmul(Tmod, Uf)
    Vmod = np.matmul(Tmod, Vf)
    for i in range(nm):
        if PC[i] >= 1e-15:
            Umod[i, :] /= (nt*PC[i])
            Vmod[i, :] /= (nt*PC[i])
        else:
            Umod[i, :] *= 0.
            Vmod[i, :] *= 0.

    if check_POD:

        print('Check if the spatial modes are orthonormal')
        Cx = inner_prod(Umod, Umod, wU)
        Cx += inner_prod(Vmod, Vmod, wV)
        # plt.figure()
        # plt.imshow(Cx)
        # plt.colorbar()
        # plt.xlabel('modes')
        # plt.ylabel('modes')
        # plt.title('Orthogonality of spatial modes')
        # plt.show()

        print('Check if the temporal modes are orthogonal')
        Ct = np.matmul(Tmod, Tmod.transpose())
        for m in range(nm):
            Ct[m, m] /= (nt*PC[m])
        # plt.figure()
        # plt.imshow(Ct)
        # plt.colorbar()
        # plt.xlabel('modes')
        # plt.ylabel('modes')
        # plt.title('Orthogonality of temporal modes (rescaled by PCs)')
        # plt.show()

        print('Check POD reconstruction')
        Uerr = np.matmul(Tmod.transpose(), Umod)
        Verr = np.matmul(Tmod.transpose(), Vmod)
        Uerr -= Uf
        Verr -= Vf
        print('     Max. zonal error = ', abs(Uerr).flatten().max())
        print('Max. meridional error = ', abs(Verr).flatten().max())

    return PC, Um, Vm, Umod, Vmod, Tmod


def spot_POD1F(U, wU, nm, check_POD):
    "Snapshot POD of 2D velocity data"

    [nt, nu] = np.shape(U)
    # nv = np.shape(V[0, :])

    # Detrending
    Um = np.mean(U, axis=0)

    Uf = U - Um

    # Build temporal correlation matrix
    C = inner_prod(Uf, Uf, wU)
    C /= nt

    # Solve eigen problem
    eigval, eigvec = linalg.eigh(C)
    # eigval = eigval.real
    # eigvec = eigvec.real
    ids = np.argsort(eigval)            # sort eigvals
    ids = ids[::-1]                     # reverse order
    PC = eigval[ids]
    Tmod = eigvec[:, ids].transpose()   # sort temporal modes
    del C, eigval, eigvec

    # Scaling temporal modes
    for i in range(nm):
        Tmod[i, :] *= np.sqrt(nt*np.maximum(PC[i], 0.))

    # Build spatial modes
    Umod = np.matmul(Tmod, Uf)
    for i in range(nm):
        Umod[i, :] /= (nt*PC[i])

    if check_POD:

        print('Check if the spatial modes are orthonormal')
        Cx = inner_prod(Umod, Umod, wU)
        # plt.figure()
        # plt.imshow(Cx)
        # plt.colorbar()
        # plt.xlabel('modes')
        # plt.ylabel('modes')
        # plt.title('Orthogonality of spatial modes')
        # plt.show()

        print('Check if the temporal modes are orthogonal')
        Ct = np.matmul(Tmod, Tmod.transpose())
        for m in range(nm):
            Ct[m, m] /= (nt*PC[m])
        # plt.figure()
        # plt.imshow(Ct)
        # plt.colorbar()
        # plt.xlabel('modes')
        # plt.ylabel('modes')
        # plt.title('Orthogonality of temporal modes (rescaled by PCs)')
        # plt.show()

        print('Check POD reconstruction')
        Uerr = np.matmul(Tmod.transpose(), Umod)
        Uerr -= Uf
        print('     Max. zonal error = ', abs(Uerr).flatten().max())

    return PC, Um, Umod, Tmod


def inner_prod(U, V, W):
    ''' Inner product of U and V associated with weights W '''

    res = np.matmul(np.matmul(U, W), V.transpose())

    return res


def mat3D_reshape(A):
    ''' Reshaping a 3D matrix '''
    [nx, ny, nz] = np.shape(A)
    A = A.reshape((nx * ny * nz)) 

    return A


def mat2D_reshape(A):
    ''' Reshaping a 3D matrix '''
    [nx, ny] = np.shape(A)
    A = A.reshape((nx * ny)) 

    return A