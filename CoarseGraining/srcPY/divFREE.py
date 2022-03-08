#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 22:04:49 2022

@author: ftucciar
"""
import sys
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as linalg
from netCDF4 import Dataset


def Stiffness_Assembly(x, y, triang, D):
    """

      Subroutine: Stiffness_Assembly

         Purpose: Assembles the stiffness matrix of the problem

      INPUT variables

         x,y         : Nodal coordinates
         triang      : Topological set of the mesh
         D           : Diffusion matrix

      Internal Variables

         a,b,c       : Coefficients of the interpolant linear polynomial

                           a_i = x_j*y_k - x_k*y_j
                           b_i = y_j - y_k
                           c_i = x_k - x_j

         Area_e      : Area of the element
         H_e         : Local (3x3) stiffness matrix

      OUTPUT variables

         matK        : Stiffness Matrix

    """
#   Set initial parameters ----------------------------------------------
    ne = np.size(triang, 0)
    n = np.size(x, 0)
    Area_e = np.zeros([ne])

#   Initialize Stiffness matrix -----------------------------------------
    rows, cols = n, n
    matK = sps.coo_matrix((rows, cols))

#   Looping the elements ------------------------------------------------
    for iel in range(ne):

        # Inizialization ------------------------------------------------
        x_e = np.zeros(3)
        y_e = np.zeros(3)
        a = np.zeros(3)
        b = np.zeros(3)
        c = np.zeros(3)
        H_e = np.zeros([3, 3])

        # Inquiring topological pattern ---------------------------------
        x_e = x[triang[iel, :]]
        y_e = y[triang[iel, :]]

        # Coefficients of the interpolant -------------------------------
        a[0] = x_e[1]*y_e[2] - x_e[2]*y_e[1]
        a[1] = x_e[2]*y_e[0] - x_e[0]*y_e[2]
        a[2] = x_e[0]*y_e[1] - x_e[1]*y_e[0]
        b[0] = y_e[1] - y_e[2]
        b[1] = y_e[2] - y_e[0]
        b[2] = y_e[0] - y_e[1]
        c[0] = x_e[2] - x_e[1]
        c[1] = x_e[0] - x_e[2]
        c[2] = x_e[1] - x_e[0]

        # Area of the element -------------------------------------------
        Area_e[iel] = (a[0] + (x_e[0] * b[0]) + (y_e[0] * c[0])) / 2

        # Local stiffness matrix build ----------------------------------
        H_e = (D[0][0] * np.outer(b, b) +
               D[1][1] * np.outer(c, c))/(4 * Area_e[iel])

        #  Sparsity pattern of the matrix -------------------------------
        row = np.repeat(triang[iel, :], 3)
        col = np.tile(triang[iel, :], 3)
        H_e = H_e.reshape(3*len(H_e))

        # Global stiffness matrix assembly ------------------------------
        matK = matK + sps.coo_matrix((H_e, (row, col)), shape=(n, n))

    return matK


def RHS_Assembly(x, y, triang, rhs):
    """

      Subroutine: RHS_Assembly

         Purpose: Assembles the right-hand side for the problem

      INPUT variables

         x,y         : Nodal coordinates
         triang      : Topological set of the mesh
         rhs         : External forcing function

      OUTPUT variables

         f           : Right hand side
         areanod     : Afference area of each node

    """
#   Set initial parameters ----------------------------------------------
    ne = np.size(triang, 0)
    n = np.size(x, 0)

#   Evaluation of the afference area ------------------------------------
    areanod = np.zeros([n])

    for k in range(ne):
        sort = np.sort(triang[k, :])
        C = np.array([[1, x[sort[0]], y[sort[0]]],
                      [1, x[sort[1]], y[sort[1]]],
                      [1, x[sort[2]], y[sort[2]]]])
        area = 0.5*np.abs(np.linalg.det(C))
        areanod[sort] = areanod[sort] + area/3

#   Evaluation of the afference area ------------------------------------
    f = np.multiply(rhs, areanod)

    return f, areanod


def DirPen(nodes, values, Stiff, rhs, penalty=10e5):
    """

      Subroutine: DirPen

         Purpose: Applies the DIRICHLET boundary conditions with
                  penalty approach

      INPUT variables

         nodes       : Boundary nodes
         values      : Values to impose
         Stiff       : Stiffness matrix of the problem
         rhs         : Right hand side of the problem
         penalty     : Value for the penalty coefficient

      OUTPUT variables

         f           : Updated right hand side
         Stiff       : Updated stiffness matrix

    """
#   Application to the right hand side ----------------------------------
    rhs[nodes] = values*penalty

#   Application to the Stiffness Matrix ---------------------------------
    diag = Stiff.diagonal()
    diag[nodes] = penalty
    Stiff.setdiag(diag)

    return Stiff, f


# mesh Inputs -----------------------------------------------------------
IDmesh = 1000
path = ''
meshName = 'mesh'+str(IDmesh)
fullpath = path

bound = np.loadtxt(fullpath + meshName + '/dirnod.dat').astype(int)
topol = np.loadtxt(fullpath + meshName + '/mesh.dat').astype(int)
coord = np.loadtxt(fullpath + meshName + '/xy.dat')

print('Loading informations from:')
print(fullpath)

# Preparing the data: The dataset was thought to be used in FORTRAN95
#   and/or MATLAB so the elements are indexed starting from 1 and not
#   0 (C,C++,Python standard), hence a little preparing is needed, by
#   scaling the data .topol and .bound all by one
topol = topol-np.ones([np.size(topol, 0), np.size(topol, 1)]).astype(int)
bound = bound-np.ones([np.size(bound, 0)]).astype(int)
dirbnd = np.zeros(np.size(bound))

# Phisiscal parameter of the problem
D = [[1.00, 0.00],
     [0.00, 1.00]]

# Numerical parameters
Penalty = 10e15

# netCDF Inputs ---------------------------------------------------------
base_dir = ''  # '/Volumes/LaCie/Nemo/Data_Tuccia/R27/'
subs_dir = ['100-102y/', '102-104y/', '104-106y/', '106-108y/', '108-110y/']
infile = 'ocref_r3.nc'
ingrid = 'domain_cfg_out.nc'

# Outputs ---------------------------------------------------------------
sfx = ['u', 'v', 'w']
outfile = base_dir + 'trial.nc'


# Opening netCDF files --------------------------------------------------
fin = Dataset(base_dir + subs_dir[0] + infile, 'r')
grid = Dataset(ingrid, 'r')

# Get domain dimension from domain file ---------------------------------
nx = grid.variables['jpiglo'][0]
ny = grid.variables['jpjglo'][0]
nz = grid.variables['jpkglo'][0]

# Evaluation of the number of times instants ----------------------------
nt = 0
for s in subs_dir:
    file1 = base_dir + s + infile
    fin = Dataset(file1, 'r')
    tmp = fin.dimensions['time_counter'].size
    nt += tmp
    fin.close()

# Save output data
# ----------------
print('Preparing output file: ', outfile)
fout = Dataset(outfile, 'w', format='NETCDF4')
# Create dimensions
fout.createDimension('time', nt)
fout.createDimension('z', nz)
fout.createDimension('y', ny)
fout.createDimension('x', nx)
# Create variables
udf = fout.createVariable('udf', 'f4', ('time', 'z', 'y', 'x',))
vdf = fout.createVariable('vdf', 'f4', ('time', 'z', 'y', 'x',))
# Add attributes
udf.long_name = '2D divfree zonal velocity'
vdf.long_name = '2D divfree meridional velocity'
udf.units = 'm/s'
vdf.units = 'm/s'

# U-points (for integration) --------------------------------------------
e1u = grid.variables['e1u'][0, :, :]
e2u = grid.variables['e2u'][0, :, :]
# V-points (for integration) --------------------------------------------
e1v = grid.variables['e1v'][0, :, :]
e2v = grid.variables['e2v'][0, :, :]
# W-points (for integration) --------------------------------------------
e1w = grid.variables['e1t'][0, :, :]
e2w = grid.variables['e2t'][0, :, :]
# f-points (for integration) --------------------------------------------
e1f = grid.variables['e1f'][0, :, :]
e2f = grid.variables['e2f'][0, :, :]
v3 = np.multiply(e1f, e2f)
# Mid-point rule area at u and v points ---------------------------------
dAu = np.multiply(e1u, e2u)
dAv = np.multiply(e1v, e2v)
dAw = np.multiply(e1w, e2w)
# Free memory
del e1w, e2w

# Allocation of empty variables -----------------------------------------
# ----------------------------
wpu = np.diag(dAu.reshape((nx * ny)))
wpv = np.diag(dAv.reshape((nx * ny)))
wpw = np.diag(dAw.reshape((nx * ny)))
# Free memory -----------------------------------------------------------
del dAu, dAv, dAw

u = np.empty((nt, nz, ny, nx))
v = np.empty((nt, nz, ny, nx))
curl = np.empty((nt, nz, ny, nx))

# POD procedure layer by layer
# ----------------------------
for k in range(nz - 1):
    print()
    print('Layer = ', k+1)
    print('--------')
    print('Collecting data')
    # Allocate temporary matrices
    u = np.empty((0, ny, nx))
    v = np.empty((0, ny, nx))
    # Append velocity data from every .nc file
    for s in subs_dir:
        file1 = base_dir + s + infile
        print('Opening file: ', file1)
        fin = Dataset(file1, 'r')
        tmp = fin.variables['ur'][:, k, :, :].copy()
        u = np.append(u, tmp, axis=0)
        tmp = fin.variables['vr'][:, k, :, :].copy()
        v = np.append(v, tmp, axis=0)
        del tmp
        fin.close()
    print('Compute curl of 2D field')
    v1 = np.multiply(e2v, v)
    v2 = np.multiply(e1u, u)
    curl[:, k, :-1, :-1] = (v1[:, :-1, 1:] - v1[:, :-1, :-1]) - \
        (v2[:, 1:, :-1] - v2[:, :-1, :-1])
    curl[:, k, :, :] = np.divide(curl[:, k, :, :], v3)
    del v1, v2
    print('Inverting the laplacian to get the streamfunction')
    # Building the Stiffness Matrix
    H = Stiffness_Assembly(coord[:, 0], coord[:, 1], topol, D)
    for t in range(nt):

        # Building the Right Hand Side
        values = curl[t, k, :, :].flatten('C')
        f, a = RHS_Assembly(coord[:, 0], coord[:, 1], topol, values)

        # Imposing the boundary conditions
        H, f = DirPen(bound, dirbnd, H, f)

        # Solving the linear system
        sol = linalg.cg(H, f)[0]

        # Reshape solution
        psimat = np.reshape(sol, (ny, nx), order='C')

        # Compute the gradient
        [dxpsi, dypsi] = np.gradient(psimat)

        # Compute the velocity field
        udf[t, k, 1:, :-1] = np.divide(psimat[1:, :-1] -
                                       psimat[:-1, :-1], e2u[1:, :-1])
        vdf[t, k, 1:, :-1] = -np.divide(psimat[:-1, 1:] -
                                        psimat[:-1, :-1], e1v[:-1, 1:])

        print(t)

# Close file
fout.close()