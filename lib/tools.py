#! /usr/bin/env python

import numpy as np
from molmod.units import parse_unit

__all__ = [
    'global_translation', 'global_rotation', 'calc_angles', 'statistics',
    'fitpar', 'has_15_bonded', 'get_atoms_within_3bonds', 'matrix_squared_sum'
]

def global_translation(coords):
    Natoms = len(coords)
    ones = np.ones(Natoms, float)
    zeros = np.zeros(Natoms, float)
    VTx = np.concatenate(np.array([ones, zeros, zeros]).transpose())/np.sqrt(Natoms)
    VTy = np.concatenate(np.array([zeros, ones, zeros]).transpose())/np.sqrt(Natoms)
    VTz = np.concatenate(np.array([zeros, zeros, ones]).transpose())/np.sqrt(Natoms)
    return VTx, VTy, VTz


def global_rotation(coords):
    #rotations (Rx = matrix of Rotation around x-axis minus the identity matrix,
    #           VRx = Vector of Rotation in x direction)
    Natoms = len(coords)
    com = coords.sum(axis=0)/coords.shape[0]
    Rz = np.array([
        [ 0.0,-1.0, 0.0],
        [ 1.0, 0.0, 0.0],
        [ 0.0, 0.0, 0.0]
    ])
    Ry = np.array([
        [ 0.0, 0.0, 1.0],
        [ 0.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0]
    ])
    Rx = np.array([
        [ 0.0, 0.0, 0.0],
        [ 0.0, 0.0, 1.0],
        [ 0.0,-1.0, 0.0]
    ])
    VRx = np.dot(coords-com, Rx.transpose()).reshape([3*Natoms])
    VRy = np.dot(coords-com, Ry.transpose()).reshape([3*Natoms])
    VRz = np.dot(coords-com, Rz.transpose()).reshape([3*Natoms])
    U, S, Vt = np.linalg.svd( np.array([VRx, VRy, VRz]).transpose() )
    VRx = U.transpose()[0]
    VRy = U.transpose()[1]
    VRz = U.transpose()[2]
    return VRx, VRy, VRz


def calc_angles(v, refs):
    angles = []
    for ref in refs:
        cangle = np.dot(ref, v)/(np.linalg.norm(ref)*np.linalg.norm(v))
        if cangle<-1.0: cangle=-1.0
        if cangle>1.0: cangle=1.0
        angles.append( np.arccos(cangle) )
    return angles


def statistics(data):
    if data is None:
        return None, None, None
    else:
        N = len(data)
        if N==0:
            return None, None, None
        elif N==1:
            return data[0], 0.0, 1
        elif N>1:
            mean = data.sum()/N
            std = np.sqrt( ((data-mean)**2).sum()/(N-1) )
            return mean, std, N


def fitpar(xs, ys, rcond=1e-3):
    '''
        Fit a parabola to xs and ys:

            ys[:] = a*xs[:]^2 + b*xs[:] + c

        Returns a, b and c.
    '''
    assert len(xs)==len(ys)
    D = np.ones([len(xs),3], float)
    for i, x in enumerate(xs):
        D[i,0] = x**2
        D[i,1] = x
    sol, res, rank, sing = np.linalg.lstsq(D, ys, rcond=rcond)
    return sol

def has_15_bonded(system):
    assert hasattr(system, 'nlist')
    if not hasattr(system, 'diheds'):
        return False
    for dihed in system.diheds:
        neigh0 = system.nlist[dihed[0]]
        neigh3 = system.nlist[dihed[3]]
        if len(neigh0)>1 or len(neigh3)>1:
            return True
    return False

def get_atoms_within_3bonds(system):
    pairs = []
    for bond in system.bonds:
        pairs.append([bond[0], bond[1]])
    for bend in system.bends:
        pairs.append([bend[0], bend[2]])
    for dihed in system.diheds:
        pairs.append([dihed[0], dihed[3]])
    return np.array(pairs)

def matrix_squared_sum(A, B):
    '''
        Calculate the sum of the product of all matrix elements

            sum(Aij*Bij, i, j)
    '''
    tmp = np.dot(A.T, B)
    dim = len(tmp)
    return sum([tmp[i,i] for i in xrange(dim)])
