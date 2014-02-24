# -*- coding: utf-8 -*-
#QuickFF is a code to quickly derive accurate force fields from ab initio input.
#Copyright (C) 2012 - 2014 Louis Vanduyfhuys <Louis.Vanduyfhuys@UGent.be>
#Steven Vandenbrande <Steven.Vandenbrande@UGent.be>, 
#Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center for Molecular Modeling
#(CMM), Ghent University, Ghent, Belgium; all rights reserved unless otherwise
#stated.
#
#This file is part of QuickFF.
#
#QuickFF is free software; you can redistribute it and/or
#modify it under the terms of the GNU General Public License
#as published by the Free Software Foundation; either version 3
#of the License, or (at your option) any later version.
#
#QuickFF is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, see <http://www.gnu.org/licenses/>
#
#--

import numpy as np
from molmod.units import deg
from molmod.molecular_graphs import HasNumNeighbors

__all__ = [
    'global_translation', 'global_rotation', 'calc_angles', 'statistics',
    'fitpar', 'matrix_squared_sum', 'find_opdist_patterns',
]

def global_translation(coords):
    '''
        A function to generate vectors that represent global translations
        of a system.

        **Arguments**

        coords
            a (N,3) numpy array describing the system that has to be translated
    '''
    Natoms = len(coords)
    ones = np.ones(Natoms, float)
    zeros = np.zeros(Natoms, float)
    VTx = np.concatenate(np.array([ones, zeros, zeros]).transpose())/np.sqrt(Natoms)
    VTy = np.concatenate(np.array([zeros, ones, zeros]).transpose())/np.sqrt(Natoms)
    VTz = np.concatenate(np.array([zeros, zeros, ones]).transpose())/np.sqrt(Natoms)
    return VTx, VTy, VTz


def global_rotation(coords):
    '''
        A function to generate vectors that represent global translations
        of a system. Rx is a matrix of rotatino around the x-axis minus
        the identity matrix. VRx is a vector of rotation around x-axis.

        **Arguments**

        coords
            a (N,3) numpy array describing the system that has to be translated
    '''
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
    U = np.linalg.svd( np.array([VRx, VRy, VRz]).transpose() )[0]
    VRx = U.transpose()[0]
    VRy = U.transpose()[1]
    VRz = U.transpose()[2]
    return VRx, VRy, VRz


def calc_angles(v, refs):
    '''
        Function to calculate the angles between v and each of the refs

        **Arguments**

        v
            a (D) numpy array defining a D-dimensional vector

        ref
            a (N, D) numpy array containing N D-dimensional vectors
    '''
    angles = []
    for ref in refs:
        cangle = np.dot(ref, v)/(np.linalg.norm(ref)*np.linalg.norm(v))
        if cangle < -1.0:
            cangle = -1.0
        if cangle > 1.0:
            cangle = 1.0
        angles.append( np.arccos(cangle) )
    return angles


def statistics(data):
    '''
        Function to calculate mean, standard deviation and dimension of data

        **Arguments**

        data
            a (N) numpy array containing the statistical input data
    '''
    if data is None:
        return None, None, None
    else:
        N = len(data)
        if N == 0:
            return None, None, None
        elif N == 1:
            return data[0], 0.0, 1
        elif N > 1:
            mean = data.sum()/N
            std = np.sqrt( ((data-mean)**2).sum()/(N-1) )
            return mean, std, N


def fitpar(xs, ys, rcond=1e-3):
    '''
        Fit a parabola to the samples (xs, ys):

            ys[:] = a*xs[:]^2 + b*xs[:] + c

        Returns the parabola parameters a, b and c.

        **Arguments**

        xs
            a (N) numpy array containing the x values of the samples

        ys
            a (N) numpy array containing the x values of the samples
    '''
    assert len(xs)==len(ys)
    D = np.ones([len(xs), 3], float)
    for i, x in enumerate(xs):
        D[i, 0] = x**2
        D[i, 1] = x
    sol = np.linalg.lstsq(D, ys, rcond=rcond)[0]
    return sol

def matrix_squared_sum(A, B):
    '''
        Calculate the sum of the product of all matrix elements

            sum(Aij*Bij, i, j)

        **Arguments**

        A
            a (M,N) numpy array

        B
            a (M,N) numpy array
    '''
    tmp = np.dot(A.T, B)
    dim = len(tmp)
    return sum([tmp[i, i] for i in xrange(dim)])

def find_opdist_patterns(graph):
    '''
        Find patterns of 4 atoms where 3 border atoms are bonded to the same
        central atom. The central atom cannot have any other neighbor except
        for these 3 border atoms. Returns a list of 4-tuples in which the first
        three atoms are the border atoms and the last atom is the central atom.

        **Arguments**

        graph
            An instance of the MolecularGraph class from the MolMod package.
            More info on this class can be found in the
            `MolMod documentation <http://molmod.github.io/molmod/reference/basic.html#molmod-molecular-graphs-molecular-graphs>`_.
    '''
    opdists = []
    for atom in xrange(len(graph.numbers)):
        if HasNumNeighbors(3)(atom, graph):
            neighs = tuple(graph.neighbors[atom])
            opdists.append([neighs[0], neighs[1], neighs[2], atom])
    return opdists
