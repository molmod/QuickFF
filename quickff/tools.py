# -*- coding: utf-8 -*-
# QuickFF is a code to quickly derive accurate force fields from ab initio input.
# Copyright (C) 2012 - 2015 Louis Vanduyfhuys <Louis.Vanduyfhuys@UGent.be>
# Steven Vandenbrande <Steven.Vandenbrande@UGent.be>,
# Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center for Molecular Modeling
# (CMM), Ghent University, Ghent, Belgium; all rights reserved unless otherwise
# stated.
#
# This file is part of QuickFF.
#
# QuickFF is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# QuickFF is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
#--

'''Various convenient tools.
'''

import numpy as np

from molmod.units import deg, angstrom
from molmod.periodic import periodic as pt
from molmod.molecular_graphs import HasNumNeighbors
from molmod.io import FCHKFile
from yaff import DihedCos, Chebychev1, Chebychev2, Chebychev3, Chebychev4, Chebychev6

from quickff.io import VASPRun

__all__ = [
    'global_translation', 'global_rotation', 'statistics', 'fitpar',
    'guess_ffatypes', 'get_vterm', 'read_abinitio', 'boxqp',
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
7
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


def guess_ffatypes(system, level):
    '''
       A method to guess atom types. This will overwrite ffatypes
       that are already defined in the system.

       **Arguments:**

       level
            A string used for guessing atom types:

                * low:     based on atomic number
                * medium:  based on atomic number and number of neighbors
                * high:    based on atomic number, number of neighbors and atomic number of neighbors
                * highest: based on index in the molecule
    '''
    if system.ffatypes is not None:
        raise ValueError('Atom types are already defined in the system.')
    if system.bonds is None:
        system.detect_bonds()
    if level == 'low':
        atypes = np.array([pt[number].symbol for number in system.numbers])
    elif level == 'medium':
        atypes = []
        for index, number in enumerate(system.numbers):
            nind = system.neighs1[index]
            sym = pt[system.numbers[index]].symbol.upper()
            atype = '%s%i' % (sym, len(nind))
            atypes.append(atype)
    elif level == 'high':
        atypes = []
        for index, number in enumerate(system.numbers):
            nind = system.neighs1[index]
            nsym = sorted([
                pt[system.numbers[neigh]].symbol.lower() for neigh in nind
            ])
            sym = pt[system.numbers[index]].symbol.upper()
            if len(nsym)==1:
                atype = '%s1_%s' % (sym, nsym[0])
            elif len(nsym)==2:
                atype = '%s2_%s%s' % (sym, nsym[0], nsym[1])
            else:
                atype = '%s%i' % (sym, len(nind))
                num_c = sum([1.0 for sym in nsym if sym == 'c'])
                num_n = sum([1.0 for sym in nsym if sym == 'n'])
                num_o = sum([1.0 for sym in nsym if sym == 'o'])
                if num_c > 0: atype += '_c%i' % num_c
                if num_n > 0: atype += '_n%i' % num_n
                if num_o > 0: atype += '_o%i' % num_o
            atypes.append(atype)
    elif level == 'highest':
        atypes = np.array([
            '%s%i' % (pt[n].symbol, i) for i, n in enumerate(system.numbers)
        ])
    else:
        raise ValueError('Invalid level, recieved %s' % level)
    system.ffatypes = atypes
    system._init_derived_ffatypes()
    _average_charges_ffatypes(system)

def _average_charges_ffatypes(system):
    'Take the average of the atomic charges and radii over the force field atom types'
    if system.charges is not None:
        for iffatype, ffatype in enumerate(system.ffatypes):
            mask = system.ffatype_ids==iffatype
            system.charges[mask] = np.mean(system.charges[mask])
    if system.radii is not None:
        for iffatype, ffatype in enumerate(system.ffatypes):
            mask = system.ffatype_ids==iffatype
            system.radii[mask] = np.mean(system.radii[mask])


def get_vterm( pars, indexes):
    # A torsion term with multiplicity m and rest value either 0 or pi/m
    # degrees, can be treated as a polynomial in cos(phi). The code below
    # selects the right polynomial.
    if pars[2] == 0.0 and pars[0] == 1:
        ic = DihedCos(*indexes)
        return Chebychev1(pars[1], ic, sign=-1)
    elif abs(pars[2] - np.pi/1)<1e-6 and pars[0] == 1:
        ic = DihedCos(*indexes)
        return Chebychev1(pars[1], ic, sign=1)
    elif pars[2] == 0.0 and pars[0] == 2:
        ic = DihedCos(*indexes)
        return Chebychev2(pars[1], ic, sign=-1)
    elif abs(pars[2] - np.pi/2)<1e-6 and pars[0] == 2:
        ic = DihedCos(*indexes)
        return Chebychev2(pars[1], ic, sign=1)
    elif pars[2] == 0.0 and pars[0] == 3:
        ic = DihedCos(*indexes)
        return Chebychev3(pars[1], ic, sign=-1)
    elif abs(pars[2] - np.pi/3)<1e-6 and pars[0] == 3:
        ic = DihedCos(*indexes)
        return Chebychev3(pars[1], ic, sign=1)
    elif pars[2] == 0.0 and pars[0] == 4:
        ic = DihedCos(*indexes)
        return Chebychev4(pars[1], ic, sign=-1)
    elif abs(pars[2] - np.pi/4)<1e-6 and pars[0] == 4:
        ic = DihedCos(*indexes)
        return Chebychev4(pars[1], ic, sign=1)
    elif pars[2] == 0.0 and pars[0] == 6:
        ic = DihedCos(*indexes)
        return Chebychev6(pars[1], ic, sign=-1)
    elif abs(pars[2] - np.pi/6)<1e-6 and pars[0] == 6:
        ic = DihedCos(*indexes)
        return Chebychev6(pars[1], ic, sign=1)
    else:
        return None


def read_abinitio(fn):
    '''Wrapper to read all information from an ab initio calculation that
    QuickFF needs. Currently Gaussian .fchk and VASP .xml files are supported.
    '''
    extension = fn.split('.')[-1]
    if extension=='fchk':
        fchk = FCHKFile(fn)
        numbers = fchk.fields.get('Atomic numbers')
        energy = fchk.fields.get('Total Energy')
        coords = fchk.fields.get('Current cartesian coordinates').reshape([len(numbers), 3])
        grad = fchk.fields.get('Cartesian Gradient').reshape([len(numbers), 3])
        hess = fchk.get_hessian().reshape([len(numbers), 3, len(numbers), 3])
        masses = None
        rvecs = None
        pbc = [0,0,0]
    elif extension=='xml':
        vasprun = VASPRun(fn,field_labels=['hessian','gradient'])
        numbers = vasprun.fields['numbers']
        coords = vasprun.fields['pos_init']
        energy = vasprun.fields['energies'][0]
        grad = vasprun.fields['gradient'][0]
        hess = vasprun.fields['hessian'].reshape((len(numbers),3,len(numbers),3 ))
        masses = vasprun.fields['masses']
        rvecs = vasprun.fields['rvecs_init']
        pbc = [1,1,1]
    else: raise NotImplementedError
    return numbers, coords, energy, grad, hess, masses, rvecs, pbc

def boxqp(A, B, bndl, bndu, x0, threshold=1e-12, status=False):
    '''Minimize the function
            1/2*xT.A.x - B.x
    subject to
            bndl < x < bndu (element-wise)
    This minimization is performed using a projected gradient method with
    step lengths computed using the Barzilai-Borwein method.
    See 10.1007/s00211-004-0569-y for a description.

    **Arguments**
        A       (n x n) NumPy array appearing in cost function
        B       (n) NumPy array appearing in cost function
        bndl    (n) NumPy array giving lower boundaries for the variables
        bndu    (n) NumPy array giving upper boundaries for the variables
        x0      (n) NumPy array providing an initial guess

    **Optional Arguments**
        threshold   Criterion to consider the iterations converged
        status      Return also the number of iterations performed
    '''
    # Check that boundaries make sense
    assert np.all(bndl<bndu), "Some lower boundaries are higher than upper boundaries"
    # Check that matrix A is positive definite
    def project(x):
        '''Project x on to the box of constraints'''
        x[x<bndl] = bndl[x<bndl]
        x[x>bndu] = bndu[x>bndu]
        return x
    def gradient(x):
        return np.dot(A,x) - B
    def stopping(x):
        q = gradient(x)
        mask = x==bndl
        q[mask] = np.amin(np.asarray([q[mask],[0.0]*np.sum(mask)]), axis=0)
        mask = x==bndu
        q[mask] = np.amax(np.asarray([q[mask],[0.0]*np.sum(mask)]), axis=0)
        return np.linalg.norm(q)
    # Bootstrapping alpha
    alpha = 0.1
    g0 = gradient(x0)
    x1 = project(x0-alpha*g0)
    gstop = np.linalg.norm(gradient(x1))
    converged = False
    nit = 0
    while converged is False:
        nit += 1
        # New gradient
        g1 = gradient(x1)
        # Compute new step length
        s = x1 - x0
        y = g1 - g0
        alpha = np.dot(s,s)/np.dot(s,y)
        # Update old values
        x0 = x1
        g0 = g1
        # Compute new values
        x1 = project(x1-alpha*g1)
        if stopping(x1)/gstop < threshold:
            converged = True
    if status: return x1, nit
    else: return x1
