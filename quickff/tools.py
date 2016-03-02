# -*- coding: utf-8 -*-
# QuickFF is a code to quickly derive accurate force fields from ab initio input.
# Copyright (C) 2012 - 2016 Louis Vanduyfhuys <Louis.Vanduyfhuys@UGent.be>
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

from molmod.units import deg, angstrom
from molmod.periodic import periodic as pt
from yaff import Chebychev1, Chebychev2, Chebychev3, Chebychev4, Chebychev6

import numpy as np, math

__all__ = [
    'global_translation', 'global_rotation', 'fitpar',
    'boxqp', 'guess_ffatypes', 'term_sort_atypes', 'get_multiplicity',
    'get_restvalue', 'get_ei_radii', 'digits'
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


def boxqp(A, B, bndl, bndu, x0, threshold=1e-9, status=False):
    '''
        Minimize the function
        
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


def guess_ffatypes(system, level):
    '''
       A method to guess atom types. This will overwrite ffatypes
       that are already defined in the system.

       **Arguments:**

       system
            A yaff system instance
       
       level
            A string used for guessing atom types:

                * low:     based on atomic number
                * medium:  based on atomic number and number of neighbors
                * high:    based on atomic number, number of neighbors and atomic number of neighbors
                * highest: based on index in the molecule
    '''
    if system.ffatypes is not None:
        raise ValueError('Atom types are already defined in the system.')
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
    system.ffatype_ids = np.zeros(len(system.numbers), int)
    system.ffatypes = []
    for i, atype in enumerate(atypes):
        if atype not in system.ffatypes:
            system.ffatypes.append(atype)
        system.ffatype_ids[i] = system.ffatypes.index(atype)
    system.ffatypes = np.array(system.ffatypes)


def term_sort_atypes(ffatypes, indexes, kind):
    '''
        Routine to sort the atoms defined in indexes to give consistent term
        names. This routine returns the sorted atom indexes as well as the 
        corresponding atom types.
    '''
    atypes = [ffatypes[i] for i in indexes]
    if kind in ['bond', 'dist', 'bend', 'angle']:
        if atypes[-1]<atypes[0] \
        or (atypes==atypes[::-1] and indexes[-1]<indexes[0]) :
            sorted_indexes = indexes[::-1]
            sorted_atypes = atypes[::-1]
        else:
            sorted_indexes = indexes
            sorted_atypes = atypes
    elif kind in ['dihed', 'dihedral', 'torsion']:
        if atypes[-1]<atypes[0] \
        or (atypes[-1]==atypes[0] and atypes[-2]<atypes[1]) \
        or (atypes==atypes[::-1] and indexes[-1]<indexes[0]):
            sorted_indexes = indexes[::-1]
            sorted_atypes = atypes[::-1]
        else:
            sorted_indexes = indexes
            sorted_atypes = atypes
    elif kind in ['opdist', 'oopdist']:
        pairs = sorted(zip(indexes[:3], atypes[:3]), key=lambda x: x[1]+str(x[0]))
        sorted_indexes = [index for index, atype in pairs]
        sorted_indexes.append(indexes[3])
        sorted_atypes = [atype for index, atype in pairs]
        sorted_atypes.append(atypes[3])
    return tuple(sorted_indexes), tuple(sorted_atypes)


def get_multiplicity(n1, n2):
    'Routine to estimate m from local topology'
    if   set([n1,n2])==set([4,4]): return 3
    elif set([n1,n2])==set([3,4]): return 6
    elif set([n1,n2])==set([2,4]): return 3
    elif set([n1,n2])==set([3,3]): return 2
    elif set([n1,n2])==set([2,3]): return 2
    elif set([n1,n2])==set([2,2]): return 1
    else:                          return None


def get_restvalue(values, m, thresshold=5*deg):
    '''
        Get a rest value of 0.0, 360/(2*m) or None depending on the given
        equilbrium values
    '''
    rv = None
    per = 360*deg/m
    for value in values:
        x = value % per
        if abs(x)<=thresshold or abs(per-x)<thresshold:
            if rv is not None and rv!=0.0:
                return None
            elif rv is None:
                rv = 0.0
        elif abs(x-per/2.0)<thresshold:
            if rv is not None and rv!=per/2.0:
                return None
            elif rv is None:
                rv = per/2.0
        else:
            return None
    return rv


def get_ei_radii(numbers):
    '''
        Routine to return atomic radii for use in the Gaussian charge
        distribution. These radii are computed according to the procedure of
        Chen et al.:
        
        First the Slater exponent is computed from the hardness using the
        formula of Rappe and Goddard (hardness of Pearson and Parr is used)
        
        Next the gaussian exponent alpha is fitted by minimizing the
        L2-difference between the between the homonuclear Coulomb integral over
        Slater orbitals and over Gaussian orbitals.
    '''
    radii = {
        'H' : 0.7309*angstrom, 
        'Li': 1.2951*angstrom, 'B' : 1.2020*angstrom, 'C' : 1.1646*angstrom,
        'N' : 1.1039*angstrom, 'O' : 1.1325*angstrom, 'F' : 1.1097*angstrom, 
        'Na': 1.7093*angstrom, 'Al': 1.6744*angstrom, 'Si': 1.6378*angstrom,
        'P' : 1.5730*angstrom, 'S' : 1.6022*angstrom, 'Cl': 1.5789*angstrom,
    }
    values = np.zeros(len(numbers), float)
    for i, number in enumerate(numbers):
        symbol = pt[number].symbol
        if not symbol in radii.keys():
            raise ValueError('No electrostatic Gaussian radii found for %s' %symbol)
        values[i] = radii[symbol]
    return values


def digits(x,n):
    """
        returns a string representation of x with exactly n digits if possible.
    """
    if len(str(x))==n: return str(x)
    sign = np.sign(x)
    x = float(abs(x))
    if np.isnan(x):
        return ''
    if abs(x) < 1e-6:
        return "." + "0"*(n-1)
    if sign<1: n -= 1
    i = int(x)
    r = x-i
    if i==0:
        if sign<0:
            return '-'+str(r).lstrip('0')[:n]
        else:
            return str(r).lstrip('0')[:n]
    if r==0:
        return str(int(i*sign))[:n]
    if len(str(i))>=(n-1):
        return str(int(i*sign))
    ndig = n - len(str(i))-1
    if sign<0:
        return '-%i.%s' %(i, str(r).lstrip('0.')[:ndig])
    else:
        return '%i.%s' %(i, str(r).lstrip('0.')[:ndig])
