# -*- coding: utf-8 -*-
# QuickFF is a code to quickly derive accurate force fields from ab initio input.
# Copyright (C) 2012 - 2019 Louis Vanduyfhuys <Louis.Vanduyfhuys@UGent.be>
# Steven Vandenbrande <Steven.Vandenbrande@UGent.be>,
# Jelle Wieme <Jelle.Wieme@UGent.be>,
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

from __future__ import print_function, absolute_import

from molmod.units import deg, angstrom, centimeter
from molmod.constants import lightspeed
from molmod.periodic import periodic as pt

from yaff import Chebychev1, Chebychev2, Chebychev3, Chebychev4, Chebychev6

from quickff.log import log

import numpy as np, math

__all__ = [
    'global_translation', 'global_rotation', 'fitpar',
    'boxqp', 'set_ffatypes', 'term_sort_atypes', 'get_multiplicity',
    'get_restvalue', 'get_ei_radii', 'digits', 'average', 'chebychev',
    'project_negative_freqs'
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


def fitpar(xs, ys, rcond=-1):
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
        D[i, 0] = 0.1*x**2
        D[i, 1] = x
    sol, res, rank, svals = np.linalg.lstsq(D, ys, rcond=rcond)
    sol[0] *= 0.1
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


def set_ffatypes(system, level):
    '''
       A method to guess atom types. This will overwrite ffatypes
       that are already defined in the system.

       **Arguments:**

       system
            A yaff system instance

       level
            If level is a string containing comma's, it is assumed to be an
            ordered list containing the atom type of each atom in the system.
            Otherwise, level is assumed to be a string defining how to guess
            the atom types from the local topology. Possible levels are:

                * low:     based on atomic number
                * medium:  based on atomic number and number of neighbors
                * high:    based on atomic number, number of neighbors and atomic number of neighbors
                * highest: based on index in the molecule

    '''
    if system.ffatypes is not None:
        raise ValueError('Atom types are already defined in the system.')
    if ',' in level:
        atypes = level.split(',')
    elif level == 'low':
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
            nsyms = sorted([
                pt[system.numbers[neigh]].symbol.lower() for neigh in system.neighs1[index]
            ])
            sym = pt[system.numbers[index]].symbol.upper()
            if len(nsyms)==1:
                atype = '%s1_%s' % (sym, nsyms[0])
            elif len(nsyms)==2:
                atype = '%s2_%s%s' % (sym, nsyms[0], nsyms[1])
            else:
                atype = '%s%i' % (sym, len(system.neighs1[index]))
                neighs = {}
                for nsym in nsyms:
                    if nsym=='h': continue
                    if nsym in list(neighs.keys()):
                        neighs[nsym] += 1
                    else:
                        neighs[nsym] = 1
                for nsym, nnum in neighs.items():
                    atype += '_%s%i' %(nsym, nnum)
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
    else:                          return np.nan


def get_restvalue(values, m, thresshold=20*deg, mode=1):
    '''
        Get a rest value of 0.0, 360/(2*m) or None depending on the given
        equilbrium values.

        For mode=0, the rest value is:

            0, if all 'values modulo per' are in the interval
            [0,thresshold] U [per-thresshold,per]

            per/2, if all 'values modulo per' are in the interval
            [per/2-thresshold,per/2+thresshold]

            None, in all other cases

        For mode=1, the rest value is determined as follows:

            first the values are folded in the interval [0,per/2] by first
            taking the module with per and then mirroring values in [per/2,per]
            on [0,per/2]. Next, the mean and std of the folded values are
            computed. If the std is larger then the thresshold, the values are
            considered to be too scattered and no rest value can be computed
            (None is returned). If the std is small enough, the rest value will
            be determined based on the mean. If the mean is close enough to 0,
            the rest value will be 0. If the mean is close enough to per/2, the
            rest value will be per/2. In all other cases no rest value will be
            computed and None is returned.

    '''
    rv = None
    per = 360*deg/m
    if mode==0:
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
    elif mode==1:
        folded = np.zeros(len(values), float)
        for i, value in enumerate(values):
            new = value % per
            if new>0.5*per: folded[i] = per - new
            else:           folded[i] = new
        mean = folded.mean()
        std = folded.std()
        assert 0.0<=mean and mean<=0.5*per
        if std<thresshold:
            if mean<=0.25*per: rv=0.0
            elif 0.5*per-mean<=0.25*per: rv=0.5*per
        return rv
    else:
        raise NotImplementedError('Mode %i in get_restvalue is not supported' %mode)


def get_ei_radii(numbers):
    '''
        Routine to return atomic radii for use in the Gaussian charge
        distribution. These radii are computed according to the procedure of
        Chen and Slater:

        First the Slater exponent is computed from the hardness using the
        formula of Rappe and Goddard (hardness of Pearson and Parr is used)

        Next the gaussian exponent alpha is fitted by minimizing the
        L2-difference between the between the homonuclear Coulomb integral over
        Slater orbitals and over Gaussian orbitals.
    '''
    radii = {
        'H' : 0.7308*angstrom,
        'Li': 1.2951*angstrom, 'B' : 1.2020*angstrom, 'C' : 1.1703*angstrom,
        'N' : 1.1048*angstrom, 'O' : 1.1325*angstrom, 'F' : 1.1096*angstrom,
        'Na': 1.7093*angstrom, 'Mg': 1.6155*angstrom, 'Al': 1.6742*angstrom,
        'Si': 1.6376*angstrom, 'P' : 1.5727*angstrom, 'S' : 1.6011*angstrom,
        'Cl': 1.5798*angstrom, 'Ca': 1.6541*angstrom, 'Sc': 2.0559*angstrom,
        'Ti': 2.0502*angstrom, 'V' : 2.0654*angstrom, 'Cr': 2.0692*angstrom,
        'Mn': 2.0323*angstrom, 'Fe': 2.0695*angstrom, 'Co': 2.0377*angstrom,
        'Ni': 2.0579*angstrom, 'Cu': 2.0573*angstrom, 'Zn': 1.9896*angstrom,
        'Ga': 2.0820*angstrom, 'Br': 2.0088*angstrom,
    }

    values = np.zeros(len(numbers), float)
    for i, number in enumerate(numbers):
        symbol = pt[number].symbol
        if not symbol in list(radii.keys()):
            raise ValueError('No electrostatic Gaussian radii found for %s' %symbol)
        values[i] = radii[symbol]
    return values


def digits(x,n):
    """
        returns a string representation of x with exactly n digits if possible.
    """
    if np.isnan(x): return ''
    if len(str(x))==n: return str(x)
    sign = np.sign(x)
    x = float(abs(x))
    if x < 0.5*10**(-(n-1)):
        return "." + "0"*(n-1)
    if sign<0: n -= 1
    i = int(x)
    r = x-int(x)
    if i==0:
        if sign<0:
            return '-'+str(r)[1:1+n]
        else:
            s = '%f' %r
            return s[1:1+n]
    if r==0:
        return str(int(i*sign))[:n]
    if len(str(i))>=(n-1):
        return str(int(i*sign))
    ndig = n - len(str(i))-1
    if sign<0:
        return '-%i.%s' %(i, str(r)[2:2+ndig])
    else:
        return '%i.%s' %(i, str(r)[2:2+ndig])

def average(data, ffatypes, fmt='full', verbose=False):
    '''
        Average the atomic parameters stored in data over atoms of the same atom
        type.

        **Arguments**

        data
            a list or numpy array containing the data

        ffatypes
            a listor numpy array containing the atom types. Should have equal
            length as data

        **Keywork arguments**

        fmt
            Should be either full, dict or sort. In case of full, the result will be
            returned as an numpy array of equal length as data and ffatypes. In
            case of dict, the result will be returned as a dictionairy of the
            following format:

                {ffatype0: value0, ffatype1: value1, ...}

            in which value0, ... is the mean value of the given ffatype. In case
            of sort, a dictionairy of the following format will be returned:

                {ffatype0: values0, ffatype1: values1, ...}

            in which values0, ... is a list of the values for the given ffatype.
    '''
    data_atypes = {}
    for value, ffatype in zip(data,ffatypes):
        if ffatype in list(data_atypes.keys()):
            data_atypes[ffatype].append(value)
        else:
            data_atypes[ffatype] = [value]
    if fmt=='sort':
        output = {}
        for ffatype, data in data_atypes.items():
            output[ffatype] = np.array(data)
    elif fmt=='full':
        output = np.zeros(len(data))
        printed = []
        for i, ffatype in enumerate(ffatypes):
            std = np.array(data_atypes[ffatype]).std()
            if not std < 1e-2 and ffatype not in printed:
                print('WARNING: charge of atom type %s has a large std: %.3e' %(ffatype, std))
                printed.append(ffatype)
            output[i] = np.array(data_atypes[ffatype]).mean()
    elif fmt=='dict':
        output = {}
        for ffatype, values in data_atypes.items():
            output[ffatype] = np.array(values).mean()
    else:
        raise IOError('Format %s not supported, should be full or dict' %fmt)
    if verbose:
        print('Averaged Atomic Charges:')
        print('------------------------')
        for ffatype, values in data_atypes.items():
            print('  %4s    % .3f +- % .3f (N=%i)' %(ffatype, np.array(values).mean(), np.array(values).std(), len(values)))
        print('')
    return output

def charges_to_bcis(charges, ffatypes, bonds, constraints={}, verbose=True):
    '''
        Transform atomic charges to bond charge increments, by definition 2
        bci's between different pairs of atoms but with identical pairs of
        atom types will be equal. Bci's will be returned as a dictionairy
        containing tuples of the format (ffatype0.ffatype1, bci_value)

        **Arguments**

        charges
            a (N,) list or numpy array containing the charges

        ffatypes
            a (N,) list or numpy array with the atom type of each atom in the
            system

        bonds
            a (B,2) list or numpy array for each bond in the system

        **Keyword Arguments**

        constraints
            a dictionairy of format (master, [(slave0,sign0), (slave1,sign1), ...])

        verbose
            increase verbosity
    '''
    assert len(charges)==len(ffatypes)
    #construct list of bond types and signs, the signs are related to the
    #direction of the bci of a certain bond. We want that bonds of type
    #A.B have a bci with equal magnitude as a bond of type B.A. Therefore, we
    #only store the bci of A.B (alphabetically) and also store a sign, which is
    #1.0 for bond A.B and -1 for bond B.A
    btypes = ['',]*len(bonds)
    signs = np.zeros([len(bonds)], float)
    for i, bond in enumerate(bonds):
        ffatype0, ffatype1 = ffatypes[bond[0]], ffatypes[bond[1]]
        if ffatype0<ffatype1:
            btypes[i] = '%s.%s' %(ffatype0,ffatype1)
            signs[i] = 1.0
        else:
            btypes[i] = '%s.%s' %(ffatype1,ffatype0)
            signs[i] = -1.0
    #decompile constraints
    masterof = {}
    for m, s in constraints.items():
        types = m.split('.')
        if types[0]>=types[1]: m = '.'.join(types[::-1])
        for slave, sign in s:
            types = slave.split('.')
            if types[0]>=types[1]: slave = '.'.join(types[::-1])
            assert slave not in list(masterof.keys()), \
                'Slave %s has multiple masters in constraints' %slave
            masterof[slave] = (m, sign)
    masterlist = []
    for btype in btypes:
        if btype in list(masterof.keys()) : continue
        if btype in masterlist: continue
        masterlist.append(btype)
    for master in masterlist:
        assert not master in list(masterof.keys()), 'master %s encountered in slaves' %master
    for slave in list(masterof.keys()):
        assert not slave in masterlist, 'slave %s encountered in masters' %slave
    if verbose:
        print('Master-slaves relations:')
        print('------------------------')
        for slave in list(masterof.keys()):
            print(slave, masterof[slave])
        if len(list(masterof.keys()))==0:
            print('(None)')
        print('')
    #construct the matrix to convert bci's to charges
    #matrix[i,n] is the contribution to charge i from bci n
    #bci p_AB is a charge transfer from B to A, hence qA+=p_AB and qB-=p_AB
    matrix = np.zeros([len(charges), len(masterlist)], float)
    for i, (btype, bond, sign) in enumerate(zip(btypes, bonds, signs)):
        if btype in masterlist:
            index = masterlist.index(btype)
            sign_switch = 1.0
        elif btype in list(masterof.keys()):
            master, sign_switch = masterof[btype]
            index = masterlist.index(master)
        else:
            raise ValueError('No master found for bond %s of type %s' %(bond, btype))
        matrix[bond[0],index] +=  sign*sign_switch
        matrix[bond[1],index] += -sign*sign_switch
    #solve the set of equations q=M.t with q the full array of atomic charges
    #and t the array of bci masters
    bcis, res, rank, svals = np.linalg.lstsq(matrix, charges, rcond=1e-6)
    #print statistics if required
    if verbose:
        print('Fitting SQ to charges:')
        print('----------------------')
        print('    sing vals = ', svals)
        if min(svals)>0:
            print('    cond numb = ', max(svals)/min(svals))
        else:
            print('    cond numb = inf')
        print('')
        print('Resulting split charges:')
        print('------------------------')
        for btype, sq in zip(masterlist, bcis):
            print('  %10s    % .3f' %(btype, sq))
        print('')
        print('Statistics of the BCI charges:')
        print('------------------------------')
        apriori_values  = average(charges, ffatypes, 'sort')
        aposteriori_values = average(np.dot(matrix, bcis), ffatypes, 'sort')
        print(' %10s |  %6s +- %5s (%2s)  |  %6s +- %5s (%2s)  | %9s ' %('Atype', '<Qin>', 'std', 'N', '<Qbci>', 'std', 'N', 'RMSD'))
        print('    '+'-'*71)
        sums = np.array([0.0, 0.0, 0.0])
        for atype, qins in apriori_values.items():
            qouts = aposteriori_values[atype]
            print(' %10s |  % 6.3f +- %5.3f (%2i)  |  % 6.3f +- %5.3f (%2i)  | % 9.6f ' %(atype,
                qins.mean(), qins.std(), len(qins),
                qouts.mean(), qouts.std(), len(qouts),
                np.sqrt(((qouts-qins)**2).mean())
            ))
            sums += np.array([qins.sum(), qouts.sum(), ((qouts-qins)**2).sum()])
        print('    '+'-'*71)
        print(' %10s |  % 9.6f             |  % 9.6f             | % 9.6f ' %('ALL',
            sums[0],
            sums[1],
            np.sqrt(sums[2]/len(ffatypes))
        ))
        print('')
    #construct output dictionnary containing also bci's of slaves
    result = dict((btype, bci) for btype, bci in zip(masterlist, bcis))
    for slave, (master, sign) in masterof.items():
        result[slave] = result[master]*sign
    return result

def chebychev(m, x):
    if m==0:
        return 1
    elif m==1:
        return x
    else:
        return 2.0*x*chebychev(m-1,x)-chebychev(m-2,x)

def project_negative_freqs(hessian, masses, thresshold=0.0):
    N = len(masses)
    sqrt_mass_matrix = np.diag(np.sqrt((np.array([masses, masses, masses]).T).ravel()))
    isqrt_mass_matrix = np.linalg.inv(sqrt_mass_matrix)
    matrix = np.dot(isqrt_mass_matrix, np.dot(hessian.reshape([3*N,3*N]), isqrt_mass_matrix))
    #diagonalize
    if ((matrix-matrix.T)<1e-6*lightspeed/centimeter).all():
        evals, evecs = np.linalg.eigh(matrix)
    else:
        evals, evecs = np.linalg.eig(matrix)
    log.dump('20 lowest frequencies [1/cm] before projection:')
    log.dump(str(evals[:4]/(lightspeed/centimeter)))
    log.dump(str(evals[4:8]/(lightspeed/centimeter)))
    log.dump(str(evals[8:12]/(lightspeed/centimeter)))
    log.dump(str(evals[12:16]/(lightspeed/centimeter)))
    log.dump(str(evals[16:20]/(lightspeed/centimeter)))
    #set negative eigenvalues to zero
    evals[evals<thresshold] = 0.0
    projected_matrix = np.dot(evecs, np.dot(np.diag(evals), evecs.T))
    projected_hessian = np.dot(sqrt_mass_matrix, np.dot(projected_matrix, sqrt_mass_matrix))
    #dump freqs after projection as check
    evals, evecs = np.linalg.eigh(projected_matrix)
    log.dump('20 lowest frequencies [1/cm] after projection:')
    log.dump(str(evals[:4]/(lightspeed/centimeter)))
    log.dump(str(evals[4:8]/(lightspeed/centimeter)))
    log.dump(str(evals[8:12]/(lightspeed/centimeter)))
    log.dump(str(evals[12:16]/(lightspeed/centimeter)))
    log.dump(str(evals[16:20]/(lightspeed/centimeter)))
    return projected_hessian.reshape([N, 3, N, 3])
