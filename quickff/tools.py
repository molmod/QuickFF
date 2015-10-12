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
    'guess_ffatypes', 'get_vterm', 'read_abinitio'
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


def dump_charges_yaff(system, fn, mode='w', scales=[1.0,1.0,1.0]):
    '''
        Write or append charges to a file in Yaff format.
        Assumes that charges and radii have been averaged over ffatypes.
    '''
    f = open(fn, mode)
    print >> f, ''
    print >> f, '# Fixed charges'
    print >> f, '# ============='
    print >> f, ''
    print >> f, "# Mathematical form: q_A = q_0A + sum'_B p_BA"
    print >> f, '# where q0_A is the reference charge of atom A. It is mostly zero, sometimes a'
    print >> f, '# non-zero integer. The total charge of a fragment is the sum of all reference'
    print >> f, '# charges. The parameter p_BA is the charge transfered from B to A. Such charge'
    print >> f, '# transfers are only carried out over bonds in the FF topology.'
    print >> f, '# The charge on an atom is modeled as a Gaussian distribution. The spread on the'
    print >> f, '# Gaussian is called the radius R. When the radius is set to zero, point charges'
    print >> f, '# will be used instead of smeared charges.'
    print >> f, ''
    print >> f, 'FIXQ:UNIT Q0 e'
    print >> f, 'FIXQ:UNIT P e'
    print >> f, 'FIXQ:UNIT R angstrom'
    print >> f, 'FIXQ:SCALE 1 %3.1f' % scales[0]
    print >> f, 'FIXQ:SCALE 2 %3.1f' % scales[1]
    print >> f, 'FIXQ:SCALE 3 %3.1f' % scales[2]
    print >> f, 'FIXQ:DIELECTRIC 1.0'
    print >> f, ''
    print >> f, '# Atomic parameters'
    print >> f, '# ----------------------------------------------------'
    print >> f, '# KEY        label  Q_0A              R_A'
    print >> f, '# ----------------------------------------------------'
    added = []
    for atype, q, radius in zip(system.ffatype_ids, system.charges, system.radii):
        if system.ffatypes[atype] not in added:
            print >> f, 'FIXQ:ATOM %8s % 13.10f  %12.10f' % (system.ffatypes[atype], q, radius/angstrom)
            added.append(system.ffatypes[atype])
    f.close()

def dump_vdw_yaff(system, fn, vdw, mode='w'):
    '''
        Write or append van der Waals parameters to a file in Yaff format.
    '''
    raise NotImplementedError
    if vdw.pot.kind.startswith('MM3Buckingham'):
        pot_id = 'MM3'
    elif vdw.pot.kind.startswith('LennartJones'):
        pot_id = 'LJ'
    else:
        raise ValueError('VDWPart in model has unsupported pot_kind: %s' %vdw.kind)
    f = open(fn, mode)
    print >> f, '# van der Waals'
    print >> f, '#=============='
    print >> f, '# The following mathemetical form is supported:'
    print >> f, '#  - MM3:   EPSILON*(1.84e5*exp(-12*r/SIGMA)-2.25*(SIGMA/r)^6)'
    print >> f, '#  - LJ:    4.0*EPSILON*((SIGMA/r)^12 - (SIGMA/r)^6)'
    print >> f, '#  - PAULI: A*exp(-B*r)'
    print >> f, ''
    print >> f, '%s:UNIT SIGMA angstrom' %(pot_id)
    print >> f, '%s:UNIT EPSILON kjmol' %(pot_id)
    print >> f, '%s:SCALE 1 %4.2f' % (pot_id, vdw.scales[0])
    print >> f, '%s:SCALE 2 %4.2f' % (pot_id, vdw.scales[1])
    print >> f, '%s:SCALE 3 %4.2f' % (pot_id, vdw.scales[2])
    print >> f, ''
    print >> f, '# -------------------------------------------'
    print >> f, '# KEY    ffatype  SIGMA         EPSILON'
    print >> f, '# -------------------------------------------'
    added = []
    for atype, sigma, epsilon in zip(self.ffatypes, self.sigmas, self.epsilons):
        if atype not in added:
            print >> f, '%s:PARS %8s % .10f % .10f' %(pot_id, atype, sigma/angstrom, epsilon/kjmol)
            added.append(atype)
    f.close()

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
