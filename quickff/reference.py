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

from __future__ import absolute_import

from yaff.pes.ff import ForceField, ForcePartValence, ForcePartPair
from yaff.pes.ext import PairPotEI
from yaff.pes.nlist import NeighborList
from yaff.pes.scaling import Scalings
from yaff.sampling.harmonic import estimate_cart_hessian

from quickff.tools import global_translation, global_rotation
from quickff.log import log

from molmod.units import angstrom

import numpy as np

__all__ = ['SecondOrderTaylor', 'YaffForceField', 'get_ei_ff']


class Reference(object):
    '''
        Abstract class for a model for reference data. A reference instance
        should be able to return energy, gradient and hessian for a given set
        of coordinates.
    '''

    def __init__(self, name):
        self.name = name

    def energy(self, coords):
        raise NotImplementedError

    def gradient(self, coords):
        raise NotImplementedError

    def hessian(self, coords):
        raise NotImplementedError


class SecondOrderTaylor(Reference):
    '''
        Second-order Taylor expansion model, can be used for the ab initio input.
    '''

    def __init__(self, name, coords=None, energy=0.0, grad=None, hess=None, pbc=[0,0,0]):
        log.dump('Initializing Second order taylor reference for %s' %name)
        self.coords0 = coords.copy()
        self.energy0 = energy
        self.grad0 = grad.copy()
        self.hess0 = hess.copy()
        assert np.all(np.array(pbc)==pbc[0]) and pbc[0] in [0,1], "PBC should be either all 0 or all 1"
        self.pbc = pbc
        self.phess0 = self._get_phess()
        super(SecondOrderTaylor, self).__init__(name)

    def update(self, coords=None, grad=None, hess=None, pbc=None):
        '''
            Method to update one or more of the attributes. The meaning of the
            optional arguments is the same as with the constructor __init__
        '''
        if pbc is not None:
            assert np.all(np.array(pbc)==pbc[0]) and pbc[0] in [0,1], "PBC should be either all 0 or all 1"
            self.pbc == pbc
        if coords is not None:
            if self.coords0 is not None:
                assert self.coords0.shape == coords.shape
            self.coords0 = coords
        if grad is not None:
            if self.grad0 is not None:
                assert self.grad0.shape == grad.shape
            self.grad0 = grad
        if hess is not None:
            if self.hess0 is not None:
                assert self.hess0.shape == hess.shape
            self.hess0 = hess
            self.phess0 = self._get_phess()

    def _get_phess(self):
        '''
            Constuct a hessian from which the translational and rotational
            degrees of freedom have been projected out.
        '''
        hess = self.hess0.copy().reshape([np.prod(self.coords0.shape), np.prod(self.coords0.shape)])
        VTx, VTy, VTz = global_translation(self.coords0)
        VRx, VRy, VRz = global_rotation(self.coords0)
        if np.all(np.array(self.pbc)==0):
            U, S, Vt = np.linalg.svd(
                np.array([VTx, VTy, VTz, VRx, VRy, VRz]).transpose()
            )
            nproj = len([s for s in S if s>1e-6])
            if nproj==5:
                log.dump('Only 5 out of the 6 trans-rot vectors were linearly independent. If the molecule is not linear, something went wrong!')
            elif not nproj==6:
                raise RuntimeError('Only %i of the 6 trans-rot vectors were linearly independent. Something went wrong!')
        elif np.all(np.array(self.pbc)==1):
            U = np.linalg.svd(
                np.array([VTx, VTy, VTz]).transpose()
            )[0]
            nproj = 3
        proj = np.dot(U[:, nproj:], U[:, nproj:].T)
        PHP = np.dot(proj, np.dot(hess, proj))
        return PHP.reshape([self.coords0.shape[0], 3, self.coords0.shape[0], 3])

    @classmethod
    def from_other_model(cls, model, coords0):
        energy0 = model.energy(coords0)
        grad0 = model.gradient(coords0).copy()
        hess0 = model.hessian(coords0).copy()
        name = model.name + ' (Harmonic)'
        if 'rvecs' in model.ff.system.__dict__:
            pbc = [1,1,1,]
        else:
            pbc = [0,0,0]
        return cls(name, coords=coords0, energy=energy0, grad=grad0, hess=hess0, pbc=pbc)

    def energy(self, coords):
        '''
            Compute the energy for the given positions
        '''
        assert np.all(coords.shape==self.coords0.shape)
        ndof = np.prod(self.coords0.shape)
        dx = (coords - self.coords0).reshape([ndof])
        energy = self.energy0 + np.dot(self.grad0.reshape([ndof]), dx)
        energy += 0.5*np.dot(dx, np.dot(self.phess0.reshape([ndof,ndof]), dx))
        return energy

    def gradient(self, coords):
        '''
            Compute the gradient for the given positions
        '''
        assert np.all(coords.shape==self.coords0.shape)
        ndof = np.prod(self.coords0.shape)
        dx = (coords - self.coords0).reshape([ndof])
        grad = self.grad0.reshape([ndof]) + np.dot(self.phess0.reshape([ndof,ndof]), dx)
        return grad.reshape(self.coords0.shape)

    def hessian(self, coords):
        '''
            Compute the hessian for the given positions
        '''
        assert np.all(coords.shape==self.coords0.shape)
        return self.phess0.copy()


class YaffForceField(Reference):
    '''
        A model object based on a YAFF force field. Such a model can be used
        to account for non-covalent interactions or residual covalent
        interactions.
    '''
    def __init__(self, name, ff):
        log.dump('Initializing Yaff force field reference for %s' %name)
        self.ff = ff
        Reference.__init__(self, name)

    def energy(self, coords):
        self.ff.update_pos(coords.copy())
        return self.ff.compute()

    def gradient(self, coords):
        gpos = np.zeros(coords.shape, float)
        self.ff.update_pos(coords.copy())
        energy = self.ff.compute(gpos=gpos)
        return gpos.reshape(coords.shape)

    def hessian(self, coords):
        self.ff.update_pos(coords.copy())
        hess = estimate_cart_hessian(self.ff)
        natoms = len(coords)
        return hess.reshape([natoms, 3, natoms, 3])


def get_ei_ff(name, system, charges, scales, radii=None, average=True, pbc=[0,0,0]):
    '''
        A routine to construct a Yaff force field for the electrostatics

        **Arguments**

        name
            A string for the name of the force field. This name will show in
            possible plots visualizing contributions along perturbation
            trajectories

        system
            A Yaff System instance representing the system

        charges
            A numpy array containing the charge of each separate atom

        scales
            A list of 4 floats, representing scale1, scale2, scale3, scale4,
            i.e. the electrostatic scaling factors

        **Optional Arguments**

        radii
            A numpy array containing the gaussian radii of each separate atom.
            If this argument is omitted, point charges are used.

        average
            If set to True, the charges and radii will first be averaged over
            atom types. This is True by default.
    '''
    if not (np.array(pbc)==0).all():
        raise NotImplementedError('Periodic system not implemented in get_ei_ff')
    if average:
        qs = {}
        rs = {}
        ffatypes = [system.ffatypes[i] for i in system.ffatype_ids]
        for i, atype in enumerate(ffatypes):
            if atype in qs: qs[atype].append(charges[i])
            else: qs[atype] = [charges[i]]
            if radii is not None:
                if atype in rs: rs[atype].append(radii[i])
                else: rs[atype] = [radii[i]]
        for i, atype in enumerate(ffatypes):
            charges[i] = np.array(qs[atype]).mean()
            if radii is not None:
                radii[i] = np.array(rs[atype]).mean()
    if radii is None: radii = np.zeros(len(system.ffatype_ids), float)
    if system.charges is None: system.charges = charges.copy()
    if system.radii is None: system.radii = radii.copy()
    pair_pot = PairPotEI(charges.astype(np.float), 0.0, 50*angstrom, None, 1.0, radii.astype(np.float))
    nlist = NeighborList(system, 0)
    scalings = Scalings(system, scale1=scales[0], scale2=scales[1], scale3=scales[2], scale4=scales[3])
    part = ForcePartPair(system, nlist, scalings, pair_pot)
    ff = ForceField(system, [part], nlist=nlist)
    return YaffForceField(name, ff)
