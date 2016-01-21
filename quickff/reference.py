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
from yaff.pes.ff import ForceField, ForcePartValence
from yaff.sampling.harmonic import estimate_cart_hessian
from quickff.tools import global_translation, global_rotation
import numpy as np

__all__ = ['SecondOrderTaylor', 'YaffForceField']


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
            U = np.linalg.svd(
                np.array([VTx, VTy, VTz, VRx, VRy, VRz]).transpose()
            )[0]
            nproj = 6
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
        grad0 = model.gradient(coords0)
        hess0 = model.hessian(coords0)
        name = model.name + ' (Harmonic)'
        return cls(name, coords=coords0, energy=energy0, grad=grad0, hess=hess0, pbc=model.pbc)
            
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
        return self.phess0


class YaffForceField(Reference):
    '''
        A model object based on a YAFF force field. Such a model can be used
        to account for non-covalent interactions or residual covalent
        interactions.
    '''
    def __init__(self, name, ff):
        self.ff = ff
        super(YaffForceField, self).__init__(name)
    
    def energy(self, coords):
        self.ff.update_pos(coords.copy())
        return self.ff.compute()
    
    def gradient(self, coords):
        gpos = np.zeros(coords.shape, float)
        self.ff.update_pos(coords.copy())
        energy = self.ff.compute(gpos=gpos)
        return gpos.reshape(coords.shap)
    
    def hessian(self, coords):
        self.ff.update_pos(coords.copy())
        hess = estimate_cart_hessian(self.ff)
        natoms = len(coords)
        return hess.reshape([natoms, 3, natoms, 3])
