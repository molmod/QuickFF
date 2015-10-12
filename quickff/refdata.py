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

'''Representation of the ab initio reference data that serve as input.
'''


import numpy as np
from quickff.tools import global_translation, global_rotation

__all__ = ['ReferenceData']

class ReferenceData(object):
    '''
        A class to store the reference data from which the force field will be
        derived. This contains the difference of the ab initio data and the
        non-covalent interactions.
    '''
    def __init__(self, coords, energy=0.0, grad=None, hess=None, ncff=None, pbc=[0,0,0]):
        '''
        **Arguments**

        coords
            A (N,3) numpy array containing the cartesian coordinates of all
            atoms in equilibrium.

        **Optional Arguments**

        energy
            A float giving the energy in equilibrium

        grad
            A (N,3) numpy array containing the cartesian gradient of all
            atoms in equilibrium. Ideally, all elements of grad are zero.

        hess
            A (N,3,N,3) numpy array containing the cartesian Hessian in
            equilibrium.

        ncff
            A Yaff ForceField instance representing the non-covalent
            interactions.

        pbc
            Periodic boundaries along cell vectors. Should be either all zeros
            (non-periodic) or all ones (3D periodic)
        '''
        self.coords = coords.copy()
        self.energy = energy
        self.grad = grad
        self.hess = hess
        self.ncff = ncff
        self.pbc = pbc
        assert np.all(np.array(self.pbc)==self.pbc[0]), "PBC should be either all 0 or all 1"
        if hess is not None:
            self.phess = self._get_phess(self.pbc)
        else:
            self.phess = None
        if self.ncff is not None:
            # Check that there is no valence part present in the non-covalent ff
            if 'valence' in [part.name for part in self.ncff.parts]:
                raise UserWarning, "The ncff contains covalent terms!"
            # Compute force-field gradient and hessian
            self.ncff.update_pos(self.coords)
            self.ff_grad = np.zeros(self.ncff.system.pos.shape,float)
            self.ff_hessian = np.zeros((np.prod(self.ncff.system.pos.shape),np.prod(self.ncff.system.pos.shape)),float)
            self.ncff.compute(gpos=self.ff_grad, hess=self.ff_hessian)

    def update(self, coords=None, grad=None, hess=None):
        '''
            Method to update one or more of the attributes. The meaning of the
            optional arguments is the same as with the constructor
            :meth:`quickff.refdata.ReferenceData`
        '''
        if coords is not None:
            if self.coords is not None:
                assert self.coords.shape == coords.shape
            self.coords = coords
        if grad is not None:
            if self.grad is not None:
                assert self.grad.shape == grad.shape
            self.grad = grad
        if hess is not None:
            if self.hess is not None:
                assert self.hess.shape == hess.shape
            self.hess = hess
            self.phess = self._get_phess()

    def calc_energy(self, coords, harmonic=False, ai=True, ff=True):
        '''
            Compute the reference energy for the given positions

            **Optional Arguments**

            harmonic
                Use a harmonic approximation to compute the non-covalent force-
                field contribution.

            ai
                Compute the ab initio contribution

            ff
                Compute the non-covalent force-field contribution
        '''
        assert np.all(coords.shape==self.coords.shape)
        ndof = np.prod(self.coords.shape)
        dx = (coords - self.coords).reshape([ndof])
        energy, grad, hess = 0.0, np.zeros(ndof), np.zeros((ndof, ndof))
        # Ab initio part
        if ai:
            energy += self.energy
            grad[:] = self.grad.reshape([ndof])
            hess[:] = self.phess.reshape([ndof, ndof])
        # Subtract non-covalent ff contribution
        if ff and self.ncff is not None:
            if harmonic:
                grad[:] -= self.ff_grad.reshape([ndof])
                hess[:] -= self.ff_hess.reshape([ndof, ndof])
            else:
                self.ncff.update_pos(coords)
                energy -= self.ncff.compute()
        energy += np.dot(grad, dx)
        energy += 0.5*np.dot(dx, np.dot(hess, dx))
        return energy

    def _check(self, natoms=None):
        'Internal consistency check'
        assert self.coords is not None
        assert self.grad is not None
        assert self.hess is not None
        assert len(self.coords.shape) == 2
        if natoms is not None:
            assert self.coords.shape[0] == natoms
        assert self.coords.shape[1] == 3
        assert self.grad.shape == self.coords.shape
        assert len(self.hess.shape) == 4
        assert self.hess.shape[0] == self.coords.shape[0]
        assert self.hess.shape[1] == 3
        assert self.hess.shape[2] == self.coords.shape[0]
        assert self.hess.shape[3] == 3

    def _get_phess(self, pbc):
        '''
            Constuct a hessian from which the translational and rotational
            degrees of freedom have been projected out.
        '''
        hess = self.hess.copy().reshape([np.prod(self.coords.shape), np.prod(self.coords.shape)])
        VTx, VTy, VTz = global_translation(self.coords)
        VRx, VRy, VRz = global_rotation(self.coords)
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
        return PHP.reshape([self.coords.shape[0], 3, self.coords.shape[0], 3])
