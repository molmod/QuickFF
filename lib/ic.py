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

import numpy as np

__all__ = ['IC']

class IC(object):
    '''
        A class functioning as a wrapper around the MolMod ic functions.
        This class should be able to calculate the value of an internal
        coordinate and its first and second order derivatives toward the
        cartesian coordinates.
    '''
    def __init__(self, name, indexes, icf, qunit='au', kunit='kjmol/au**2'):
        '''
            **Arguments**

            name
                a string defining a unique name for the ic

            indexes
                a list of indices refering to the relevant atoms of the system

            icf
                a method from the `molmod.ic <http://molmod.github.io/molmod/reference/algo.html#module-molmod.ic>`_
                module to calculate the value, gradient and hessian of the ic

            **Optional Arguments**

            qunit
                a string describing the conversion of the unit of the ic value.
                More info regarding possible strings can be found in the
                `MolMod documentation <http://molmod.github.io/molmod/reference/const.html#module-molmod.units>`_.

            kunit
                a string describing the conversion of the unit of the force
                constant. More info regarding possible strings can be found in
                the `MolMod documentation <http://molmod.github.io/molmod/reference/const.html#module-molmod.units>`_.
        '''
        self.name = name
        self.indexes = indexes
        self.icf = icf
        self.qunit = qunit
        self.kunit = kunit

    def value(self, coords):
        '''
            Calculate the value of the internal coordinate

            **Arguments**

            coords
                a (N,3) numpy array containing the cartesian coordinates
                of all atoms.
        '''
        _rs = coords[self.indexes]
        return self.icf(_rs, deriv=0)[0]

    def grad(self, coords):
        '''
            Calculate the first order derivative of the internal coordinate

            **Arguments**

            coords
                a (N,3) numpy array containing the cartesian coordinates
                of all atoms.
        '''
        _rs = coords[self.indexes]
        _q_deriv = self.icf(_rs, deriv=1)[1]
        grad = np.zeros(coords.shape, float)
        grad[self.indexes] = _q_deriv
        grad = grad.reshape([3*len(coords)])
        return grad

    def hess(self, coords):
        '''
            Calculate the second order derivative of the internal coordinate

            **Arguments**

            coords
                a (N,3) numpy array containing the cartesian coordinates
                of all atoms.
        '''
        Natoms = len(coords)
        _rs = coords[self.indexes]
        _q_hess = self.icf(_rs, deriv=2)[2]
        hess = np.zeros([Natoms, 3, Natoms, 3], float)
        for i_1, index1 in enumerate(self.indexes):
            for i_2, index2 in enumerate(self.indexes):
                hess[index1, :, index2, :] = _q_hess[i_1, :, i_2, :]
        hess = hess.reshape([3*Natoms, 3*Natoms])
        return hess
