#!/usr/bin/env python
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

from nose.tools import assert_raises
import numpy as np

from yaff import ForcePartValence, ForceField, PairPotEI, NeighborList,\
     ForcePartPair, Scalings

from quickff.refdata import *
from quickff.context import *
from quickff.tools import read_abinitio
from common import *

def test_covalent_warning():
    system = get_water()
    part = ForcePartValence(system)
    ff = ForceField(system, [part])
    # Check that the user can not pass a force field with covalent terms to
    # construct the reference data
    with assert_raises(UserWarning):
        refdata = ReferenceData(system.pos, ncff=ff)


def test_pbc_error():
    system = get_water()
    part = ForcePartValence(system)
    ff = ForceField(system, [part])
    # Check that the user can not pass pbc that do not match the ff pbc
    with assert_raises(AssertionError):
        refdata = ReferenceData(system.pos, ncff=ff, pbc=[1,1,0])


def test_hessian_projection_benzene():
    fn = context.get_fn('systems/benzene/gaussian.fchk' )
    numbers, coords, energy, grad, hess, masses, rvecs, pbc = read_abinitio(fn)
    refdata = ReferenceData(coords, grad=grad, hess=hess)
    u0, v0 = np.linalg.eigh(refdata.hess.reshape((coords.shape[0]*3,coords.shape[0]*3)))
    u1, v1 = np.linalg.eigh(refdata.phess.reshape((coords.shape[0]*3,coords.shape[0]*3)))
    u0 = np.sort(u0)
    u1 = np.sort(u1)
    # Benzene has 12 atoms, so 36 dofs...
    assert u0.shape[0]==36
    assert u1.shape[0]==36
    # ...but 3 are global translations and 3 are global vibrations.
    # We thus expect 6 small eigenvalues (in absolute value)
    assert np.sum(np.abs(u0/u0[7])<1e-2)==6
    # For the projected hessian, these 6 values have to be really small
    assert np.sum(np.abs(u1/u1[7])<1e-10)==6
    # All the other eigenvalues should not change alot
    assert np.all(np.abs(u0[6:]-u1[6:])<1e-6)


def test_hessian_projection_mil53():
    fn = context.get_fn('systems/mil53al_np/vasprun.xml')
    numbers, coords, energy, grad, hess, masses, rvecs, pbc = read_abinitio(fn)
    refdata = ReferenceData(coords, grad=grad, hess=hess, pbc=[1,1,1])
    u0, v0 = np.linalg.eigh(refdata.hess.reshape((coords.shape[0]*3,coords.shape[0]*3)))
    u1, v1 = np.linalg.eigh(refdata.phess.reshape((coords.shape[0]*3,coords.shape[0]*3)))
    u0 = np.sort(u0)
    u1 = np.sort(u1)
    # Mil53 unit cell has 76 atoms, so 228 dofs...
    assert u0.shape[0]==228
    assert u1.shape[0]==228
    # ...but 3 are global translations.
    # We thus expect 3 small eigenvalues (in absolute value)
    assert np.sum(np.abs(u0/u0[4])<1e-2)==3
    # For the projected hessian, these 6 values have to be really small
    assert np.sum(np.abs(u1/u1[4])<1e-10)==3
    # All the other eigenvalues should not change alot
    assert np.all(np.abs(u0[3:]-u1[3:])<1e-5)


def test_calc_energy_water():
    # Construct a force field for water with electrostatic contributions
    system = get_water()
    pair_pot = PairPotEI(system.charges, 0.0, 100.0)
    nlist = NeighborList(system)
    part = ForcePartPair(system, nlist, Scalings(system,1.0,1.0,1.0), pair_pot)
    nlist.update()
    ff = ForceField(system, [part])
    # Calculate gradient and hessian in equilibrium
    ndof = system.natom*3
    grad, hess = np.zeros((system.natom,3)), np.zeros((ndof, ndof))
    E0 = ff.compute(gpos=grad, hess=hess)
    refdata = ReferenceData(system.pos, energy=E0, grad=grad, hess=hess.reshape((system.natom,3,system.natom,3)), ncff=ff)
    # Check that different contributions are added correctly
    E1 = refdata.calc_energy(system.pos, ai=True, ff=False)
    E2 = refdata.calc_energy(system.pos, ai=False, ff=True)
    E3 = refdata.calc_energy(system.pos)
    assert E1+E2==E3
    assert E1==E0
    assert E2==-E0
    # Now displace the atoms a bit
    system.pos[:] += np.random.normal(0.0,0.02,system.pos.shape)
    E4 = ff.compute()
    E5 = refdata.calc_energy(system.pos, ai=True, ff=False)
    E6 = refdata.calc_energy(system.pos, ai=False, ff=True)
    E7 = refdata.calc_energy(system.pos)
    assert np.abs(E4/E5-1.0)<0.1
    assert np.abs(E5+E6-E7)<1e-15
    assert E4==-E6
