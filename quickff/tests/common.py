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

from molmod.units import parse_unit, angstrom, kjmol
import numpy as np, os

from yaff import System

from quickff.context import context

__all__ = ['get_water', 'get_ethanol']

def get_system(molecule, atypes_level='high', ei_path=None):
    moldir = context.get_fn('systems/%s' %molecule)
    if ei_path is not None:
        system = System.from_files(
            [os.path.join(moldir, 'gaussian.fchk'), os.path.join(moldir, 'gaussian_wpart.h5')],
            ei_path=ei_path
        )
    else:
        system = System.from_files([os.path.join(moldir, 'gaussian.fchk')])
    system.guess_ffatypes(atypes_level)
    system.determine_ics_from_topology()
    return system

def get_water():
    coords = np.array([
        [ 0.000000000,  0.763382315, -0.468300621],
        [ 0.000000000, -0.000000000,  0.117075156],
        [-0.000000000, -0.763382315, -0.468300621],
    ])*angstrom
    numbers = np.array([1,8,1])
    bonds = np.array([[0,1],[1,2]])
    fcharges = np.array([0.5, -1, 0.5]) #'formal' charges
    fradii = 0.5*np.array([2.571,3.118,2.571])/2.0**(1.0/6.0)*angstrom #0.5 rescaled UFF-vdW minima
    system = System(numbers, coords, charges=fcharges, radii=fradii)
    return system

def get_ethanol():
    coords = np.array([
        [ 1.207471, -0.223281,  0.000000],
        [-0.075681,  0.553519, -0.000000],
        [ 1.250572, -0.854108,  0.888799],
        [ 1.250572, -0.854109, -0.888799],
        [ 2.055418,  0.462705, -0.000000],
        [-0.147118,  1.198973, -0.883292],
        [-1.097876, -0.419344,  0.000000],
        [-0.147118,  1.198973,  0.883292],
        [-1.917880,  0.080549,  0.000000],
    ])*angstrom
    numbers = np.array([6, 6, 1, 1, 1, 1, 8, 1, 1])
    epsilons = np.array([0.439,0.439,0.184,0.184,0.184,0.184,0.251,0.184,0.184])*kjmol  #UFF
    sigmas = np.array([3.431,3.431,2.571,2.571,2.571,2.571,3.118,2.571,2.571])*angstrom #UFF
    fcharges = np.array([-0.3, 0.3, 0.1, 0.1, 0.1, 0.1, -1.0, 0.1, 0.5]) #'formal' charges
    fradii = 0.5*sigmas/2.0**(1.0/6.0) #0.5 rescaled UFF-vdW minima
    return coords, numbers, fcharges, fradii, epsilons, sigmas
