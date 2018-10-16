#!/usr/bin/env python
# -*- coding: utf-8 -*-
# QuickFF is a code to quickly derive accurate force fields from ab initio input.
# Copyright (C) 2012 - 2018 Louis Vanduyfhuys <Louis.Vanduyfhuys@UGent.be>
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

from yaff import System
from molmod.io.fchk import FCHKFile

from quickff.reference import SecondOrderTaylor, get_ei_ff
from quickff.program import BaseProgram
from quickff.settings import Settings
from quickff.tools import get_ei_radii
from quickff.log import log

import h5py as h5

#define class for deriving the force field
class Program(BaseProgram):
    def run(self):
        with log.section('PROGRAM', 2):
            #deriving diagonal force field
            self.do_pt_generate()
            self.do_pt_estimate()
            self.average_pars()
            self.do_hc_estimatefc(['HC_FC_DIAG'])
            #adding and fitting cross terms
            self.do_cross_init()
            self.do_hc_estimatefc(['HC_FC_CROSS_ASS', 'HC_FC_CROSS_ASA', 'HC_FC_CROSS_DSS', 'HC_FC_CROSS_DSD', 'HC_FC_CROSS_DAA', 'HC_FC_CROSS_DAD'], logger_level=1)
            #write output
            self.make_output()


#initialize settings
settings = Settings(
    #general settings
    fn_yaff      = 'pars_cov.txt',
    fn_sys       = 'system.chk',
    plot_traj    = 'All', 
    xyz_traj     = True , 
    fn_traj      = 'trajectories.pp',
    log_level    = 'high',
    program_mode = 'DeriveFF',
    #algorithm settings
    do_hess_mass_weighting = True,
    do_hess_negfreq_proj   = True,
    do_cross_svd           = True,
    pert_traj_tol          = 1e-9,
    cross_svd_rcond        = 1e-9,
    #FF expression settings
    do_bonds      = True,
    do_bends      = True,
    do_dihedrals  = True,
    do_oops       = True,
    do_cross_ASS  = True,
    do_cross_ASA  = True,
    bond_term     = 'BONDHARM',
    bend_term     = 'BENDAHARM',
    do_squarebend = False,
    do_bendclin   = False,
    do_sqoopdist_to_oopdist = False,
)

#initialize System, ab initio input and FF references
with log.section('INIT', 1, timer='Initializing'):
    #load Gaussian Formatted Checkpoint file
    fchk = FCHKFile('gaussian.fchk')
    numbers = fchk.fields.get('Atomic numbers')
    energy = fchk.fields.get('Total Energy')
    coords = fchk.fields.get('Current cartesian coordinates').reshape([len(numbers), 3])
    grad = fchk.fields.get('Cartesian Gradient').reshape([len(numbers), 3])
    hess = fchk.get_hessian().reshape([len(numbers), 3, len(numbers), 3])

    #Construct Yaff System file
    system = System(numbers, coords)
    system.detect_bonds()
    system.set_standard_masses()
    
    #define atom types
    rules = [
        ('H', '1 & =1%8'), #hydrogen atom with one oxygen neighbor
        ('O', '8 & =2%1'), #oxygen atom with two hydrogen neighbors
    ]
    system.detect_ffatypes(rules)

    #Construct a QuickFF SecondOrderTaylor object containing the AI reference
    ai = SecondOrderTaylor('ai', coords=coords, energy=energy, grad=grad, hess=hess)

    #construct electrostatic force field from HE charges in gaussian_wpart.h5
    with h5.File('gaussian_mbis.h5') as f:
        charges = f['charges'][:]
    radii = get_ei_radii(numbers)
    scales = [1.0, 1.0, 1.0, 1.0]
    ff_ei = get_ei_ff('EI', system, charges, scales, radii=radii, average=True, pbc=[0,0,0])

#initialize and run program
program = Program(system, ai, settings, ffrefs=[ff_ei])
program.run()
