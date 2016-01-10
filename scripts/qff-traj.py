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

from optparse import OptionParser

from molmod.units import *

from yaff import System, ForceField
from yaff import log

from quickff.refdata import ReferenceData
from quickff.program import Program
from quickff.tools import guess_ffatypes, read_abinitio

def parser():
    usage = "%prog [options] icname fn"
    description = 'This script is part of QuickFF. It will generate the '+\
                  'perturbation trajectory for the ic defined in icname '+\
                  'and plot the energy contributions. The system is defined '+\
                  'in the files fns.'
    parser = OptionParser(usage=usage, description=description)
    parser.add_option(
        '--chk-fn', default=None,
        help='Read Yaff system from this .chk file. If not given, the Yaff '+\
             'system is constructed from the ab initio data.'
    )
    parser.add_option(
        '--ncff-fn', default=None,
        help='Generate a non-covalent force field from this Yaff parameter file.'
    )
    parser.add_option(
        '--atypes-level', default=None,
        help='Assign atom types according to level ATYPES_LEVEL. Low will '     +\
             'assign atom types based only on atom number. Medium will assign ' +\
             'atom types based on atom number and the number of neighbors. '    +\
             'High will assign atom types based on atom number, number of '     +\
             'neighbors and the nature of the neighbors. Highest will assign '  +\
             'atom types based on atom index. None will not guess but use the ' +\
             'atom types defined in the input files fns. [default=%default]'
    )
    parser.add_option(
        '--yaff-output', default=False, action='store_true',
        help = 'Provide Yaff screen output. ' +\
               '[default=False]'
    )
    parser.add_option(
        '--start', default=None, type=float,
        help='Defines the smallest value (in atomic untis) of the perturbed '
             'ic to generate the perturbation trajectory.'
    )
    parser.add_option(
        '--end', default=None, type=float,
        help='Defines the highest value (in atomic untis) of the perturbed '
             'ic to generate the perturbation trajectory.'
    )
    parser.add_option(
        '--steps', default=51, type=int,
        help='Defines the number of steps to generate the perturbation '
             'trajectory. [default=%default]'
    )
    options, args = parser.parse_args()
    icname = args[0]
    fn = args[1]
    return icname, fn, options

def main():
    #Parse args
    icname, fn, options = parser()
    if not options.yaff_output: log.set_level(log.silent)
    numbers, coords, energy, grad, hess, masses, rvecs, pbc = read_abinitio(fn)
    # Setup system
    if options.chk_fn is not None:
        system = System.from_file(options.chk_fn)
        #TODO Check that there is some correspondence to ab initio data
    else:
        system = System(numbers, coords, masses=masses, rvecs=rvecs)
        system.detect_bonds()
    if options.atypes_level is not None:
        guess_ffatypes(system, options.atypes_level)
    if system.ffatypes is None: raise NotImplementedError
    # Read non-covalent force field
    if options.ncff_fn is not None:
        if pbc[0]==0: rcut = 50*angstrom
        else: rcut = 15*angstrom
        ff = ForceField.generate(system, options.ncff_fn, rcut=rcut)
    else:
        ff = None
    # Construct reference data
    refdata = ReferenceData(coords, energy, grad, hess, ncff=ff, pbc=pbc)
    # Construct QuickFF program
    program = Program(system, refdata, fn_traj='trajectories.pps')
    program.plot_pt(icname, start=options.start, end=options.end, steps=options.steps)

if __name__=='__main__':
    main()
