#!/usr/bin/env python
# -*- coding: utf-8 -*-
#QuickFF is a code to quickly derive accurate force fields from ab initio input.
#Copyright (C) 2012 - 2014 Louis Vanduyfhuys <Louis.Vanduyfhuys@UGent.be>
#Steven Vandenbrande <Steven.Vandenbrande@UGent.be>, 
#Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center for Molecular Modeling
#(CMM), Ghent University, Ghent, Belgium; all rights reserved unless otherwise
#stated.
#
#This file is part of QuickFF.
#
#QuickFF is free software; you can redistribute it and/or
#modify it under the terms of the GNU General Public License
#as published by the Free Software Foundation; either version 3
#of the License, or (at your option) any later version.
#
#QuickFF is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, see <http://www.gnu.org/licenses/>
#
#--

from optparse import OptionParser

from quickff.system import System
from quickff.model import Model
from quickff.program import Program
from quickff.fftable import FFTable

def parser():
    usage = "%prog [options] icname fns"
    description = 'This script is part of QuickFF. It will generate the '+\
                  'perturbation trajectory for the ic defined in icname '+\
                  'and plot the energy contributions. The system is defined '+\
                  'in the files fns.'
    parser = OptionParser(usage=usage, description=description)
    parser.add_option(
        '--ei-model', default='HarmPoint',
        help='Defines the potential used for the electrostatic interactions. '  +\
             'Can be CoulPoint, CoulGauss, HarmPoint, HarmGauss or Zero. If '   +\
             'CoulPoint/CoulGauss is chosen, the exact Coulombic potential '    +\
             'between point/gaussian charges will be used to evaluate EI '      +\
             'interactions. If HarmPoint/HarmGauss is chosen, a second order '  +\
             'Taylor expansion of the Coulomb potential is used. Harmonic is a '+\
             'lot faster and should already give accurate results. '            +\
             '[default=%default]'
    )
    parser.add_option(
        '--ei-scales', default='0.0,0.0,1.0', type=str,
        help='Defines the scaling rule for the electrostatic interactions. '    +\
             'Three comma-separated floats are required. The first one sets '   +\
             'the scale for atoms separated by 1 bond, the second for atoms '   +\
             'separated by 2 bonds etc ... By default, all interactions are '   +\
             'left unscaled [default=%default]'
    )
    parser.add_option(
        '--ei-path', default=None,
        help='Defines the path in the HDF5 file which contains a dataset '      +\
             '`EI_PATH/charges` from which the atomic charges will be '         +\
             'extracted.'
    )
    parser.add_option(
        '--vdw-model', default='Zero',
        help='Defines the potential used for van der Waals interactions. '      +\
             'Possible choices are: LJ, Harmonic and Zero. LJ implies using '   +\
             'the Lennart-Jones potential. Harmonic implies approximating the ' +\
             'LJ potential by means of a second order Taylor expansion. Zero '  +\
             'implies ignoring electrostatic interactions. [default=%default]'
    )
    parser.add_option(
        '--vdw-scales', default='0.0,0.0,1.0', type=str,
        help='Defines the scaling rule for the van der Waals interactions. '    +\
             'Three comma-separated floats are required. The first one sets '   +\
             'the scale for atoms separated by 1 bond, the second for atoms '   +\
             'separated by 2 bonds etc ... [default=%default]'
    )
    parser.add_option(
        '--vdw-path', default=None,
        help='Defines the path in the HDF5 file which contains 2 dataset: '     +\
             '`VDW_PATH/epsilons` and `VDW_PATH/sigmas` from which the atomic ' +\
             'vdW parameters will be extracted.'
    )
    parser.add_option(
        '--vdw-from', default='none',
        help='Defines from which force field to extract vdW parameters. If '    +\
             'this value is anythin else then None, the values extracted from ' +\
             'a HDF5 file will be overwritten. Currently only UFF is supported '+\
             '[default=%default]'
    )
    parser.add_option(
        '--atypes-level', default=None,
        help='Assign atom types according to level ATYPES_LEVEL. LOW will '     +\
             'assign atom types based only on atom number. MEDIUM will assign ' +\
             'atom types based on atom number and the number of neighbors. '    +\
             'HIGH will assign atom types based on atom number, number of '     +\
             'neighbors and the nature of the neighbors. HIGHEST will assign '  +\
             'atom types based on atom index. NONE will not guess but use the ' +\
             'atom types defined in the input files fns. [default=%default]'
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
    options.ei_scales = [float(x) for x in options.ei_scales.split(',')]
    options.vdw_scales = [float(x) for x in options.vdw_scales.split(',')]
    icname = args[0]
    fns = args[1:]
    return icname, fns, options

def main():
    #Parse args
    icname, fns, options = parser()
    #Setup system, model and program
    system = System.from_files(
        fns, ei_path=options.ei_path, vdw_path=options.vdw_path
    )
    if options.atypes_level is not None:
        system.guess_ffatypes(options.atypes_level)
    system.determine_ics_from_topology()
    if options.vdw_model.lower() != 'zero':
        if options.vdw_from.lower() == 'uff':
            system.read_uff_vdw()
        elif options.vdw_from.lower() != 'none':
            raise ValueError('Unsupported value for vdw_from, recieved %s' %options.vdw_from)
    model = Model.from_system(
        system,
        ei_scales=options.ei_scales, ei_pot_kind=options.ei_model,
        vdw_scales=options.vdw_scales, vdw_pot_kind=options.vdw_model,
    )
    program = Program(system, model, fn_traj='trajectories.pps')
    program.plot_pt(icname, start=options.start, end=options.end, steps=options.steps)

if __name__=='__main__':
    main()
