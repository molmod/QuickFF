#! /usr/bin/env python
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

from quickff.program import __all__ as allowed_programs, __dict__ as program_modes
from quickff.log import log, version
from quickff.tools import guess_ffatypes
from quickff.reference import SecondOrderTaylor, YaffForceField
from quickff.io import read_abinitio
from quickff.paracontext import paracontext

from yaff import System, ForceField, Cell

from molmod.units import angstrom
from molmod.io.chk import load_chk

from optparse import OptionParser, OptionGroup
import sys, os

def parse():
    usage = '%prog [options] fn1 [fn2 [...]]'
    description  = 'This script will apply QuickFF to derive a covalent force '
    description += 'field for the given system. The arguments fn1, fn2, ... '
    description += 'represent all input files that specify the system and the '
    description += 'ab initio reference data. Files later in the list '
    description += 'overwrite earlier files'
    parser = OptionParser(usage=usage, description=description)
    parser.add_option(
        '--scoop', default=False, action='store_true',
        help='Flag to enable parallelisation using SCOOP. With SCOOP, the '
             'command to run QuickFF is slightly different, the absolute path '
             'to quickff.py should be used. For example, to run on 4 cores:\n'
             'python -m scoop -n4 /path/to/%prog --scoop [options] fns\n'
             '[default=%default]'
    )
    parser.add_option(
        '-s', '--silent', default=False, action='store_true',
        help='Swith of all logging completely, overwrites all other verbosity '
             'options. [default=%default]'
    )
    parser.add_option(
        '-v', '--verbose', default=False, action='store_true',
        help='Increases verbosity, is overwriten if SILENT or VERY_VERBOSE is '
             'switched on. [default=%default]'
    )
    parser.add_option(
        '-V', '--very-verbose', default=False, action='store_true',
        help='Increases verbosity to highest level, is overwriten if SILENT is '
             'switched on. [default=%default]'
    )
    parser.add_option(
        '-l', '--logfile', default=None,
        help='Redirect logger output to a file with the given name. '
             '[default=%default]'
    )
    parser.add_option(
        '--version', default=False, action='store_true',
        help='Print QuickFF version number and exit.'
    )
    #General settings options
    settings = OptionGroup(parser, 'General', 'General specifications of the program')
    settings.add_option(
        '-m', '--program-mode', default='DeriveNonDiagFF',
        help='Specify the program mode which defines the set of instructions '
             'that will be executed. Allowed strings are the names of the '
             'program classes defined in quickff/program.py, which are: '
             ''+', '.join([prog for prog in allowed_programs if not prog=='BaseProgram'])+' (case sensitive). [default=%default]'
    )
    settings.add_option(
        '--fn-traj', default=None,
        help='Name for the file to read/write the perturbation trajectories '
             'from/to. If the given file exists, the trajectories are read '
             'from the file. Otherwise, the trajectories are written to the '
             'given file. [default=%default]'
    )
    settings.add_option(
        '--only-traj', default=None,
        help='Construct the perturbation trajectory only for the terms with '+\
             'the given basenames. This options is only applied in the ' +\
             'MakeTrajectories program.'
    )
    settings.add_option(
        '-e', '--ener-traj', default=False, action='store_true',
        help='Plot the various energy contributions along the perturbation '
             'trajectories to. [default=%default]'
    )
    settings.add_option(
        '-x', '--xyz-traj', default=False, action='store_true',
        help='Write the perturbation trajectories in XYZ format. '
             '[default=%default]'
    )
    settings.add_option(
        '--suffix', default='',
        help = "Suffix that will be added to all output files. [default='']"
    )
    parser.add_option_group(settings)
    #Force field options
    ff = OptionGroup(parser, 'Force Field', 'All options related to the definition and derivation of the force field.')
    ff.add_option(
        '--ei', default=None,
        help='A Yaff parameters file defining the electrostatic contribution '
             'of the force field. [default=%default]'
    )
    ff.add_option(
        '--vdw', default=None,
        help='A Yaff parameters file defining the van der Waals contribution '
             'of the force field. [default=%default]'
    )
    ff.add_option(
        '--covres', default=None,
        help='A Yaff parameters file defining a residual contribution to the '
             'covalent part of the force field. [default=%default]'
    )
    parser.add_option_group(ff)
    #System options
    system = OptionGroup(parser, 'System', 'All options related to the definition of the system.')
    system.add_option(
        '--ffatypes', default=None,
        help='Assign atom types in the system. If FFATYPES is a string equal to '
             'either low,medium, high or highest, atom types will be assigned '
             'through the automatic built-in detection (see documentation). If '
             'FFATYPES is None, the atom types are assumed to be defined in the '
             'input files. [default=%default]'
    )
    parser.add_option_group(system)
    options, args = parser.parse_args()
    if options.version:
        print version
        sys.exit()
    assert len(args)>0, 'No input files found.'
    return options, args

def main():
    options, fns = parse()
    #define logger
    if options.silent:
        log.set_level('silent')
    else:
        if options.very_verbose:
            log.set_level('highest')
        elif options.verbose:
            log.set_level('high')
        if options.logfile is not None and isinstance(options.logfile, str):
            log.write_to_file(options.logfile)
    with log.section('QFF', 1, timer='Initializing'):
        log.dump('Initializing system')
        #read system and ab initio reference
        system = None
        energy = 0.0
        grad = None
        hess = None
        rvecs = None
        for fn in fns:
            if fn.endswith('.fchk') or fn.endswith('.xml'):
                numbers, coords, energy, grad, hess, masses, rvecs, pbc = read_abinitio(fn)
                if system is None:
                    system = System(
                        numbers, coords, rvecs=rvecs, charges=None, radii=None,
                        masses=masses
                    )
                else:
                    system.pos = coords.copy()
                    system.cell = Cell(rvecs)
                    system.numbers = numbers.copy()
                    if masses is not None: system.masses = masses.copy()
                    system._init_derived()
            elif fn.endswith('.chk'):
                sample = load_chk(fn)
                if 'energy' in sample.keys():       energy = sample['energy']
                if 'grad' in sample.keys():         grad = sample['grad']
                elif 'gradient' in sample.keys():   grad = sample['gradient']
                if 'hess' in sample.keys():         hess = sample['hess']
                elif 'hessian' in sample.keys():    hess = sample['hessian']
                if system is None:
                    system = System.from_file(fn)
                else:
                    if 'pos' in sample.keys():          system.pos = sample['pos']
                    elif 'coords' in sample.keys():     system.pos = sample['coords']
                    if 'rvecs' in sample.keys():        system.cell = Cell(sample['rvecs'])
                    elif 'cell' in sample.keys():       system.cell = Cell(sample['cell'])
                    if 'bonds' in sample.keys():        system.bonds = sample['bonds']
                    if 'ffatypes' in sample.keys():     system.ffatypes = sample['ffatypes']
                    if 'ffatype_ids' in sample.keys():  system.ffatype_ids = sample['ffatype_ids']
                    system._init_derived()
            else:
                raise NotImplementedError('File format for %s not supported' %fn)
        assert system is not None, 'No system could be defined from input'
        assert grad is not None, 'No ab initio gradient found in input'
        assert hess is not None, 'No ab initio hessian found in input'
        #complete the system information
        if system.bonds is None: system.detect_bonds()
        if system.masses is None: system.set_standard_masses()
        if system.ffatypes is None:
            if options.ffatypes in ['low', 'medium', 'high', 'highest']:
                guess_ffatypes(system, options.ffatypes)
            elif options.ffatypes is not None:
                raise NotImplementedError('Guessing atom types from %s not implemented' %options.ffatypes)
            else:
                raise AssertionError('No atom types defined')
        #construct ab initio reference
        ai = SecondOrderTaylor('ai', coords=system.pos.copy(), energy=energy, grad=grad, hess=hess, pbc=pbc)
        #detect a priori defined contributions to the force field
        refs = []
        if options.ei is not None:
            if rvecs is None:
                ff = ForceField.generate(system, options.ei, rcut=50*angstrom)
            else:
                ff = ForceField.generate(system, options.ei, rcut=20*angstrom, alpha_scale=3.2, gcut_scale=1.5, smooth_ei=True)
            refs.append(YaffForceField('EI', ff))
        if options.vdw is not None:
            ff = ForceField.generate(system, options.vdw, rcut=20*angstrom)
            refs.append(YaffForceField('vdW', ff))
        if options.covres is not None:
            ff = ForceField.generate(system, options.covres)
            refs.append(YaffForceField('Cov res', ff))
    #define quickff program
    assert options.program_mode in allowed_programs, \
        'Given program mode %s not allowed. Choose one of %s' %(
            options.program_mode,
            ', '.join([prog for prog in allowed_programs if not prog=='BaseProgram'])
        )
    mode = program_modes[options.program_mode]
    only_traj = 'PT_ALL'
    if options.only_traj is not None: only_traj = options.only_traj.split(',')
    program = mode(
        system, ai, ffrefs=refs,
        fn_traj=options.fn_traj, only_traj=only_traj,
        plot_traj=options.ener_traj, xyz_traj=options.xyz_traj,
        suffix=options.suffix
    )
    #run program
    program.run()

#Set the parallel context if scoop is enabled. This should be performed outside
#the __main__ block to ensure the set the context for all workers.
if '--scoop' in sys.argv[1:]:
    paracontext.use_scoop()

if __name__=='__main__':
    main()
