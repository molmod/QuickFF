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

from argparse import ArgumentParser
import sys, os

description  = '''\
This script will apply QuickFF to derive a covalent force field for the given
system from the ab initio input given in the input files.'''
#The arguments fn1, fn2, ... represent all input files that specify the system and the ab initio reference data. 

def parse():
    parser = ArgumentParser(description=description)
    parser.add_argument(
        '--version', default=False, action='store_true',
        help='Print QuickFF version number and exit.'
    )
    parser.add_argument(
        '-s', '--silent', default=False, action='store_true',
        help='Swith of all logging completely, overwrites all other verbosity '
             'options.'
    )
    parser.add_argument(
        '-v', '--verbose', default=False, action='store_true',
        help='Increases verbosity, is overwriten if SILENT or VERY_VERBOSE is '
             'switched on.'
    )
    parser.add_argument(
        '-V', '--very-verbose', default=False, action='store_true',
        help='Increases verbosity to highest level, is overwriten if SILENT is '
             'switched on.'
    )
    parser.add_argument(
        '-l', '--logfile', default=None,
        help='Redirect logger output to a file with name LOGFILE.'
    )
    parser.add_argument(
        '--scoop', default=False, action='store_true',
        help='Flag to enable parallelisation using SCOOP. With SCOOP, the '
             'command to run QuickFF is slightly different, the absolute path '
             'to quickff.py should be used. For example, to run on 4 cores: '
             'python -m scoop -n4 /path/to/%(prog)s --scoop [options] fns'
    )
    #General settings options
    settings = parser.add_argument_group(title='General QuickFF specifications')
    settings.add_argument(
        '-m', '--program-mode', default='DeriveNonDiagFF',
        choices=[prog for prog in allowed_programs if not prog=='BaseProgram'],
        help='Specify the program mode which defines the set of instructions '
             'that will be executed.'
    )
    settings.add_argument(
        '--fn-traj', default=None,
        help='Read/write the perturbation trajectories from/to FN_TRAJ. If the '
             'given file exists, the trajectories are read from the file. '
             'Otherwise, the trajectories are written to the given file.'
    )
    settings.add_argument(
        '--only-traj', default=None,
        help='Construct the perturbation trajectory only for the terms with '+\
             'the given basenames. This options is only applied in the ' +\
             'MakeTrajectories program.'
    )
    settings.add_argument(
        '-e', '--ener-traj', default=False, action='store_true',
        help='Plot the various energy contributions along the perturbation '
             'trajectories to.'
    )
    settings.add_argument(
        '-x', '--xyz-traj', default=False, action='store_true',
        help='Write the perturbation trajectories in XYZ format. '
    )
    settings.add_argument(
        '--suffix', default='',
        help = 'Suffix that will be added to all output files.'
    )
    #Force field options
    ff = parser.add_argument_group(title='Options related to the definition and derivation of the force field')
    ff.add_argument(
        '--ei', default=None,
        help='A Yaff parameters file defining the electrostatic contribution '
             'of the force field.'
    )
    ff.add_argument(
        '--vdw', default=None,
        help='A Yaff parameters file defining the van der Waals contribution '
             'of the force field.'
    )
    ff.add_argument(
        '--covres', default=None,
        help='A Yaff parameters file defining a residual contribution to the '
             'covalent part of the force field.'
    )
    #System options
    system = parser.add_argument_group(title='Options related to the definition of the system')
    system.add_argument(
        '--ffatypes', default=None,
        choices=['None','low','medium','high','highest'],
        help='Assign atom types in the system through the automatic built-in '
             'detection (see documentation). By default (or if None is given), '
             'the atom types are assumed to be defined in the input files. '
             '[default=%(default)s]'
    )
    #Input files
    parser.add_argument(
        'fn', nargs='+',
        help='Input file name that specify the system and ab initio reference '
             'data. Multiple file names are allowed, but at least one should '
             'be given. Files later in the list overwrite information from '
             'earlier files. Allowed file formats are MolMod checkpoint files '
             '(file.chk), Gaussian formatted checkpoint files (file.fchk) '
             'and VASP xml files (file.xml).  '
    )
    args = parser.parse_args()
    if args.version:
        print version
        sys.exit()
    if args.ffatypes.lower()=='none':
        args.ffatypes = None
    return args

def main():
    args = parse()
    #define logger
    if args.silent:
        log.set_level('silent')
    else:
        if args.very_verbose:
            log.set_level('highest')
        elif args.verbose:
            log.set_level('high')
        if args.logfile is not None and isinstance(args.logfile, str):
            log.write_to_file(args.logfile)
    with log.section('QFF', 1, timer='Initializing'):
        log.dump('Initializing system')
        #read system and ab initio reference
        system = None
        energy = 0.0
        grad = None
        hess = None
        rvecs = None
        for fn in args.fn:
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
            if args.ffatypes in ['low', 'medium', 'high', 'highest']:
                guess_ffatypes(system, args.ffatypes)
            elif args.ffatypes is not None:
                raise NotImplementedError('Guessing atom types from %s not implemented' %args.ffatypes)
            else:
                raise AssertionError('No atom types defined')
        #construct ab initio reference
        ai = SecondOrderTaylor('ai', coords=system.pos.copy(), energy=energy, grad=grad, hess=hess, pbc=pbc)
        #detect a priori defined contributions to the force field
        refs = []
        if args.ei is not None:
            if rvecs is None:
                ff = ForceField.generate(system, args.ei, rcut=50*angstrom)
            else:
                ff = ForceField.generate(system, args.ei, rcut=20*angstrom, alpha_scale=3.2, gcut_scale=1.5, smooth_ei=True)
            refs.append(YaffForceField('EI', ff))
        if args.vdw is not None:
            ff = ForceField.generate(system, args.vdw, rcut=20*angstrom)
            refs.append(YaffForceField('vdW', ff))
        if args.covres is not None:
            ff = ForceField.generate(system, args.covres)
            refs.append(YaffForceField('Cov res', ff))
    #define quickff program
    assert args.program_mode in allowed_programs, \
        'Given program mode %s not allowed. Choose one of %s' %(
            args.program_mode,
            ', '.join([prog for prog in allowed_programs if not prog=='BaseProgram'])
        )
    mode = program_modes[args.program_mode]
    only_traj = 'PT_ALL'
    if args.only_traj is not None: only_traj = args.only_traj.split(',')
    program = mode(
        system, ai, ffrefs=refs,
        fn_traj=args.fn_traj, only_traj=only_traj,
        plot_traj=args.ener_traj, xyz_traj=args.xyz_traj,
        suffix=args.suffix
    )
    #run program
    program.run()

#Set the parallel context if scoop is enabled. This should be performed outside
#the __main__ block to ensure the set the context for all workers.
if '--scoop' in sys.argv[1:]:
    paracontext.use_scoop()

if __name__=='__main__':
    main()
