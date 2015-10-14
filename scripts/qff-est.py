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

from yaff import System, ForceField
from yaff import log

from quickff.refdata import ReferenceData
from quickff.program import Program
from quickff.tools import guess_ffatypes, read_abinitio
from quickff.paracontext import *

def parser():
    usage = "%prog [options] fn"
    description = 'This script is part of QuickFF. It calculates the harmonic'+\
                  'ff parameters. The file fn should contain ab initio '+\
                  'reference data. Currently Gaussian .fchk and VASP .xml '+\
                  'files are supported.'
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
        '--ic-ids', default='all',
        help="A comma-separated list of identifiers specifying which icnames "  +\
             "should be included in the Valence Part. Each identifier can "     +\
             "be a specific IC name such as 'bond/C3_cc.H1_c' or can be one "   +\
             "of the following strings: 'bonds', 'angles', 'diheds', "          +\
             "'opdists' or 'all', in which case all bonds, angles, ... "        +\
             "will be included. [default=%default]"
    )
    parser.add_option(
        '--refine-rvs', default=False, action='store_true',
        help = 'Refine the rest values after the second step of force-field ' +\
               'development' +\
               '[default=False]'
    )
    parser.add_option(
        '--suffix', default='',
        help = "Suffix that will be added to all output files. [default='']"
    )
    parser.add_option(
        '--fn-traj', default=None,
        help='Load/Save perturbation trajectories to/from cPickled file.'
    )
    parser.add_option(
        '--scoop', default=False, action='store_true',
        help = 'Use scoop distribute tasks over workers. Note that the main ' +\
               'program has to be started with -m scoop for this to work, e.g. ' +\
               'python -m scoop -n4 ~/bin/qff-est.py --scoop ...' +\
               '[default=False]'
    )
    parser.add_option(
        '--linear', default=False, action='store_true',
        help = 'Use a linearized model to avoid construction of perturbation '+\
               'trajectories ' +\
               '[default=False]'
    )
    parser.add_option(
        '--strain-taylor', default=False, action='store_true',
        help = 'Use a Taylor expansion to approximate the strain function used '+\
               'to generate trajectories ' +\
               '[default=False]'
    )
    parser.add_option(
        '--yaff-output', default=False, action='store_true',
        help = 'Provide Yaff screen output. ' +\
               '[default=False]'
    )
    options, fn = parser.parse_args()
    options.ic_ids = options.ic_ids.split(',')
    return fn[0], options

def main(fn, options):
    if not options.yaff_output: log.set_level(log.silent)
    if options.linear and options.refine_rvs:
        raise UserWarning, "The options --linear and --refine-rvs are incompatible, choose one of both"
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
    if system.ffatypes is None: raise UserWarning, "No force-field atom types defined!"
    # Read non-covalent force field
    if options.ncff_fn is not None:
        ff = ForceField.generate(system, options.ncff_fn)
    else:
        ff = None
    # Construct reference data
    refdata = ReferenceData(coords, energy, grad, hess, ncff=ff, pbc=pbc)
    # Construct QuickFF program
    program = Program(system, refdata, fn_traj=options.fn_traj, skip_ics=options.ic_ids,
        refineq=options.refine_rvs,strain_taylor=options.strain_taylor)
    # Run the program
    qff = program.run(linear_model=options.linear)
    #Make output
    qff.dump_ffit2('pars_ffit2%s.txt' % options.suffix, mode='w')
    qff.dump_yaff('pars_yaff%s.txt' % options.suffix, mode='w')
    if options.ncff_fn is not None:
        with open('pars_yaff%s.txt' % options.suffix, mode='a') as fout:
            with open(options.ncff_fn, mode='r') as fin:
                for line in fin: fout.write(line)

#Use scoop if requested. This has to be outside of __main__ to set the
#context for all workers
#Parse args
fns, options = parser()
if options.scoop:
    paracontext.use_scoop()

if __name__=='__main__':
    main(fns, options)
