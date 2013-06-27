#! /usr/bin/env python

import numpy as np
from optparse import OptionParser
from molmod.units import *

from quickff.system import System
from quickff.model import Model
from quickff.program import Program
from quickff.fftable import FFTable

def parser():
    usage = "%prog [options] fns"
    description = 'This script is part of QuickFF. It calculates the harmonic'+\
                  'ff parameters a system specified in the files fns.'
    parser = OptionParser(usage=usage, description=description)
    parser.add_option(
        '--eirule', default=-1, type=int,
        help='Defines the exclusion rule for the electrostatic interactions.' +\
             'If set to -1, no electrostatic interactions are taken into '    +\
             'account during force field fitting. If the rule is set to '     +\
             '0,1,2,3: pairs seperated by less then or equal to 0,1,2,3 bonds'+\
             'respectively are excluded from ei interactions. [default=-1]'
    )
    parser.add_option(
        '--atypes-level', default=None,
        help='Assign atom types according to level ATYPES_LEVEL. LOW will'    +\
             'assign atom types based only on atom number. MEDIUM will assign'+\
             'atom types based on atom number and the number of neighbors.'   +\
             'HIGH will assign atom types based on atom number, number of'    +\
             'neighbors and the nature of the neighbors. HIGHEST will assign' +\
             'atom types based on atom index. NONE will not guess but use the'+\
             'atom types defined in the input files fns. [default=%default]'
    )
    options, fns = parser.parse_args()
    return fns, options

def main():
    #Parse args
    fns, options = parser()
    #Setup system, model and program
    system = System.from_files(fns)
    if options.atypes_level is not None:
        system.guess_ffatypes(options.atypes_level)
    system.determine_ics_from_topology()
    model = Model.from_system(system, eirule=options.eirule)
    program = Program(system, model)
    #Run program
    ff = program.run()
    #Make output
    ff.dump_ffit2('pars_ffit2.txt', mode='w')
    ff.dump_yaff('pars_yaff.txt', mode='w')
    system.dump_charges_yaff('pars_yaff.txt', options.eirule, mode='a')
    system.dump('system.chk')

if __name__=='__main__':
    main()
