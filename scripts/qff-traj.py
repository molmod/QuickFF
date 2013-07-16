#! /usr/bin/env python

from optparse import OptionParser

from quickff.system import System
from quickff.model import Model
from quickff.program import Program
from quickff.fftable import FFTable

def parser():
    usage = "%prog [options] icname fns"
    description = 'This script is part of QuickFF. It will generate the'+\
                  'perturbation trajectory for the ic defined in icname'+\
                  'and plot the energy contributions. The system is defined'+\
                  'in the files fns'
    parser = OptionParser(usage=usage, description=description)
    parser.add_option(
        '--eirule', default=-1, type=int,
        help='Defines the exclusion rule for the electrostatic interactions. '  +\
             'If set to -1, no electrostatic interactions are taken into '      +\
             'account during force field fitting. If the rule is set to '       +\
             '0,1,2,3: pairs seperated by less then or equal to 0,1,2,3 bonds ' +\
             'respectively are excluded from ei interactions. [default=-1]'
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
    options, args = parser.parse_args()
    icname = args[0]
    fns = args[1:]
    return icname, fns, options

def main():
    icname, fns, options = parser()
    system = System.from_files(fns)
    if options.atypes_level is not None:
        system.guess_ffatypes(options.atypes_level)
    system.determine_ics_from_topology()
    model = Model.from_system(system, eirule=options.eirule)
    program = Program(system, model)
    program.plot_pt(icname)

if __name__=='__main__':
    main()
