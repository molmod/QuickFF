#! /usr/bin/env python

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
        '--ei-model', default='Harm',
        help='Defines the potential used for the electrostatic interactions. '  +\
             'Possible choices are: Coulomb, Harmonic and Zero. Coulomb '       +\
             'implies using the coulomb potential. Harmonic implies '           +\
             'approximating the coulomb potential by means of a second order '  +\
             'Taylor expansion. Zero implies ignoring electrostatic '           +\
             'interactions. [default=%default]'
    )
    parser.add_option(
        '--ei-scales', default='1.0,1.0,1.0', type=str,
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
        '--vdw-model', default='Harm',
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
             '`EI_PATH/epsilons` and `EI_PATH/sigmas` from which the atomic '   +\
             'vdW parameters will be extracted.'
    )
    parser.add_option(
        '--vdw-from', default=None,
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
        '--nb-exclusion', default=None,
        help='A Regular Expression to match atom types that are to be excluded '+\
             'in the non-bonded interactions [default=%default].'
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
        elif optins.vdw_from.lower is not None:
            raise ValueError('Unsupported value for vdw_from, recieved %s' %options.vdw_from)
    model = Model.from_system(
        system,
        ei_scales=options.ei_scales, ei_pot_kind=options.ei_model,
        vdw_scales=options.vdw_scales, vdw_pot_kind=options.vdw_model,
    )
    program = Program(system, model, fns_traj='trajectories.pps')
    program.plot_pt(icname, start=options.start, end=options.end, steps=options.steps)

if __name__=='__main__':
    main()
