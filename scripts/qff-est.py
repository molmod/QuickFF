#! /usr/bin/env python

import numpy as np, os, time
from optparse import OptionParser
from molmod.units import *

from quickff.system import System
from quickff.perturbation import *
from quickff.fftable import FFTable
from quickff.ffit import ZFitProgram, MFitProgram

def parser():
    usage = "%prog [options] icnames system"
    description = """This script is part of QuickFF. It calculates the harmonic ff parameters for the internal coordinates defined in <icnames> in the system defined in <system>. <icnames> is a list of strings (format = 'kind/attypes') defining the intenal coordinates icname, e.g. [bond/O.H, bend/H.O.H] in which O and H are atom types. The user can also write 'all' for all ics. <system> is a gaussian fchk file or molmod chk file."""
    parser = OptionParser(usage=usage, description=description)
    parser.add_option(
        '--psf', default=None,
        help='A MolMod PSF file from which the topology and atom types are taken. ' \
            +'If this argument is not present, the topology is taken from the geometry and the atom symbols are used as atom types.'
    )
    parser.add_option(
        '--cost-strain', default=1.0, type=float,
        help='The weight of the strain cost in the total cost function for calculating the perturbation on the geometry. The strain cost expresses for a given geometry the deviation of the internal coordinates from their equilibrium value. [default = %default]'
    )
    parser.add_option(
        '--cost-energy', default=0.0, type=float,
        help='The weight of the energy cost in the total cost function for calculating the perturbation on the geometry. The energy is calculated as the second order Taylor expression using the forces and hessian in equilibrium [default = %default]'
    )
    parser.add_option(
        '--relative-amplitude', default=0.05, type=float,
        help='Defines the relative amplitude of the change in the ic value relative to the ic equilibrium value. [default=%default]'
    )
    parser.add_option(
        '--ei-rule', default=-1, type=int,
        help='Defines the exclusion rule for the electrostatic interactions. If set to -1, no electrostatic interactions are taken into account during force field fitting. If the rule is set to 0,1,2,3: pairs seperated by less then or equal to 0,1,2,3 bonds respectively are excluded from ei interactions. [default=-1]'
    )
    parser.add_option(
        '--charges', default=None, 
        help='Overwrite the charges of the system. The charges are comma seperated and the order should be identical to the order of atoms in the system file. Alternatively, a hipart charge txt file can also be used to specify the charges. If no charges are defined but electrostatics are switched on (ei-rule>-1), the charges are taken from the system file if possible, if not, an error is raised. [default=%default]'
    )
    parser.add_option(
        '--atypes-level', default='None', 
        help='Overwrite the atom types according to level ATYPES_LEVEL. Low will choose atom types based only on atom number, medium will choose atom types based on local topology, high will choose atom types based on atom index and None will not guess atom types. [default=%default]'
    )
    parser.add_option(
        '--no-remove', default=True, dest='remove', action='store_false', 
        help='Do not remove intermediate files and directories.'
    )
    parser.add_option(
        '--only-system', default=False, dest='only_system', action='store_true', 
        help='Only construct a QuickFF chk and Yaff chk file containing the system.'
    )
    options, args = parser.parse_args()
    if options.atypes_level=='None': options.atypes_level=None
    icnames = args[:-1]
    fn_chk = args[-1]
    return icnames, fn_chk, options

def main():
    icnames, fn_chk, options = parser()
    system = System(fn_chk, fn_psf=options.psf, guess_atypes_level=options.atypes_level, charges=options.charges)
    if icnames==['all']:
        icnames = system.icnames
    else:
        for icname in icnames:
            if not icname in system.icnames:
                raise ValueError('Illegal icname in command arguments, recieved %s' %icname)
    system.define_models(eirule=options.ei_rule)
    pt = RelaxedGeometryPT(system, energy_penalty=options.cost_energy, strain_penalty=options.cost_strain, dq_rel=options.relative_amplitude)
    
    if options.only_system:
        system.dump_sample_qff('system_qff.chk')
        system.dump_sample_yaff('system_yaff.chk')
        return

    fftab_init = pt.estimate(icnames)
    fftab_init.print_screen()
    fftab_init.dump_pars_ffit2('int-pars.txt')
    if options.ei_rule>-1:
        zfit = ZFitProgram(system)
        charges = zfit.run()
        system.sample['ac'] = charges
    system.dump_sample_qff('int-system.chk')
    
    mfit = MFitProgram(system, icnames)
    fftab_fine = mfit.run()
    fftab_fine.print_screen()
    
    system.dump_sample_qff('system_qff.chk')
    system.dump_sample_yaff('system_yaff.chk')
    if options.remove:
        os.system('rm -r int-system.chk int-pars.txt out')
    print 'SYSTEM TIMER:', time.ctime()

if __name__=='__main__':
    main()
