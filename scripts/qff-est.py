#! /usr/bin/env python

import numpy as np, os
from optparse import OptionParser
from molmod.units import *

from quickff.system import System
from quickff.perturbation import estimate
from quickff.fftable import FFTable
from quickff.ffit import ZFitProgram, MFitProgram

def parser():
    usage = "%prog [options] icnames system"
    description = """This script is part of QuickFF. It calculates the harmonic ff parameters for the internal coordinates defined in <icnames> in the system defined in <system>. <icnames> is a list of strings (format = 'kind/attypes') defining the intenal coordinates icname, e.g. [bond/O.H, bend/H.O.H] in which O and H are atom types. The user can also write 'bond/all', 'bend/all' and/or 'dihedral/all' for defining all bonds, bends and dihedrals or write 'all' for all ics. <system> is a gaussian fchk file or molmod chk file."""
    parser = OptionParser(usage=usage, description=description)
    parser.add_option(
        '--psf', default=None,
        help='A MolMod PSF file from which the topology and atom types are taken. ' \
            +'If this argument is not present, the topology is taken from the geometry and the atom symbols are used as atom types.'
    )
    parser.add_option(
        '--coupling', default=None,
        help='Specify whether ics of the same type should be coupled or not. Possible couplings are symmetric and None. [default=%default]'
    )
    parser.add_option(
        '--free-depth', default=0, type=int,
        help='Defines the depth of the layer of atoms that are free. All other atoms are constrained to their original position with a spring. [default=%default]'
    )
    parser.add_option(
        '--spring', default=10.0, type=float,
        help='Defines the strength of the constraining spring [in kjmol/A^2]. [default=%default]'
    )
    parser.add_option(
        '--ei-rule', default=-1, type=int,
        help='Defines the exclusion rule for the electrostatic interactions. If set to -1, no electrostatic interactions are taken into account during force field fitting. If the rule is set to 0,1,2,3: pairs seperated by less then or equal to 1,2,3 bonds respectively are excluded from ei interactions. [default=-1]'
    )
    parser.add_option(
        '--charges', default=None, 
        help='Specify the charges of the system, this is necessary if a ei-rule is specified different then -1. The charges are comma seperated and the order should be identical to the order in the sample file. If the charges are not specified but ei-rule is larger then -1, the charges are taken from the psf if present, or from the sample file if present. Alternatively, a hipart charge txt file can also be used to specify the charges. [default=%default]'
    )
    parser.add_option(
        '--no-remove', default=True, dest='remove', action='store_false', 
        help='Do not remove intermediate files and directories.'
    )
    parser.add_option(
        '--only-system', default=False, dest='only_system', action='store_true', 
        help='Only construct MolMod chk file containing the system sample from system file and possible psf file.'
    )
    options, args = parser.parse_args()
    options.spring = options.spring*kjmol/angstrom**2
    icnames = args[:-1]
    fn_chk = args[-1]
    return icnames, fn_chk, options

def main():
    icnames, fn_chk, options = parser()
    if options.ei_rule==-1:
        eikind = 'Zero'
    else:
        eikind = 'Harmonic'
    system = System('system', fn_chk, fn_psf=options.psf, eikind=eikind, eirule=options.ei_rule, charges=options.charges)
    system.find_ic_patterns(icnames)
    
    if options.only_system:
        system.dump_sample('system.chk')
    else:
        fftab_init = estimate(system, coupling=options.coupling, free_depth=options.free_depth, spring=options.spring)
        fftab_init.print_screen()
        fftab_init.dump_pars_ffit2('int-pars.txt')
        system.dump_sample('int-system.chk')
        if options.ei_rule>-1:
            zfit = ZFitProgram(system)
            charges = zfit.run()
            system.sample['ac'] = charges
            system.sample['charges'] = charges
            system.dump_sample('int-system.chk')
        mfit = MFitProgram(system)
        fftab_fine = mfit.run()
        fftab_fine.print_screen()
        system.dump_sample('system.chk')
        if options.remove:
            os.system('rm -r int-system.chk int-pars.txt out')

if __name__=='__main__':
    main()
