#! /usr/bin/env python

import numpy as np
from optparse import OptionParser
from molmod.units import *

from quickff.system import System
from quickff.perturbation import estimate
from quickff.fftable import FFTable
from quickff.ffit import FFitProgram

def parser():
    usage = "%prog [options] icnames system"
    description = """This script is part of QuickFF. It calculates the harmonic ff parameters for the internal coordinates defined in <icnames> in the system defined in <system>. <icnames> is a list of strings (format = 'kind/attypes') defining the intenal coordinates icname, e.g. [bond/O.H, bend/H.O.H] in which O and H are atom types. The user can also write 'bond/all', 'bend/all' and/or 'dihedral/all' for defining all bonds, bends and dihedrals or write 'all' for all ics. <system> is a gaussian fchk file or molmod chk file."""
    parser = OptionParser(usage=usage, description=description)
    parser.add_option(
        '--psf', default=None,
        help='A MolMod PSF file from which the topology and atom types are taken. ' \
            +'If this argument is not present, it is assumed that all necessary info is already in the system file.'
    )
    parser.add_option(
        '--coupling', default=None,
        help='Specify whether ics of the same type should be coupled or not. Possible couplings are: symmetric. [default=%default]'
    )
    parser.add_option(
        '--free-depth', default=0, type=int,
        help='Defines the depth of the layer of atoms that are free. All other atoms are constrained to their original position with a spring. [default=%default]'
    )
    parser.add_option(
        '--spring', default=10.0*kjmol/angstrom**2, type=float,
        help='Defines the strength of the constraining spring [in kjmol/A^2]. [default=%default]'
    )
    options, args = parser.parse_args()
    options.spring = options.spring*kjmol/angstrom**2
    icnames = args[:-1]
    fn_chk = args[-1]
    return icnames, fn_chk, options

def main():
    icnames, fn_chk, options = parser()
    system = System('system', fn_chk)
    if options.psf is not None:
        system.topology_from_psf(options.psf)
    system.find_ic_patterns(icnames)
    system.dump_sample('system.chk')
    fftab_init = estimate(system, coupling=options.coupling, free_depth=options.free_depth, spring=options.spring)
    fftab_init.print_screen()
    fftab_init.dump_pars_ffit2('pars_init.txt')
    program = FFitProgram(system)
    fftab_fine = program.run()
    fftab_fine.print_screen()

if __name__=='__main__':
    main()
