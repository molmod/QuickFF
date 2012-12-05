#! /usr/bin/env python

import numpy as np, matplotlib.pyplot as pp
from optparse import OptionParser
from molmod.units import *

from quickff.system import System
from quickff.perturbation import calculate_perturbation, analyze_perturbation, plot_perturbation
from quickff.tools import fitpar
from quickff.evaluators import *

def parser():
    usage = "%prog [options] icname system"
    description = """This script is part of QuickFF. It calculates the perturbation trajectory for the internal coordinate defined in <icname> in the system defined in <system>. <icname> has the format 'kind/attypes', e.g. bond/O.H in which O and H are atom types. <system> is a gaussian fchk file or molmod chk file."""
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
    parser.add_option(
        '--relative-amplitude', default=0.05, type=float,
        help='Defines the relative amplitude of the change in the ic value relative to the ic equilibrium value. [default=%default]'
    )
    parser.add_option(
        '--eunit', default='kjmol', type=str,
        help='The energy unit used. [default=%default]'
    )
    parser.add_option(
        '--qunit', default='au', type=str,
        help='The ic unit used. [default=%default]'
    )
    parser.add_option(
        '--kunit', default='None', type=str,
        help='The unit used for the force constant. [default=eunit/qunit^2]'
    )
    options, args = parser.parse_args()
    if not len(args)==2:
        raise IOError('Too many input argumets: expected 2, recieved %i' %(len(args)))
    icname = args[0]
    fn_sys = args[1]
    options.spring = options.spring*kjmol/angstrom**2
    if options.kunit=='None':
        options.kunit = '%s/%s**2' %(options.eunit, options.qunit)
    return icname, fn_sys, options

def main():
    icname, fn_sys, options = parser()
    system = System('sys', fn_sys)
    if options.psf is not None:
        system.topology_from_psf(options.psf)
    system.find_ic_patterns([icname])    
    plot_perturbation(system, icname, coupling=options.coupling, free_depth=options.free_depth, spring=options.spring)

if __name__=='__main__':
    main()  