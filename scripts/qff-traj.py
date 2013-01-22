#! /usr/bin/env python

import numpy as np, matplotlib.pyplot as pp
from optparse import OptionParser
from molmod.units import *

from quickff.system import System
from quickff.perturbation import *
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
        '--cost-ic', default=1.0, type=float,
        help='The weight of the ic cost in the total cost function for calculating the perturbation on the geometry. The ic cost expresses for a given geometry the deviation of the internal coordinates from their equilibrium value. [default = %default]'
    )
    parser.add_option(
        '--cost-energy', default=1.0, type=float,
        help='The weight of the energy cost in the total cost function for calculating the perturbation on the geometry. The energy is calculated as the second order Taylor expression using the forces and hessian in equilibrium [default = %default]'
    )
    parser.add_option(
        '--relative-amplitude', default=0.10, type=float,
        help='Defines the relative amplitude of the change in the ic value relative to the ic equilibrium value. [default=%default]'
    )
    parser.add_option(
        '--ei-rule', default=-1, type=int,
        help='Defines the exclusion rule for the electrostatic interactions. If set to -1, no electrostatic interactions are taken into account during force field fitting. If the rule is set to 0,1,2,3: pairs seperated by less then or equal to 1,2,3 bonds respectively are excluded from ei interactions. [default=10]'
    )
    parser.add_option(
        '--charges', default=None, 
        help='Specify the charges of the system, this is necessary if a ei-rule is specified different then -1. The charges are comma seperated and the order should be identical to the order in the sample file. If the charges are not specified but ei-rule is larger then -1, the charges are taken from the psf if present, or from the sample file if present. Alternatively, a hipart charge txt file can also be used to specify the charges. [default=%default]'
    )
    parser.add_option(
        '--atypes-level', default='medium', 
        help='Overwrite the atom types according to level ATYPES_LEVEL. Low will choose atom types based only on atom number, medium will choose atom types based on local topology and high will choose atom types based on atom index. [default=%default]'
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
        raise IOError('Too many input argumets: expected 2, recieved %i: ' %(len(args))), args
    icname = args[0]
    fn_chk = args[1]
    if options.kunit=='None':
        options.kunit = '%s/%s**2' %(options.eunit, options.qunit)
    return icname, fn_chk, options

def main():
    icname, fn_chk, options = parser()
    system = System(fn_chk, fn_psf=options.psf, guess_atypes_level=options.atypes_level, charges=options.charges)
    system.define_models(eirule=options.ei_rule)
    system.find_ic_patterns(['all']) 
    pt = RelaxedGeometryPT(energy_penalty=options.cost_energy, ic_penalty=options.cost_ic, dq_rel=options.relative_amplitude)
    pt.plot_single(system, icname, qunit=options.qunit, kunit=options.kunit, eunit=options.eunit)

if __name__=='__main__':
    main()  
