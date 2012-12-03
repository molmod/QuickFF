#! /usr/bin/env python

import numpy as np
from optparse import OptionParser
from molmod.units import *
from molmod.io.chk import dump_chk

from quickff.system import System
from quickff.estimater import estimate
from quickff.fftable import FFTable

def parser():
    usage = "%prog [options] fchk"
    description = """This script is part of QuickFF. It converts a gaussian fchk file into a molmod chk checkpoint file."""
    parser = OptionParser(usage=usage, description=description)
    parser.add_option(
        '--psf', default=None,
        help='A MolMod PSF file from which the topology and atom types are taken. If this argument is not present, atom symbols are used as atom types and the topology is estimated based on the geometry.'
    )
    parser.add_option(
        '--out', default='system.chk',
        help='The name of the MolMod chk output file. [default=%default]'
    )
    options, args = parser.parse_args()
    fn_fchk = args[0]
    return fn_fchk, options

def main():
    fn_fchk, options = parser()
    sample = System.load_sample(fn_fchk)
    if options.psf is not None:
        System.topology_from_psf(sample, options.psf)
    dump_chk(options.out, sample)

if __name__=='__main__':
    main()
