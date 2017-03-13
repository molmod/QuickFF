#! /usr/bin/env python
# -*- coding: utf-8 -*-
# QuickFF is a code to quickly derive accurate force fields from ab initio input.
# Copyright (C) 2012 - 2016 Louis Vanduyfhuys <Louis.Vanduyfhuys@UGent.be>
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

from argparse import ArgumentParser
import os

import h5py as h5

from molmod.io.chk import load_chk
from yaff import System

from quickff.tools import guess_ffatypes, get_ei_radii, average, charges_to_bcis
from quickff.io import read_abinitio, make_yaff_ei, read_bci_constraints


description  = '''\
This script reads atomic charges from an input file and makes a Yaff parameters file
suitable for the QuickFF option --ei.'''


def parse_args():
    parser = ArgumentParser(description=description)
    parser.add_argument(
        '-v', '--verbose', default=False, action='store_true',
        help='Increase verbosity of the script [default=%(default)s].'
    )
    parser.add_argument(
        '--ffatypes', default=None,
        choices=['None','list_of_atypes', 'low','medium','high','highest'],
        help='Assign atom types in the system by parsing an ordered list of'
             'atom types as argument or through the automatic built-in '
             'detection (see documentation). By default (or if None is given), '
             'the atom types are assumed to be defined in the input files. '
             '[default=%(default)s]'
    )
    parser.add_argument(
        '--gaussian', default=False, action='store_true',
        help='Use gaussian smeared charges. The radii are taken from the input '
             'file fn_in (from dataset /path/radii for HDF5 file or from label '
             '`radii` for CHK file) if the data is present, otherwise the radii '
             'are estimated according to the procedure of Chen et al. See '
             '``quickff.tools.get_ei_radii`` for more info.'
    )
    parser.add_argument(
        '--bci', default=False, action='store_true',
        help='Convert averaged atomic charges to bond charge increments, i.e. '
             'charge transfers along the chemical bonds in the system. In this '
             'way, one is certain of ending up with a globally neutral system '
             'even if bcis from different systems are combined. This option '
             'requires the definition of bonds, if the bonds cannot be read '
             'from the system file, they will be estimated from the interatomic '
             'distances using the detect_bonds routine in the Yaff System '
             'class. [default=%(default)s]'
    )
    parser.add_argument(
        '--bci-constraints', default=None,
        help='A file containing constraints for the charge to bci fit in a '
             'master: slave0,slave1,...: sign format. A new line should be used '
             'for each master and the format is insensitive towards spaces.'
             'Sign should be 1.0 or -1.0 indicating wheter or not a sign switch '
             'should be introduced when mapping the slaves to the master.'
    )
    parser.add_argument(
        'fn_sys',
        help='Any file from which the system can be extracted (MolMod CHK, Gaussian '
             'FCHK, XYZ, ...).'
    )
    parser.add_argument(
        'charges',
        help='The atomic charges to be used. This argument has the form fn_charges:path, '
             'where fn_charges is a filename that contains the charges and path refers '
             'to a location within that file where the charges can be found. If '
             'fn_charges is an HDF5 file, path is the location of a dataset containing '
             'the charges, e.g. \'/charges\'. If fn_charges is a MolMod CHK file, path '
             'is the label of the dataset that contains the atomic charges.'
    )
    parser.add_argument(
        'fn_out', default='pars_ei.txt', nargs='?',
        help='Name of the Yaff file to write the parameters to. [default=%(default)s]'
    )
    args = parser.parse_args()
    if args.charges.count(':') != 1:
        parser.error('The argument charges must contain exactly one colon.')
    return args


def main():
    args = parse_args()

    # Load system file
    if args.fn_sys.endswith('.fchk'):
        numbers, coords, energy, grad, hess, masses, rvecs, pbc = read_abinitio(args.fn_sys, do_hess=False)
        system = System(numbers, coords, rvecs=None, charges=None, radii=None, masses=masses)
        system.detect_bonds()
    else:
        system = System.from_file(args.fn_sys)

    # Guess atom types if needed
    if args.ffatypes is not None:
        guess_ffatypes(system, args.ffatypes)
    ffatypes = [system.ffatypes[i] for i in system.ffatype_ids]

    # Load atomic charges
    fn_charges, _, path = args.charges.partition(':')
    if fn_charges.endswith('.h5'):
        with h5.File(fn_charges, 'r') as f:
            if not path in f:
                raise IOError('Given HDF5 file %s does not contain a dataset %s' % (fn_charges, path))
            charges = f[path][:]
            radii = None
            if args.gaussian:
                path_radii = os.path.join(os.path.dirname(path), 'radii')
                if 'radii' in f[path]:
                    radii = average(f['%s/radii' %path][:], ffatypes, fmt='dict')
                else:
                    radii = average(get_ei_radii(system.numbers), ffatypes, fmt='dict')
    elif fn_charges.endswith('.chk'):
        sample = load_chk(fn_charges)
        if path in sample.keys():
            charges = sample[path]
        else:
            raise IOError('Given CHK file %s does not contain a dataset with label %s' % (fn_charges, path))
        radii = None
        if args.gaussian:
            if 'radii' in sample.keys():
                radii = average(sample['radii'], ffatypes, fmt='dict')
    else:
        raise IOError('Invalid extension, fn_charges should be a HDF5 or a CHK file.')

    # Derive charge parameters
    if args.bci:
        constraints = {}
        if args.bci_constraints is not None:
            constraints = read_bci_constraints(args.bci_constraints)
        bcis = charges_to_bcis(charges, ffatypes, system.bonds, constraints=constraints, verbose=args.verbose)
        make_yaff_ei(args.fn_out, None, bcis=bcis, radii=radii)
    else:
        charges = average(charges, ffatypes, fmt='dict', verbose=args.verbose)
        make_yaff_ei(args.fn_out, charges, radii=radii)


if __name__=='__main__':
    main()
