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

from quickff.tools import guess_ffatypes, get_ei_radii, average, charges_to_bcis
from quickff.io import read_abinitio, make_yaff_ei, read_bci_constraints
from yaff import System
from optparse import OptionParser

import h5py

def parse():
    usage = '%prog [options] fn_sys fn_hdf5 path'
    description  = 'This script will read the atomic charges from the input '
    description += 'arguments and make the Yaff parameters file for use with '
    description += 'the QuickFF option --ei. fn_sys is any file from which the '
    description += 'system can be extracted (MolMod CHK, Gaussian FCHK, XYZ, '
    description += '...). Path represents the path to the charges, which are '
    description += 'assumed to be stored in the file fn_hdf5 in the dataset '
    description += '/path/charges.'
    parser = OptionParser(usage=usage, description=description)
    parser.add_option(
        '--ffatypes', default=None,
        help='Assign atom types in the system. If FFATYPES is a string equal to '
             'either low,medium, high or highest, atom types will be assigned '
             'through the automatic built-in detection (see documentation). If '
             'FFATYPES is None, the atom types are assumed to be defined in the '
             'input files. [default=%default]'
    )
    parser.add_option(
        '--gaussian', default=False, action='store_true',
        help='Use gaussian smeared charges. The radii are taken from the '
             'dataset /path/radii if it exists, otherwise the radii are '
             'estimated according to the procedure of Chen et al. See '
             '``quickff.tools.get_ei_radii`` for more info.'
    )
    parser.add_option(
        '--bci', default=False, action='store_true',
        help='Convert averaged atomic charges to bond charge increments, i.e. '
             'charge transfers along the chemical bonds in the system. In this '
             'way, one is certain of ending up with a globally neutral system '
             'even if bcis from different systems are combined. This option '
             'requires the definition of bonds, if the bonds cannot be read '
             'from the system file, they will be estimated from the interatomic '
             'distances using the detect_bonds routine in the Yaff System '
             'class. [default=%default]'
    )
    parser.add_option(
        '--bci-constraints', default=None,
        help='A file containing constraints for the charge to bci fit in a '
             'master: slave0,slave1,...: sign format. A new line should be used '
             'for each master and the format is insensitive towards spaces.'
             'Sign should be 1.0 or -1.0 indicating wheter or not a sign switch '
             'should be introduced when mapping the slaves to the master.'
    )
    parser.add_option(
        '-o', '--output', default=None,
        help='Name of the Yaff file to write the parameters to.'
             ' [default=pars_ei_${PATH}.txt]'
    )
    options, args = parser.parse_args()
    assert len(args)==3, 'Exactly three arguments are required: a system file, a HDF5 file and the path to the charges.'
    return options, args

def main():
    options, args = parse()
    fn_sys, fn_hdf5, path = args
    if fn_sys.endswith('.fchk'):
        numbers, coords, energy, grad, hess, masses, rvecs, pbc = read_abinitio(fn_sys)
        system = System(numbers, coords, rvecs=None, charges=None, radii=None, masses=masses)
        system.detect_bonds()
    else:
        system = System.from_file(fn_sys)
    if options.ffatypes is not None:
        guess_ffatypes(system, options.ffatypes)
    ffatypes = [system.ffatypes[i] for i in system.ffatype_ids]
    h5 = h5py.File(fn_hdf5)
    if not path in h5 or 'charges' not in h5[path]:
        raise IOError('Given file %s does not contain dataset %s/charges' %(fn_hdf5, path))
    radii = None
    if options.gaussian:
        if 'radii' in h5[path]:
            radii = average(h5['%s/radii' %path][:], ffatypes, fmt='dict')
        else:
            radii = average(get_ei_radii(system.numbers), ffatypes, fmt='dict')
    if options.output is None:
        fn_out = 'pars_ei_%s.txt' %path.replace('/', '_')
    else:
        fn_out = options.output
    if options.bci:
        constraints = {}
        if options.bci_constraints is not None:
            constraints = read_bci_constraints(options.bci_constraints)
        charges = average(h5['%s/charges' %path][:], ffatypes, fmt='full')
        bcis = charges_to_bcis(charges, ffatypes, system.bonds, constraints=constraints)
        make_yaff_ei(fn_out, None, bcis=bcis, radii=radii)
    else:
        charges = average(h5['%s/charges' %path][:], ffatypes, fmt='dict')
        make_yaff_ei(fn_out, charges, radii=radii)

if __name__=='__main__':
    main()
