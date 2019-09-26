# -*- coding: utf-8 -*-
# QuickFF is a code to quickly derive accurate force fields from ab initio input.
# Copyright (C) 2012 - 2019 Louis Vanduyfhuys <Louis.Vanduyfhuys@UGent.be>
# Steven Vandenbrande <Steven.Vandenbrande@UGent.be>,
# Jelle Wieme <Jelle.Wieme@UGent.be>,
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
from contextlib import contextmanager
from glob import glob
import os
import shutil
import tempfile

import h5py as h5

from quickff.io import read_abinitio

from quickff.io import read_abinitio
from quickff.reference import SecondOrderTaylor
from quickff.valence import ValenceFF

from yaff import System

from quickff.log import log
log.set_level('silent')

try:
    from importlib.resources import path
except ImportError:
    from importlib_resources import path


__all__ = ['log', 'read_system', 'tmpdir']

def read_system(name):
    # Load system data
    dn = 'quickff.data.systems'
    if '/' in name:
        words = name.split('/')
        dn += '.%s' %('.'.join(words[:-1]))
        name = words[-1]
    with path(dn, name) as fn:
        numbers, coords, energy, grad, hess, masses, rvecs, pbc = read_abinitio(fn)
        fns_wpart = glob(os.path.join(os.path.dirname(fn), 'gaussian_mbis.h5'))
    # Try to load charges.
    charges = None
    if len(fns_wpart) > 0:
        with h5.File(fns_wpart[0], 'r') as f:
            charges = f['charges'][:]
    # Create system object.
    system = System(numbers, coords, charges=charges)
    system.detect_bonds()
    system.set_standard_masses()
    # Load ab initio data.
    ai = SecondOrderTaylor('ai', coords=system.pos.copy(), energy=energy, grad=grad, hess=hess, pbc=pbc)
    return system, ai

@contextmanager
def tmpdir(name):
    """Create a temporary directory where output files of a test can be written to.

    The name argument can be used to obtained a recognizable test directory name. The
    temporary directory and its contents are automatically removed.
    """
    dn = tempfile.mkdtemp(name)
    try:
        yield dn
    finally:
        shutil.rmtree(dn)
