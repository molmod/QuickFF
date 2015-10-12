#!/usr/bin/env python
# -*- coding: utf-8 -*-
# QuickFF is a code to quickly derive accurate force fields from ab initio input.
# Copyright (C) 2012 - 2015 Louis Vanduyfhuys <Louis.Vanduyfhuys@UGent.be>
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
#!/usr/bin/env python
import numpy as np, os

from molmod.units import *
from molmod.constants import lightspeed
from molmod.periodic import periodic as pt
from molmod.io import FCHKFile
from yaff import System

from quickff.context import context
from quickff.program import Program
from quickff.tools import guess_ffatypes, read_abinitio
from quickff.refdata import ReferenceData

def test_water():
    fn_fchk = os.path.join(context.get_fn('systems/water'), 'gaussian.fchk')
    numbers, coords, energy, grad, hess, masses, rvecs, pbc = read_abinitio(fn_fchk)
    system = System(numbers, coords)
    guess_ffatypes(system, 'highest')
    refdata = ReferenceData(coords, energy, grad, hess)
    program = Program(system, refdata)
    trajectories = program.generate_trajectories()
    fftab = program.estimate_from_pt(trajectories, verbose=False)
    fftab = program.refine_cost(fftab, verbose=False)
