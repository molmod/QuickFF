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

import os

from quickff.io import *
from quickff.context import *

def test_mil53al_np():
    fn_vasprun = os.path.join(context.get_fn('systems/mil53al_np'), 'vasprun.xml')
    # Load basic information
    vasprun = VASPRun(fn_vasprun)
    field_labels = vasprun.fields.keys()
    assert len(field_labels)==5
    for label in ['masses', 'pos_init', 'numbers', 'rvecs_init', 'energies']:
        assert label in field_labels
    # Request loading more information
    vasprun = VASPRun(fn_vasprun, field_labels=['hessian','gradient'])
    field_labels = vasprun.fields.keys()
    assert len(field_labels)==7
    for label in ['masses', 'pos_init', 'numbers', 'rvecs_init', 'energies', 'hessian', 'gradient']:
        assert label in field_labels
    natom = len(vasprun.fields['numbers'])
    assert natom == 76
    assert vasprun.fields['hessian'].shape[0] == 3*natom
    assert vasprun.fields['hessian'].shape[1] == 3*natom
    # VASP computes the gradient for every displacement
    assert vasprun.fields['gradient'].shape[0] == 6*natom+1
    assert vasprun.fields['gradient'].shape[1] == natom
    assert vasprun.fields['gradient'].shape[2] == 3
    assert vasprun.fields['masses'].shape[0] == natom
    assert vasprun.fields['rvecs_init'].shape[0] == 3
    assert vasprun.fields['rvecs_init'].shape[1] == 3
    assert vasprun.fields['pos_init'].shape[0] == natom
    assert vasprun.fields['pos_init'].shape[1] == 3
