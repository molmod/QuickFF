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

import numpy as np

from quickff.tools import *

def test_fitpar():
    x = np.linspace(-5.0,1.0,15)
    pars0 = [3.2,1.2,0.5]
    y = pars0[0]*x*x + pars0[1]*x + pars0[2]
    pars1 = fitpar(x,y)
    assert np.all(np.abs(pars1-pars0)<1e-10)

def test_statistics():
    N0 = 1e7
    mu0 = 3.2
    sigma0 = 2.6
    x = np.random.normal(mu0, sigma0, N0)
    mu1, sigma1, N1 = statistics(x)
    assert N0==N1
    assert np.abs(sigma1/sigma0-1.0)<1e-3
    assert np.abs(mu1/mu0-1.0)<1e-3
