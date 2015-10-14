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
from nose.tools import assert_raises
import time

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

def test_boxqp_wrong_boundaries():
    n = 2
    A = np.zeros((2,2))
    B = np.zeros((2))
    bndl = np.array([1.0,5.0])
    bndu = np.array([3.0,4.0])
    x0 = np.zeros((2))
    with assert_raises(AssertionError):
        boxqp(A, B, bndl, bndu, x0)

def test_boxqp_no_boundaries():
    A = np.array([[3.0,1.0],[1.0,2.0]])
    B = np.array([5.0,6.0])
    bndl = np.array([np.NINF, np.NINF])
    bndu = np.array([np.inf, np.inf])
    x0 = np.array([-2.0,8.0])
    # BOXQP solution
    xA = boxqp(A,B,bndl, bndu, x0)
    # Exact solution
    xB = np.dot(np.linalg.inv(A),B)
    assert np.all(np.abs(xA-xB)<1e-6)

def test_boxqp_check_boundaries():
    A = np.array([[3.0,1.0],[1.0,2.0]])
    B = np.array([5.0,6.0])
    bndl = np.array([3.2,-1.0])
    bndu = np.array([8.5,np.inf])
    x0 = np.array([-2.0,8.0])
    x = boxqp(A,B,bndl, bndu, x0)
    assert np.all(bndl<=x)
    assert np.all(bndu>=x)

def test_boxqp():
    print "%7s %7s %7s | %10s  %10s  | %10s  %10s   | %10s   %10s  |%10s  %10s" % \
            ("n","na(x*)","na(x0)","RMSD BOXQP", "RMSD SLSQP","DEV BOXQP", "DEV SLSQP","TIME BOXQP", "TIME SLSQP","NIT BOXQP", "NIT SLSQP")
    print "-"*127
    # Uncomment following line to perform comparison of SLSQP minimizer and BOXQP
    # Takes quite a long time...
    #for n in [20,50,100,200,500]:
    for n in [20,50]:
        for i in xrange(3):
            na_sol = n/10 + 2*n/5*i
            for j in xrange(3):
                na_init = n/10 + 2*n/5*j
                check_boxqp(n,na_sol,na_init)
    #assert False

def check_boxqp(n, na_sol, na_init, ncond=6, ndeg=3):
    '''Check the solution of our boxqp implementation by generating a random
    problem. A comparison with SciPy SLSQP is also provided
    '''
    # Construct the problem
    A, B, bndl, bndu, x0, sol = construct_boxqp(n, na_sol, na_init)
    # Use our minimizer
    t0 = time.time()
    x, nit = boxqp(A, B, bndl, bndu, x0, status=True)
    t1 = time.time()
    # Use conventional SLSQP SciPy minimizer
    def fun(x):
        return 0.5*np.dot(np.dot(x,A),x) - np.dot(x,B)
    def jac(x):
        return np.dot(A,x) - B
    bounds = np.array([(bndl[i],bndu[i]) for i in xrange(n)])
    from scipy.optimize import minimize
    t2 = time.time()
    result = minimize(
            fun, x0, jac=jac, method='SLSQP', bounds=bounds,
            tol=1e-8, options={'disp': False}
        )
    t3 = time.time()
    rmsd0 = np.sqrt(np.mean((x-sol)**2))
    rmsd1 = np.sqrt(np.mean((result.x-sol)**2))
    dev0 = fun(x) - fun(sol)
    dev1 = fun(result.x) - fun(sol)
    print "%7d %7d %7d |   %.1e   %.1e     |   %+.1e   %+.1e    |    %7.1f      %7.1f  |  %8d    %8d" % \
        (n,na_sol,na_init,rmsd0, rmsd1, dev0, dev1, t1-t0, t3-t2, nit,result.nit)
    assert np.all(np.abs(x-sol)<1e-6)
    assert np.all(x>=bndl)
    assert np.all(x<=bndu)

def construct_boxqp(n, na_sol, na_init, ncond=8, ndeg=3):
    '''Construct a randomly generated box-constrained quadratic programming
    problem, cfr 10.1007/s00211-004-0569-y
    '''
    assert na_sol<=n
    assert na_init<=n
    # Construct the A matrix
    P = np.eye(n)
    for i in xrange(3):
        w = np.random.uniform(size=n)
        w /= np.linalg.norm(w)
        P = np.dot(P,np.eye(n)-2.0*np.outer(w,w))
    D = np.diag(np.exp(np.linspace(0,n-1,n)*ncond/(n-1)))
    A = np.dot(np.dot(P,D),P.transpose())
    # Check that A is hermitian (symmetric)
    assert np.all(np.abs(A-A.transpose())<1e-10)
    # Check that A is positive definite using Sylvester's criterion (takes a long time)
    #for i in xrange(n):
    #    assert np.linalg.det(A[:i+1,:i+1]) > 0.0
    # Check that A has the correct condition number
    assert np.abs(np.linalg.cond(A)-np.exp(ncond))<1e-6
    # Construct solution
    x = (np.random.uniform(size=n)-0.5)*2.0
    # Select active set for both solution and initial guess
    active_sol = np.random.uniform(size=n) <= 1.*na_sol/n
    active_init = np.random.uniform(size=n) <= 1.*na_init/n
    # Construct boundaries
    bndl = -np.ones(n)
    bndu = np.ones(n)
    r = np.zeros(n)
    low = np.random.uniform(size=active_sol.sum())>0.5
    high = np.logical_not(low)
    bndl[active_sol][low] = x[active_sol][low]
    r[active_sol][low] = np.power(10,-np.random.uniform(size=low.sum())*ndeg)
    bndu[active_sol][high] = x[active_sol][high]
    r[active_sol][high] = np.power(10,-np.random.uniform(size=high.sum())*ndeg)
    # Construct initial guess
    x0 = 0.5*(bndl+bndu)
    low = np.random.uniform(size=active_init.sum())>0.5
    high = np.logical_not(low)
    x0[active_init][low] = bndl[active_init][low]
    x0[active_init][high] = bndu[active_init][high]
    B = np.dot(A,x) - r
    return A, B, bndl, bndu, x0, x
