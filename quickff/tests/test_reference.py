from molmod.units import angstrom, kjmol

from common import log, read_system

from nose import SkipTest

import numpy as np

def do_taylor(name, ntests=10, eps=1e-4, gtol=1e-3*kjmol/angstrom, htol=1e-3*kjmol/angstrom**2):
    with log.section('NOSETST', 2):
        system, ref = read_system(name)
    print ''
    for i in xrange(ntests):
        npos = system.pos + np.random.normal(0.0, 0.1, system.pos.shape)*angstrom
        gref = ref.gradient(npos)
        gnum = np.zeros(3*len(system.numbers), float)
        for j in xrange(3*len(system.numbers)):
             dx = np.zeros(3*len(system.numbers), float)
             dx[j] = 1.0
             dx = dx.reshape(npos.shape)
             eplus = ref.energy(npos + eps*dx)
             emin = ref.energy(npos - eps*dx)
             gnum[j] = (eplus-emin)/(2.0*eps)
        std = gref.std()
        M = (abs(gref-gnum.reshape(gref.shape))).max()
        print '  Sample %i:    STD=%.6e    MaxDev=%.6e' %(i, std/(kjmol/angstrom), M/(kjmol/angstrom))
        assert M<gtol
    for i in xrange(ntests):
        npos = system.pos + np.random.normal(0.0, 0.1, system.pos.shape)*angstrom
        href = ref.hessian(npos)
        hnum = np.zeros([3*len(system.numbers), 3*len(system.numbers)], float)
        for k in xrange(3*len(system.numbers)):
            dx = np.zeros(3*len(system.numbers), float)
            dx[k] = 1.0
            dx = dx.reshape(npos.shape)
            gplus = ref.gradient(npos + eps*dx)
            gmin = ref.gradient(npos - eps*dx)
            hnum[k,:] = (gplus-gmin).reshape(3*len(system.numbers))/(2.0*eps)
        hnum = 0.5*(hnum+hnum.T)
        std = href.std()
        M = (abs(gref-gnum.reshape(gref.shape))).max()
        print '  Sample %i:    STD=%.6e    MaxDev=%.6e: ' %(i, std/(kjmol/angstrom**2), M/(kjmol/angstrom**2))
        assert M<htol

def test_taylor_water():
    do_taylor('water/gaussian.fchk')

def test_taylor_methane():
    do_taylor('methane/gaussian.fchk')

def test_taylor_ethene():
    do_taylor('ethene/gaussian.fchk')

def test_taylor_ethanol():
    do_taylor('ethanol/gaussian.fchk')

def test_taylor_amoniak():
    do_taylor('amoniak/gaussian.fchk')

def test_taylor_benzene():
    do_taylor('benzene/gaussian.fchk')
