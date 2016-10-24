from molmod.units import angstrom, kjmol

from common import log, read_system

from nose import SkipTest

import numpy as np

def do_taylor(name, ntests=10, eps=1e-3*angstrom, gtol=1e-6*kjmol/angstrom, htol=1e-9*kjmol/angstrom**2):
    with log.section('PROGRAM', 2):
        system, ref = read_system(name)
    print ''
    print 'Checking if analytic (A) and numerical (N) GRADIENT are eqaul.'
    print 'The used criterian is that RMSD(A-D)<%.3e kJ/(mol.A).' % (gtol/(kjmol/angstrom))
    print 'The values given below are in kJ/(mol.A).'
    for i in xrange(ntests):
        npos = system.pos + np.random.normal(0.0, 0.1, system.pos.shape)*angstrom
        ganalytic = ref.gradient(npos)
        gnumerical = np.zeros(3*len(system.numbers), float)
        for j in xrange(3*len(system.numbers)):
             dx = np.zeros(3*len(system.numbers), float)
             dx[j] = 1.0
             dx = dx.reshape(npos.shape)
             eplus = ref.energy(npos + eps*dx)
             emin = ref.energy(npos - eps*dx)
             gnumerical[j] = (eplus-emin)/(2.0*eps)
        std = ganalytic.std()
        rmsd = np.sqrt(((ganalytic-gnumerical.reshape(ganalytic.shape))**2).mean())
        print '  STD(A)=%.6e    RMSD(A-N)=%.6e' %(std/(kjmol/angstrom), rmsd/(kjmol/angstrom))
        assert rmsd<gtol
    print 'Checking if analytic (A) and numerical (N) HESSIAN are eqaul.'
    print 'The used criterian is that both STD(N-N.T) and RMSD(A-D) are smaller than %.3e kJ/(mol.A^2).' % (htol/(kjmol/angstrom**2))
    print 'The values given below are in kJ/(mol.A^2).'
    for i in xrange(ntests):
        npos = system.pos + np.random.normal(0.0, 0.1, system.pos.shape)*angstrom
        hanalytic = ref.hessian(npos)
        hnumerical = np.zeros([3*len(system.numbers), 3*len(system.numbers)], float)
        for k in xrange(3*len(system.numbers)):
            dx = np.zeros(3*len(system.numbers), float)
            dx[k] = 1.0
            dx = dx.reshape(npos.shape)
            gplus = ref.gradient(npos + eps*dx)
            gmin = ref.gradient(npos - eps*dx)
            hnumerical[k,:] = (gplus-gmin).reshape(3*len(system.numbers))/(2.0*eps)
        symm = (hnumerical-hnumerical.T).std()
        std = hanalytic.std()
        rmsd = np.sqrt(((hanalytic-hnumerical.reshape(hanalytic.shape))**2).mean())
        print '  STD(A)=%.6e    STD(N-N.T)=%.6e    RMSD(A-N)=%.6e: ' %(std/(kjmol/angstrom**2), rmsd/(kjmol/angstrom**2), rmsd/(kjmol/angstrom**2))
        assert symm<htol
        assert rmsd<htol

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
