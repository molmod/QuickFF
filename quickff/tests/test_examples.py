from molmod.units import *
from molmod.constants import lightspeed
from molmod.periodic import periodic as pt

from quickff.tools import guess_ffatypes
from quickff.program import DeriveNonDiagFF
from quickff.context import context

from yaff import System

import numpy as np, os

from common import log, read_system

def test_h2():
    #frequency of H2 stretch mode in gaussian.fchk calculation is 4416.656/cm
    #and an equilibrium bond length of 0.7442380 A. This test checks if the
    #force field predicts the same values
    r0 = 0.7442380*angstrom
    freq = (2*np.pi)*4416.65640485*lightspeed/centimeter
    mass = pt['H'].mass/2 #reduced mass for the H2 stretch mode
    #Load system, model and pert. theory and estimate ff
    with log.section('PROGRAM', 2):
        system, ai = read_system('H2/gaussian.fchk')
        guess_ffatypes(system, 'low')
        program = DeriveNonDiagFF(system, ai)
        program.do_pt_generate()
        program.do_pt_estimate()
        K_pt, rv_pt = program.valence.get_params(0, only='all')
        program.do_hc_estimatefc(['HC_FC_DIAG'])
        K_hc, rv_hc = program.valence.get_params(0, only='all')
    #print results
    print ''
    print 'AI     :    K = %.3f kjmol/A^2    q0 = %.6f A' %(mass*freq**2/(kjmol/angstrom**2), r0/angstrom)
    print 'FF (PT):    K = %.3f kjmol/A^2    q0 = %.6f A' %(K_pt/(kjmol/angstrom**2), rv_pt/angstrom)
    print 'FF (HC):    K = %.3f kjmol/A^2    q0 = %.6f A' %(K_hc/(kjmol/angstrom**2), rv_hc/angstrom)
    print ''
    #perform assertion checks
    assert abs(K_pt/(mass*freq**2)-1.0) < 1e-3
    assert abs(rv_pt/r0-1.0) < 1e-3
    assert abs(K_hc/(mass*freq**2)-1.0) < 1e-3
    assert abs(rv_hc/r0-1.0) < 1e-3
    assert abs(K_hc/K_pt-1.0) < 1e-6
    assert abs(rv_hc/rv_pt-1.0) < 1e-6
