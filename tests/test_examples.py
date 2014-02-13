from molmod.units import *
from molmod.constants import lightspeed
from molmod.periodic import periodic as pt
from molmod.ic import bond_length
from molmod.io.xyz import XYZWriter

from quickff.system import System
from quickff.model import Model
from quickff.context import context
from quickff.program import Program

import numpy as np, os

def test_h2():
    #frequency of H2 stretch mode in gaussian.fchk calculation is 4416.921/cm
    #and an equilibrium bond length of 0.7442380 A this test checks if the
    #force field predicts the same values
    r0 = 0.7442380*angstrom
    freq = (2*np.pi*lightspeed)*4416.921/centimeter
    mass = pt['H'].mass/2 #reduced mass for the H2 stretch mode
    print 'AI     :    K = %.3f kjmol/A^2    q0 = %.6f A' %(
        mass*freq**2/(kjmol/angstrom**2), r0/angstrom
    )
    #Load system, model and pert. theory
    fn_fchk = os.path.join(context.get_fn('systems/H2'), 'gaussian.fchk')
    system = System.from_files([fn_fchk])
    system.guess_ffatypes('low')
    system.determine_ics_from_topology()
    model = Model.from_system(system, ai_project=False, ei_pot_kind='Zero', vdw_pot_kind='Zero')
    program = Program(system, model)
    #estimate pars from perturbation theory
    trajectories = program.generate_trajectories(verbose=False)
    fftab = program.estimate_from_pt(trajectories, verbose=False)
    k, q0 = fftab['bond/H.H']
    print 'FF (pt):    K = %.3f kjmol/A^2    q0 = %.6f A' %(
        k/(kjmol/angstrom**2), q0/angstrom
    )
    assert abs(k/(mass*freq**2)-1.0) < 1e-3
    assert abs(q0/r0-1.0) < 1e-3
    #refine pars using Hessian cost function
    fftab = program.refine_cost(verbose=False)
    k, q0 = fftab['bond/H.H']
    print 'FF (rc):    K = %.3f kjmol/A^2    q0 = %.6f A' %(
        k/(kjmol/angstrom**2), q0/angstrom
    )
    assert abs(k/(mass*freq**2)-1.0) < 1e-3
    assert abs(q0/r0-1.0) < 1e-3
