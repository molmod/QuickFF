from molmod.units import angstrom, kjmol
import numpy as np

from quickff.model import *
from common import get_system, get_scaled_pairs


def run_taylor(molecule, pot):
    tolE_abs = 5e-2
    system = get_system(molecule, ei_path='/wpart/he')
    scaled_pairs = get_scaled_pairs(system)
    coords0 = system.ref.coords
    if pot.lower() in ['coulomb', 'ei', 'electrostatic']:
        exact = CoulombPot(system.charges, [1.0, 1.0, 1.0], scaled_pairs, coords0=coords0)
    elif pot.lower() in ['lennardjones', 'lj', 'vdw', 'vanderwaals']:
        system.read_uff_vdw()
        exact = LennardJonesPot(system.sigmas, system.epsilons, [0.0,0.0,1.0], scaled_pairs, coords0=coords0)
    else:
        raise ValueError('Invalic potential specification: %s' %pot)
    harm = HarmonicPot(coords0, 0.0, exact.calc_gradient(coords0), exact.calc_hessian(coords0))
    dxs = np.random.normal(0.0, 5e-3*angstrom, [100, 3*len(coords0)])
    exactEs = np.zeros(dxs.shape[0], float)
    harmEs = np.zeros(dxs.shape[0], float)
    for i, dx in enumerate(dxs):
        exactEs[i] = exact.calc_energy(coords0+dx.reshape([-1,3]))
        harmEs[i] = harm.calc_energy(coords0+dx.reshape([-1,3]))
    meanE = exactEs.mean()
    stdE = exactEs.std()
    print '     Deformation [A]               |   Energy [kJ/mol] (mean=%.3e  std=%.3e)' %(meanE/kjmol, stdE/kjmol)
    print '    ===================================================================================='
    print '     i    rms          max         |   exact        taylor      error_abs    error_rel  '
    print '    ------------------------------------------------------------------------------------'
    maxs = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    for i, dx in enumerate(dxs):
        current = [
            np.linalg.norm(dx)/np.sqrt(3*len(coords0)), max(abs(dx)), 
            exactEs[i], harmEs[i],
            abs(harmEs[i]-exactEs[i]), abs(harmEs[i]-exactEs[i])/abs(exactEs[i])
        ]
        print '    %2i    %.3e    %.3e   |  % .6f    % .6f    %.3e    %.3e   ' %(
            i, current[0]/angstrom, current[1]/angstrom, current[2]/kjmol, current[3]/kjmol, current[4]/kjmol, current[5]
        )
        maxs = np.array([max(prev, abs(curr)) for prev, curr in zip(maxs, current)])
    print '    ------------------------------------------------------------------------------------'
    print '    max   %.3e    %.3e   |  % .6f    % .6f    %.3e    %.3e   ' % (
        maxs[0]/angstrom, maxs[1]/angstrom, maxs[2]/kjmol, maxs[3]/kjmol, maxs[4]/kjmol, maxs[5]
    )
    print '    ===================================================================================='
    for i, dx in enumerate(dxs):
        assert abs(exactEs[i]-harmEs[i])<tolE_abs*stdE


def test_ei_water():
    run_taylor('water', 'ei')

def test_ei_ethanol():
    run_taylor('ethanol', 'ei')

def test_lj_ethanol():
    run_taylor('ethanol', 'lj')
