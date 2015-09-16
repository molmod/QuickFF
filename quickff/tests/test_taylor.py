from molmod.units import angstrom, kjmol
import numpy as np

from quickff.model import *
from common import get_system, get_scaled_pairs


def run_taylor(molecule, pot):
    tolE_abs = 5e-2
    system = get_system(molecule, ei_path='/wpart/he')
    scaled_pairs = get_scaled_pairs(system)
    coords0 = system.ref.coords
    if pot.lower() in ['coulpoint']:
        exact = CoulPointPot(system.charges, [1.0, 1.0, 1.0], scaled_pairs, coords0=coords0)
    elif pot.lower() in ['coulgauss']:
        system.read_uff_vdw()
        exact = CoulGaussPot(system.charges, 0.5*system.sigmas/2**(1.0/6.0), [1.0, 1.0, 1.0], scaled_pairs, coords0=coords0)
    elif pot.lower() in ['lj']:
        system.read_uff_vdw()
        exact = LennardJonesPot(system.sigmas, system.epsilons, [0.0,0.0,1.0], scaled_pairs, coords0=coords0)
    elif pot.lower() in ['mm3']:
        system.read_uff_vdw()
        exact = MM3BuckinghamPot(system.sigmas, system.epsilons, [0.0,0.0,1.0], scaled_pairs, coords0=coords0)
    else:
        raise ValueError('Invalic potential specification: %s' %pot)
    harm = HarmonicPot(exact.kind, coords0, 0.0, exact.calc_gradient(coords0), exact.calc_hessian(coords0))
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


def test_taylor_ei_point_water():
    run_taylor('water', 'coulpoint')

def test_taylor_ei_gauss_water():
    run_taylor('water', 'coulgauss')

def test_taylor_ei_point_ethanol():
    run_taylor('ethanol', 'coulpoint')

def test_taylor_ei_gauss_ethanol():
    run_taylor('ethanol', 'coulgauss')

def test_taylor_lj_ethanol():
    run_taylor('ethanol', 'lj')

def test_taylor_mm3_ethanol():
    run_taylor('ethanol', 'mm3')
