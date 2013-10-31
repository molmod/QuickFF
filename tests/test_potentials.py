from molmod.units import angstrom, kjmol
from molmod.minimizer import check_delta
import numpy as np

from quickff.model import *
from quickff.perturbation import RelaxedGeoPertTheory
from quickff.uff import get_uff_sigmas_epsilons

from common import get_water_coords, get_system

def get_gradient_check(pot):
    def fun(x, do_gradient=True):
        coords = x.reshape([-1, 3])
        energy = pot.calc_energy(coords)
        if do_gradient:
            gradient = pot.calc_gradient(coords)
            return energy, gradient
        else:
            return energy
    return fun

def get_hessian_check(pot, i):
    def fun(x, do_gradient=True):
        coords = x.reshape([-1, 3])
        gradient = pot.calc_gradient(coords)[i]
        if do_gradient:
            hessian = pot.calc_hessian(coords)[i,:]
            return gradient, hessian
        else:
            return gradient
    return fun

def run_ei_taylor(molecule):
    tolE = 1e-1*kjmol
    system = get_system(molecule, 'high', 'he')
    coords0 = system.ref.coords
    coul = CoulombPot(system.charges, [1.0, 1.0, 1.0],  [[], [], []], coords0=coords0.copy())
    model = Model.from_system(system, ei_scales=[1.0, 1.0, 1.0], ei_pot_kind='harm', vdw_pot_kind='zero')
    model.val.determine_dihedral_potentials(system, verbose=False)
    pt = RelaxedGeoPertTheory(system, model)
    print ''
    for icname, ics in system.ics.iteritems():
        for ic in ics:
            trajectory = pt.generate(ic)
            for i, dx in enumerate(trajectory):
                exact = coul.calc_energy(coords0+dx)
                appro = model.ei.calc_energy(coords0+dx)
                print '%s step %i: exact = %.9f kjmol    appro = %.9f kjmol' %(
                    icname, i, exact/kjmol, appro/kjmol
                )
                assert abs(exact-appro)<tolE

def test_pot_coulomb_gradient_water():
    coords0 = get_water_coords()
    pot = CoulombPot([1.0, -2.0, 1.0], [1.0, 1.0, 1.0],  [[], [], []])
    fun = get_gradient_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 9])
    check_delta(fun, coords0.ravel(), dxs)

def test_pot_coulomb_hessian_water():
    coords0 = get_water_coords()
    pot = CoulombPot([1.0, -2.0, 1.0], [1.0, 1.0, 1.0],  [[], [], []])
    for i in xrange(3*len(coords0)):
        fun = get_hessian_check(pot, i)
        dxs = np.random.normal(0.0, 1e-4, [100, 9])
        check_delta(fun, coords0.ravel(), dxs)

def test_pot_harmonic_gradient_water():
    coords0 = get_water_coords()
    coul = CoulombPot([1.0, -2.0, 1.0], [1.0, 1.0, 1.0],  [[], [], []])
    pot = HarmonicPot(coords0, coul.calc_energy(coords0), coul.calc_gradient(coords0), coul.calc_hessian(coords0))
    fun = get_gradient_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 9])
    check_delta(fun, coords0.ravel(), dxs)

def test_pot_harmonic_hessian_water():
    coords0 = get_water_coords()
    coul = CoulombPot([1.0, -2.0, 1.0], [1.0, 1.0, 1.0],  [[], [], []])
    pot = HarmonicPot(coords0, coul.calc_energy(coords0), coul.calc_gradient(coords0), coul.calc_hessian(coords0))
    for i in xrange(3*len(coords0)):
        fun = get_hessian_check(pot, i)
        dxs = np.random.normal(0.0, 1e-4, [100, 9])
        check_delta(fun, coords0.ravel(), dxs)

def test_ei_taylor_water():
    run_ei_taylor('water')

def test_ei_taylor_ethanol():
    run_ei_taylor('ethanol')
