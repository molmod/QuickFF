from molmod.units import angstrom, kjmol
from molmod.minimizer import check_delta
import numpy as np

from quickff.model import *

from common import get_water, get_ethanol

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

def get_hessian_check(pot):
    def fun(x, do_gradient=True):
        coords = x.reshape([-1, 3])
        gradient = pot.calc_gradient(coords)
        if do_gradient:
            hessian = pot.calc_hessian(coords)
            return gradient, hessian
        else:
            return gradient
    return fun

#test potentials for water
def test_coulomb_gradient_water():
    coords, numbers, fcharges = get_water()
    pot = CoulombPot(fcharges, [1.0, 1.0, 1.0],  [[], [], []])
    fun = get_gradient_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_coulomb_hessian_water():
    coords, numbers, fcharges = get_water()
    pot = CoulombPot(fcharges, [1.0, 1.0, 1.0],  [[], [], []])
    fun = get_hessian_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_harmonic_gradient_water():
    coords, numbers, fcharges = get_water()
    coul = CoulombPot(fcharges, [1.0, 1.0, 1.0],  [[], [], []])
    pot = HarmonicPot(coords, coul.calc_energy(coords), coul.calc_gradient(coords), coul.calc_hessian(coords))
    fun = get_gradient_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_harmonic_hessian_water():
    coords, numbers, fcharges = get_water()
    coul = CoulombPot(fcharges, [1.0, 1.0, 1.0],  [[], [], []])
    pot = HarmonicPot(coords, coul.calc_energy(coords), coul.calc_gradient(coords), coul.calc_hessian(coords))
    fun = get_hessian_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

#test potentials for ethanol
def test_coulomb_gradient_ethanol():
    coords, numbers, fcharges, sigmas, epsilons = get_ethanol()
    pot = CoulombPot(fcharges, [1.0, 1.0, 1.0],  [[], [], []])
    fun = get_gradient_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_coulomb_hessian_ethanol():
    coords, numbers, fcharges, sigmas, epsilons = get_ethanol()
    pot = CoulombPot(fcharges, [1.0, 1.0, 1.0],  [[], [], []])
    fun = get_hessian_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_lj_gradient_ethanol():
    coords, numbers, fcharges, sigmas, epsilons = get_ethanol()
    pot = LennartJonesPot(sigmas, epsilons, [1.0, 1.0, 1.0],  [[], [], []])
    fun = get_gradient_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_lj_hessian_ethanol():
    coords, numbers, fcharges, sigmas, epsilons = get_ethanol()
    pot = LennartJonesPot(sigmas, epsilons, [1.0, 1.0, 1.0],  [[], [], []])
    fun = get_hessian_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_harmonic_gradient_ethanol():
    coords, numbers, fcharges, sigmas, epsilons = get_ethanol()
    coul = CoulombPot(fcharges, [1.0, 1.0, 1.0],  [[], [], []])
    pot = HarmonicPot(coords, coul.calc_energy(coords), coul.calc_gradient(coords), coul.calc_hessian(coords))
    fun = get_gradient_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_harmonic_hessian_ethanol():
    coords, numbers, fcharges, sigmas, epsilons = get_ethanol()
    coul = CoulombPot(fcharges, [1.0, 1.0, 1.0],  [[], [], []])
    pot = HarmonicPot(coords, coul.calc_energy(coords), coul.calc_gradient(coords), coul.calc_hessian(coords))
    fun = get_hessian_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)
