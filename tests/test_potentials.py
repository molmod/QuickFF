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
def test_coulpoint_gradient_water():
    coords, numbers, fcharges, fei_sigmas = get_water()
    pot = CoulPointPot(fcharges, [1.0, 1.0, 1.0],  [[], [], []])
    fun = get_gradient_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_coulpoint_hessian_water():
    coords, numbers, fcharges, fei_sigmas = get_water()
    pot = CoulPointPot(fcharges, [1.0, 1.0, 1.0],  [[], [], []])
    fun = get_hessian_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_coulgauss_gradient_water():
    coords, numbers, fcharges, fei_sigmas = get_water()
    pot = CoulGaussPot(fcharges, fei_sigmas, [1.0, 1.0, 1.0],  [[], [], []])
    fun = get_gradient_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_coulgauss_hessian_water():
    coords, numbers, fcharges, fei_sigmas = get_water()
    pot = CoulGaussPot(fcharges, fei_sigmas, [1.0, 1.0, 1.0],  [[], [], []])
    fun = get_hessian_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_harmpoint_gradient_water():
    coords, numbers, fcharges, fei_sigmas = get_water()
    coul = CoulPointPot(fcharges, [1.0, 1.0, 1.0],  [[], [], []])
    pot = HarmonicPot(coords, coul.calc_energy(coords), coul.calc_gradient(coords), coul.calc_hessian(coords))
    fun = get_gradient_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_harmpoint_hessian_water():
    coords, numbers, fcharges, fei_sigmas = get_water()
    coul = CoulPointPot(fcharges, [1.0, 1.0, 1.0],  [[], [], []])
    pot = HarmonicPot(coords, coul.calc_energy(coords), coul.calc_gradient(coords), coul.calc_hessian(coords))
    fun = get_hessian_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

#test potentials for ethanol
def test_coulpoint_gradient_ethanol():
    coords, numbers, fcharges, fei_sigmas, epsilons, vdw_sigmas = get_ethanol()
    pot = CoulPointPot(fcharges, [1.0, 1.0, 1.0],  [[], [], []])
    fun = get_gradient_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_coulpoint_hessian_ethanol():
    coords, numbers, fcharges, fei_sigmas, epsilons, vdw_sigmas = get_ethanol()
    pot = CoulPointPot(fcharges, [1.0, 1.0, 1.0],  [[], [], []])
    fun = get_hessian_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_coulgauss_gradient_ethanol():
    coords, numbers, fcharges, fei_sigmas, epsilons, vdw_sigmas = get_ethanol()
    pot = CoulGaussPot(fcharges, fei_sigmas, [1.0, 1.0, 1.0],  [[], [], []])
    fun = get_gradient_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_coulgauss_hessian_ethanol():
    coords, numbers, fcharges, fei_sigmas, epsilons, vdw_sigmas = get_ethanol()
    pot = CoulGaussPot(fcharges, fei_sigmas, [1.0, 1.0, 1.0],  [[], [], []])
    fun = get_hessian_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_lj_gradient_ethanol():
    coords, numbers, fcharges, fei_sigmas, epsilons, vdw_sigmas = get_ethanol()
    pot = LennardJonesPot(vdw_sigmas, epsilons, [1.0, 1.0, 1.0],  [[], [], []])
    fun = get_gradient_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_lj_hessian_ethanol():
    coords, numbers, fcharges, fei_sigmas, epsilons, vdw_sigmas = get_ethanol()
    pot = LennardJonesPot(vdw_sigmas, epsilons, [1.0, 1.0, 1.0],  [[], [], []])
    fun = get_hessian_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_mm3_gradient_ethanol():
    coords, numbers, fcharges, fei_sigmas, epsilons, vdw_sigmas = get_ethanol()
    pot = MM3BuckinghamPot(vdw_sigmas, epsilons, [1.0, 1.0, 1.0],  [[], [], []])
    fun = get_gradient_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_mm3_hessian_ethanol():
    coords, numbers, fcharges, fei_sigmas, epsilons, vdw_sigmas = get_ethanol()
    pot = MM3BuckinghamPot(vdw_sigmas, epsilons, [1.0, 1.0, 1.0],  [[], [], []])
    fun = get_hessian_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_harmpoint_gradient_ethanol():
    coords, numbers, fcharges, fei_sigmas, epsilons, vdw_sigmas = get_ethanol()
    coul = CoulPointPot(fcharges, [1.0, 1.0, 1.0],  [[], [], []])
    pot = HarmonicPot(coords, coul.calc_energy(coords), coul.calc_gradient(coords), coul.calc_hessian(coords))
    fun = get_gradient_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_harmpoint_hessian_ethanol():
    coords, numbers, fcharges, fei_sigmas, epsilons, vdw_sigmas = get_ethanol()
    coul = CoulPointPot(fcharges, [1.0, 1.0, 1.0],  [[], [], []])
    pot = HarmonicPot(coords, coul.calc_energy(coords), coul.calc_gradient(coords), coul.calc_hessian(coords))
    fun = get_hessian_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)
