from molmod.units import angstrom, kjmol
from molmod.minimizer import check_delta
import numpy as np
from nose.tools import assert_raises
from yaff import System, NeighborList, Scalings, PairPotEI, ForcePartPair, \
    ForceField, ForcePartValence, PairPotLJ
from yaff import log

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
    coords, numbers, fcharges, fradii = get_water()
    pot = CoulPointPot(fcharges, [1.0, 1.0, 1.0],  [[], [], []])
    fun = get_gradient_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_coulpoint_hessian_water():
    coords, numbers, fcharges, fradii = get_water()
    pot = CoulPointPot(fcharges, [1.0, 1.0, 1.0],  [[], [], []])
    fun = get_hessian_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_coulgauss_gradient_water():
    coords, numbers, fcharges, fradii = get_water()
    pot = CoulGaussPot(fcharges, fradii, [1.0, 1.0, 1.0],  [[], [], []])
    fun = get_gradient_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_coulgauss_hessian_water():
    coords, numbers, fcharges, fradii = get_water()
    pot = CoulGaussPot(fcharges, fradii, [1.0, 1.0, 1.0],  [[], [], []])
    fun = get_hessian_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_harmpoint_gradient_water():
    coords, numbers, fcharges, fradii = get_water()
    coul = CoulPointPot(fcharges, [1.0, 1.0, 1.0],  [[], [], []])
    pot = HarmonicPot(coul.kind, coords, coul.calc_energy(coords), coul.calc_gradient(coords), coul.calc_hessian(coords))
    fun = get_gradient_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_harmpoint_hessian_water():
    coords, numbers, fcharges, fradii = get_water()
    coul = CoulPointPot(fcharges, [1.0, 1.0, 1.0],  [[], [], []])
    pot = HarmonicPot(coul.kind, coords, coul.calc_energy(coords), coul.calc_gradient(coords), coul.calc_hessian(coords))
    fun = get_hessian_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_nbyaff_coulomb_water():
    log.set_level(log.silent)
    coords, numbers, fcharges, fradii = get_water()
    pot = CoulPointPot(fcharges, [1.0, 1.0, 1.0],  [[], [], []])
    # Calculate energy, gradient and hessian using QuickFF
    E_QFF = pot.calc_energy(coords)
    G_QFF = pot.calc_gradient(coords)
    H_QFF = pot.calc_hessian(coords)
    # Make a Yaff system
    system = System(numbers,coords,charges=fcharges)
    nlist = NeighborList(system)
    scalings = Scalings(system, 1.0,1.0,1.0)
    # Make a Yaff pair potential for electrostatic interactions
    pair_pot = PairPotEI(system.charges,0.0,100.0*angstrom)
    part_pair = ForcePartPair(system, nlist, scalings, pair_pot)
    ff = ForceField(system, [part_pair], nlist)
    pot2 = NonbondedYaffPot(ff)
    # Calculate energy, gradient and hessian using Yaff (called from QuickFF)
    E_YAFF = pot2.calc_energy(coords)
    G_YAFF = pot2.calc_gradient(coords)
    H_YAFF = pot2.calc_hessian(coords)
    assert np.abs(E_YAFF-E_QFF) < 1e-8
    assert np.all( np.abs(G_YAFF-G_QFF)) < 1e-8
    assert np.all( np.abs(H_QFF-H_YAFF) ) < 1e-8


def test_nbyaff_valence_water():

    coords, numbers, fcharges, fradii = get_water()
    system = System(numbers, coords)
    part = ForcePartValence(system)
    ff = ForceField(system,[part])
    # We do not allow to include a covalent part in the NonbondedYaff,
    # because this is probably a mistake made by the user
    with assert_raises(UserWarning):
        nbyaff = NonbondedYaffPot(ff)
    

#test potentials for ethanol
def test_coulpoint_gradient_ethanol():
    coords, numbers, fcharges, fradii, epsilons, sigmas = get_ethanol()
    pot = CoulPointPot(fcharges, [1.0, 1.0, 1.0],  [[], [], []])
    fun = get_gradient_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_coulpoint_hessian_ethanol():
    coords, numbers, fcharges, fradii, epsilons, sigmas = get_ethanol()
    pot = CoulPointPot(fcharges, [1.0, 1.0, 1.0],  [[], [], []])
    fun = get_hessian_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_coulgauss_gradient_ethanol():
    coords, numbers, fcharges, fradii, epsilons, sigmas = get_ethanol()
    pot = CoulGaussPot(fcharges, fradii, [1.0, 1.0, 1.0],  [[], [], []])
    fun = get_gradient_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_coulgauss_hessian_ethanol():
    coords, numbers, fcharges, fradii, epsilons, sigmas = get_ethanol()
    pot = CoulGaussPot(fcharges, fradii, [1.0, 1.0, 1.0],  [[], [], []])
    fun = get_hessian_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_lj_gradient_ethanol():
    coords, numbers, fcharges, fradii, epsilons, sigmas = get_ethanol()
    pot = LennardJonesPot(sigmas, epsilons, [1.0, 1.0, 1.0],  [[], [], []])
    fun = get_gradient_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_lj_hessian_ethanol():
    coords, numbers, fcharges, fradii, epsilons, sigmas = get_ethanol()
    pot = LennardJonesPot(sigmas, epsilons, [1.0, 1.0, 1.0],  [[], [], []])
    fun = get_hessian_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_mm3_gradient_ethanol():
    coords, numbers, fcharges, fradii, epsilons, sigmas = get_ethanol()
    pot = MM3BuckinghamPot(sigmas, epsilons, [1.0, 1.0, 1.0],  [[], [], []])
    fun = get_gradient_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_mm3_hessian_ethanol():
    coords, numbers, fcharges, fradii, epsilons, sigmas = get_ethanol()
    pot = MM3BuckinghamPot(sigmas, epsilons, [1.0, 1.0, 1.0],  [[], [], []])
    fun = get_hessian_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_harmpoint_gradient_ethanol():
    coords, numbers, fcharges, fradii, epsilons, sigmas = get_ethanol()
    coul = CoulPointPot(fcharges, [1.0, 1.0, 1.0],  [[], [], []])
    pot = HarmonicPot(coul.kind, coords, coul.calc_energy(coords), coul.calc_gradient(coords), coul.calc_hessian(coords))
    fun = get_gradient_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)

def test_harmpoint_hessian_ethanol():
    coords, numbers, fcharges, fradii, epsilons, sigmas = get_ethanol()
    coul = CoulPointPot(fcharges, [1.0, 1.0, 1.0],  [[], [], []])
    pot = HarmonicPot(coul.kind, coords, coul.calc_energy(coords), coul.calc_gradient(coords), coul.calc_hessian(coords))
    fun = get_hessian_check(pot)
    dxs = np.random.normal(0.0, 1e-4, [100, 3*len(numbers)])
    check_delta(fun, coords.ravel(), dxs)


def test_nbyaff_lj_ethanol():
    log.set_level(log.silent)
    coords, numbers, fcharges, fradii, epsilons, sigmas = get_ethanol()
    pot = LennardJonesPot(sigmas, epsilons, [1.0, 1.0, 1.0],  [[], [], []])
    # Calculate energy, gradient and hessian using QuickFF
    E_QFF = pot.calc_energy(coords)
    G_QFF = pot.calc_gradient(coords)
    H_QFF = pot.calc_hessian(coords)
    # Make a Yaff system
    system = System(numbers,coords,charges=fcharges)
    nlist = NeighborList(system)
    scalings = Scalings(system, 1.0,1.0,1.0)
    # Make a Yaff pair potential for electrostatic interactions
    pair_pot = PairPotLJ(sigmas, epsilons, 100.0*angstrom)
    part_pair = ForcePartPair(system, nlist, scalings, pair_pot)
    ff = ForceField(system, [part_pair], nlist)
    pot2 = NonbondedYaffPot(ff)
    # Calculate energy, gradient and hessian using Yaff (called from QuickFF)
    E_YAFF = pot2.calc_energy(coords)
    G_YAFF = pot2.calc_gradient(coords)
    H_YAFF = pot2.calc_hessian(coords)
    assert np.abs(E_YAFF-E_QFF) < 1e-8
    assert np.all( np.abs(G_YAFF-G_QFF)) < 1e-8
    assert np.all( np.abs(H_QFF-H_YAFF) ) < 1e-8
