from molmod.units import parse_unit, angstrom, kjmol
import numpy as np, os

from quickff.system import System
from quickff.context import context

__all__ = ['get_water', 'get_ethanol', 'translate_rule', 'get_scaled_pairs']

def get_system(molecule, atypes_level='high', ei_path=None):
    moldir = context.get_fn('examples/%s' %molecule)
    if ei_scheme is not None:
        system = System.from_files(
            [os.path.join(moldir, 'gaussian.fchk'), os.path.join(moldir, 'gaussian_wpart.h5')],
            ei_path=ei_path
        )
    else:
        system = System.from_files([os.path.join(moldir, 'gaussian.fchk')])
    system.guess_ffatypes(atypes_level)
    system.determine_ics_from_topology()
    return system

def get_water():
    coords = np.array([
        [ 0.000000000,  0.763382315, -0.468300621],
        [ 0.000000000, -0.000000000,  0.117075156],
        [-0.000000000, -0.763382315, -0.468300621],
    ])*angstrom
    numbers = np.array([1,8,1])
    fcharges = np.array([0.5, -1, 0.5]) #'formal' charges
    return coords, numbers, fcharges

def get_ethanol():
    coords = np.array([
        [ 1.207471, -0.223281,  0.000000],
        [-0.075681,  0.553519, -0.000000],
        [ 1.250572, -0.854108,  0.888799],
        [ 1.250572, -0.854109, -0.888799],
        [ 2.055418,  0.462705, -0.000000],
        [-0.147118,  1.198973, -0.883292],
        [-1.097876, -0.419344,  0.000000],
        [-0.147118,  1.198973,  0.883292],
        [-1.917880,  0.080549,  0.000000],
    ])*angstrom
    numbers = np.array([6, 6, 1, 1, 1, 1, 8, 1, 1])
    fcharges = np.array([-0.3, 0.3, 0.1, 0.1, 0.1, 0.1, -1.0, 0.1, 0.5])                #'formal' charges
    sigmas = np.array([3.431,3.431,2.571,2.571,2.571,2.571,3.118,2.571,2.571])*angstrom #UFF
    epsilons = np.array([0.439,0.439,0.184,0.184,0.184,0.184,0.251,0.184,0.184])*kjmol  #UFF
    return coords, numbers, fcharges, sigmas, epsilons

def translate_rule(rule):
    if rule==-1:
        pot_kind='Zero'
        scales = [0.0,0.0,0.0]
    elif rule==0:
        pot_kind='Harm'
        scales = [1.0,1.0,1.0]
    elif rule==1:
        pot_kind='Harm'
        scales = [0.0,1.0,1.0]
    elif rule==2:
        pot_kind='Harm'
        scales = [0.0,0.0,1.0]
    elif rule==3:
        pot_kind='Harm'
        scales = [0.0,0.0,0.0]
    else:
        raise ValueError('Invalid rule %i' %rule)
    return scales, pot_kind

def get_scaled_pairs(system):
    'Generate list of atom pairs subject to scaling of non-bonding interactions'
    scaled_pairs = [[],[],[]]
    for bond in system.bonds:
        scaled_pairs[0].append([bond[0], bond[1]])
    for bend in system.bends:
        scaled_pairs[1].append([bend[0], bend[2]])
    for dihed in system.diheds:
        scaled_pairs[2].append([dihed[0], dihed[3]])
    return scaled_pairs
