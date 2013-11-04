from molmod.units import parse_unit, angstrom, kjmol
import numpy as np, os

from quickff.system import System
from quickff.context import context

__all__ = [
    'get_system', 'get_water', 'get_ethanol', 'translate_rule',
    'get_scaled_pairs', 'read_qff_out',
]

def get_system(molecule, atypes_level='high', ei_scheme='he', vdw_from='uff', verbose=False):
    moldir = context.get_fn('examples/%s' %molecule)
    system = System.from_files(
        [os.path.join(moldir, 'gaussian.fchk'), os.path.join(moldir, 'gaussian.fchk.h5')],
        ei_scheme=ei_scheme
    )
    system.guess_ffatypes(atypes_level)
    if vdw_from.lower() in ['uff']:
        system.read_uff_vdw()
    else:
        raise ValueError('Invalid value for vdw_from: %s' %vdw_from)
    system.determine_ics_from_topology()
    if verbose:
        system.print_atom_info()
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
    
def read_qff_out(fn):
    ffs = {'pt': {}, 'refined': {}}
    f = open(fn, 'r')
    section = 'skip'
    for line in f.readlines():
        words = line.split()
        if len(words)==0:
            continue
        if line.startswith('Estimating all pars'):
            section = 'pt'
            continue
        elif line.startswith('Refining force constants'):
            section = 'refined'
            continue
        elif line.startswith('~~~~~~~~~~~~~~~~~~~~~~~~'):
            section = 'skip'
        if section != 'skip':
            icname = words[0]
            k = float(words[3])*parse_unit(words[6])
            q0 = None
            if len(words)>7:
                q0 = float(words[9])*parse_unit(words[12])
            ffs[section][icname] = [k, q0]
    f.close()
    return ffs
