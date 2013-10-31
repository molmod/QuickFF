from molmod.units import parse_unit, angstrom

from quickff.system import System
from quickff.context import context

import numpy as np, os

def get_water_coords():
    coords0 = np.array([
        [ 0.000000000,  0.763382315, -0.468300621],
        [ 0.000000000, -0.000000000,  0.117075156],
        [-0.000000000, -0.763382315, -0.468300621],
    ])*angstrom
    return coords0

def get_system(molecule, atypes_level, ei_scheme):
    moldir = context.get_fn('examples/%s' %molecule)
    system = System.from_files(
        [os.path.join(moldir, 'gaussian.fchk'), os.path.join(moldir, 'gaussian.fchk.h5')],
        ei_scheme=ei_scheme
    )
    system.guess_ffatypes(atypes_level)
    system.determine_ics_from_topology()
    return system
    
def get_reference_ff(molecule, atypes_level, ei_scheme, ei_rule):
    moldir = context.get_fn('examples/%s' %molecule)
    fn = os.path.join(moldir, '%s_%s_%s.qff' %(atypes_level, ei_rule, ei_scheme) )
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
