from molmod.units import parse_unit

from quickff.system import System
from quickff.model import Model
from quickff.program import Program
from quickff.context import context

import numpy as np, os

def get_program(molecule, atypes_level, ei_scheme, ei_rule):
    moldir = context.get_fn('examples/%s' %molecule)
    system = System.from_files(
        [os.path.join(moldir, 'gaussian.fchk'), os.path.join(moldir, 'gaussian.fchk.h5')],
        ei_scheme=ei_scheme
    )
    system.guess_ffatypes(atypes_level)
    system.determine_ics_from_topology()
    print sorted(system.ics.keys())
    if ei_rule==-1:
        ei_pot_kind='Zero'
        ei_scales = [0.0,0.0,0.0]
    elif ei_rule==0:
        ei_pot_kind='Harm'
        ei_scales = [1.0,1.0,1.0]
    elif ei_rule==1:
        ei_pot_kind='Harm'
        ei_scales = [0.0,1.0,1.0]
    elif ei_rule==2:
        ei_pot_kind='Harm'
        ei_scales = [0.0,0.0,1.0]
    elif ei_rule==3:
        ei_pot_kind='Harm'
        ei_scales = [0.0,0.0,0.0]
    else:
        raise ValueError('Invalid ei-rule %i' %ei_rule)
    model = Model.from_system(system, ei_pot_kind=ei_pot_kind, ei_scales=ei_scales)
    model.val.determine_dihedral_potentials(system, verbose=False)
    program = Program(system, model)
    return program
    
def get_ref_ffs(molecule, atypes_level, ei_scheme, ei_rule):
    moldir = context.get_fn('examples/%s' %molecule)
    fn = os.path.join(moldir, '%s_%s_%s.qff' %(atypes_level, ei_rule, ei_scheme) )
    ffs = read_ff(fn)
    return ffs

def read_ff(fn):
    ffs = {'pre': {}, 'post': {}}
    f = open(fn, 'r')
    section = 'skip'
    for line in f.readlines():
        words = line.split()
        if len(words)==0:
            continue
        if line.startswith('Estimating all pars'):
            section = 'pre'
            continue
        elif line.startswith('Refining force constants'):
            section = 'post'
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
