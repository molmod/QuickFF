#! /usr/bin/env python

import sys
sys.path.append('/home/louis/Documents/Doctoraat/Code/mil53/model/metals/covalent/est/lib')

from model import *
from ic import *
from system import *

__all__=['ic_evaluator', 'energy_evaluator']

def ic_evaluator(name=None, i=-1, ic=None):
    def eval_ic(system, coords):
        if isinstance(ic, IC):
            match = ic
        else:
            assert name is not None and i>-1
            match = system.ics[name][i]
        return match.value(coords)
    return eval_ic


def energy_evaluator(model_name):
    def eval_energy(system, coords):
        model = getattr(system, model_name)
        return model.get_energy(coords)
    return eval_energy
