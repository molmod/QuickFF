#! /usr/bin/env python

import sys
sys.path.append('/home/louis/Documents/Doctoraat/Code/mil53/model/metals/covalent/est/lib')

from model import *
from ic import *
from system import *

__all__=['ic_evaluator', 'energy_evaluator']

def ic_evaluator(name, i):
    def eval_ic(system, coords):
        ic = system.ics[name][i]
        return ic.value(coords)
    return eval_ic


def energy_evaluator(model_name):
    def eval_energy(system, coords):
        model = getattr(system, model_name)
        return model.get_energy(coords)
    return eval_energy
