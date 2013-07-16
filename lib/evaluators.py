#! /usr/bin/env python

__all__ = ['eval_ic', 'eval_energy']

def eval_ic(ic):
    "Evaluator for the value of internal coordinate <<ic>>"
    def evaluator(model, coords):
        return ic.value(coords)
    return evaluator


def eval_energy(part_name):
    "Evaluator for the value of the <<part_name>> energy"
    def evaluator(model, coords):
        part_model = getattr(model, part_name)
        return part_model.calc_energy(coords)
    return evaluator
