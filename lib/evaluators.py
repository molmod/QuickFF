#! /usr/bin/env python

__all__=['eval_ic', 'eval_energy']

def eval_ic(ic):
    def eval_ic(model, coords):
        return ic.value(coords)
    return eval_ic


def eval_energy(part_name):
    def evaluator(model, coords):
        part_model = getattr(model, part_name)
        return part_model.calc_energy(coords)
    return evaluator
