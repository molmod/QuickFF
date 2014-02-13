__all__ = ['eval_ic', 'eval_energy']

def eval_ic(ic):
    '''
        Return an evaluator for the value of the given internal coordinate, i.e.
        return a method that calculates the value of the given ic for a certain
        atomic configuration.
    
        **Arguments**
        
        ic
            an instance of the class :class:`quickff.ic.IC`
    '''
    def evaluator(model, coords):
        return ic.value(coords)
    return evaluator


def eval_energy(part_name):
    '''
        Return an evaluator for the value of the `part_name` energy, i.e.
        return a method that calculates the `part_name` energy for a certain
        atomic configuration.
        
        **Arguments**
        
        part_name
            a string defining an attribute of a model instance that is needed
            to evaluate the correct part of the energy.
    '''
    def evaluator(model, coords):
        part_model = getattr(model, part_name)
        return part_model.calc_energy(coords)
    return evaluator
