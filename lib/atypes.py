#! /usr/bin/env python

from molmod.graphs import *
from molmod.molecular_graphs import *
from molmod.periodic import periodic as pt
import numpy as np

__all__ = ['assign_atypes', 'number_of_neighbors', 'number_of_X_neighbors']

filters = {
    1: [('Hc', CritAnd(HasAtomNumber(1), HasNeighborNumbers(6  ))),
        ('Ho', CritAnd(HasAtomNumber(1), HasNeighborNumbers(8  )))],
    8: [('Oa', CritAnd(HasAtomNumber(8), HasNeighborNumbers(6,1))),
        ('Ok', CritAnd(HasAtomNumber(8), HasNeighborNumbers(6  ))),
        ('Oe', CritAnd(HasAtomNumber(8), HasNeighborNumbers(6,6))),
        ('Ow', CritAnd(HasAtomNumber(8), HasNeighborNumbers(1,1)))],
}

def assign_atypes(system, level):
    if level=='low':
        return np.array([pt[number].symbol for number in system.numbers])
    elif level=='high':
        return np.array(['%s%i' %(pt[n].symbol, i) for i, n in enumerate(system.numbers)])
    elif level=='medium':
        graph = MolecularGraph(system.bonds, system.numbers)
        atypes = []
        for index, number in enumerate(system.numbers):
            if number in [6, 14]:
                atype = pt[number].symbol+str(number_of_neighbors(index,graph))
                numC = number_of_X_neighbors(index, 6, graph)
                if numC==1:
                    atype += '_c'
                elif numC>1:
                    atype += '_c%i' %numC
                numO = number_of_X_neighbors(index, 8, graph)
                if numO==1:
                    atype += '_o'
                elif numO>1:
                    atype += '_o%i' %numO
                atypes.append(atype)
            elif number in filters.keys():
                found = False
                for atype, filter in filters[number]:
                    if filter(index, graph):
                        atypes.append(atype)
                        found = True
                        break
                if not found:
                    atypes.append(pt[number].symbol)
            else:
                atypes.append(pt[number].symbol)
        return np.array(atypes)
    else:
        raise ValueError('Invalid level, recieved %s' %level)

def number_of_neighbors(index, graph):
    "Calculate the number of neighbors of the atom with index 'index'"
    result = 0
    for bond in graph.edges:
        if index in bond: result += 1
    return result

def number_of_X_neighbors(index, number, graph):
    "Calculate the number of neighbors of the atom with index 'index' that have atom number 'number'"
    result = 0
    for k,l in graph.edges:
        if (k==index and graph.numbers[l]==number) or (l==index and graph.numbers[k]==number):
            result += 1
    return result
