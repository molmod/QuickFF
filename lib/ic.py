#! /usr/bin/env python

import numpy as np

__all__ = ['IC']

class IC(object):
    '''
        A class functioning as a wrapper around the MolMod ic functions.
        This class should be able to calculate the value of an internal
        coordinate and its first and second order derivatives toward the
        cartesian coordinates.
    '''
    def __init__(self, name, indexes, icf, qunit='au', kunit='kjmol/au**2'):
        self.name = name
        self.indexes = indexes
        self.icf = icf
        self.qunit = qunit
        self.kunit = kunit

    def value(self, coords):
        'Calculate the value of the internal coordinate'
        _rs = coords[self.indexes]
        return self.icf(_rs, deriv=0)[0]

    def grad(self, coords):
        'Calculate the first order derivative of the internal coordinate'
        _rs = coords[self.indexes]
        _q_deriv = self.icf(_rs, deriv=1)[1]
        grad = np.zeros(coords.shape, float)
        grad[self.indexes] = _q_deriv
        grad = grad.reshape([3*len(coords)])
        return grad

    def hess(self, coords):
        'Calculate the second order derivative of the internal coordinate'
        Natoms = len(coords)
        _rs = coords[self.indexes]
        _q_hess = self.icf(_rs, deriv=2)[2]
        hess = np.zeros([Natoms, 3, Natoms, 3], float)
        for i_1, index1 in enumerate(self.indexes):
            for i_2, index2 in enumerate(self.indexes):
                hess[index1, :, index2, :] = _q_hess[i_1, :, i_2, :]
        hess = hess.reshape([3*Natoms, 3*Natoms])
        return hess
