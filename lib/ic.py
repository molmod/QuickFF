#! /usr/bin/env python

import numpy as np

__all__ = ['IC']

class IC(object):
    def __init__(self, name, indexes, icf, qunit='au', kunit='kjmol/au**2'):
        self.name = name
        self.indexes = indexes
        self.icf = icf
        self.qunit = qunit
        self.kunit = kunit

    def value(self, coords):
        rs = coords[self.indexes]
        q = self.icf(rs, deriv=0)[0]
        return q

    def grad(self, coords):
        Natoms = len(coords)
        rs = coords[self.indexes]
        q, q_deriv = self.icf(rs, deriv=1)
        grad = np.zeros(coords.shape, float)
        grad[self.indexes] = q_deriv
        grad = grad.reshape([3*Natoms])
        return grad

    def hess(self, coords):
        Natoms = len(coords)
        rs = coords[self.indexes]
        q, q_deriv, q_hess = self.icf(rs, deriv=2)
        hess = np.zeros([Natoms, 3, Natoms, 3], float)
        for i1, index1 in enumerate(self.indexes):
            for i2, index2 in enumerate(self.indexes):
                hess[index1, :, index2, :] = q_hess[i1, :, i2, :]
        hess = hess.reshape([3*Natoms, 3*Natoms])
        return hess

    def _get_pairs(self):
        '''
            A method to determine all pairs of atoms in this ic. This is
            usefull later for determining which atom pairs are connected
            by internal coordinates and hence which AI Hessian elements
            are relevant.
        '''
        result = []
        for i in xrange(self.ic.indexes):
            for j in xrange(i+1):
                result.append([i,j])
        return result

    bpairs = property(_get_pairs)
