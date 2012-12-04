#! /usr/bin/env python

from molmod.periodic import periodic as pt
import numpy as np

__all__ = ['IC']

class IC(object):
    def __init__(self, indexes, icf, name=None):
        self.indexes = indexes
        self.icf = icf
        self.name = name
    
    def value(self, coords):
        Natoms = len(coords)
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

    def test(self, coords, epsilon=1e-4, ntests=None, threshold=1e-5):
        Natoms = len(coords)
        if ntests is None: ntests = 3*Natoms
        oom = np.linalg.norm(coords)
        hess = self.hess(coords)
        for i in xrange(ntests):
            v = np.random.normal(size=coords.shape)
            v /= np.linalg.norm(v)
            gp = self.grad(coords+epsilon*oom*v)
            gm = self.grad(coords-epsilon*oom*v)
            num = (gp-gm)/(2.0*epsilon)
            ana = oom*np.dot( hess , v.reshape([3*Natoms]) )
            error = np.linalg.norm(num-ana)/np.linalg.norm(num)
            if error>threshold:
                print '    IC TEST : test nr. %i for %s FAILED: epsilon=%.6f  norm(num)=%.6e  error=%.6e' %(i,self.name,epsilon,np.linalg.norm(num),error)
