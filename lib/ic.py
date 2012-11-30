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
        rs = coords[self.indexes]
        q = self.icf(rs, deriv=0)[0]
        return q
    
    def grad(self, coords):
        rs = coords[self.indexes]
        q, q_deriv = self.icf(rs, deriv=1)
        grad = np.zeros(coords.shape, float)
        grad[self.indexes] = q_deriv
        grad = grad.reshape(len(coords)*3)
        return grad
    
    def hess(self, coords):
        N = coords.shape[0]
        rs = coords[self.indexes]
        q, q_deriv, q_hess = self.icf(rs, deriv=2)
        hess = np.zeros([N, 3, N, 3], float)
        for i1, index1 in enumerate(self.indexes):
            for i2, index2 in enumerate(self.indexes):
                hess[index1, :, index2, :] = q_hess[i1, :, i2, :]
        hess = hess.reshape([3*N, 3*N])
        return hess

    def test(self, coords, epsilon=1e-4, ntests=None, threshold=1e-5):
        N = len(coords)
        if ntests is None: ntests = 3*N
        oom = np.linalg.norm(coords)
        hess = self.hess(coords)
        for i in xrange(ntests):
            v = np.random.normal(size=coords.shape)
            v /= np.linalg.norm(v)
            gp = self.grad(coords+epsilon*oom*v).reshape(3*N)
            gm = self.grad(coords-epsilon*oom*v).reshape(3*N)
            num = (gp-gm)/(2.0*epsilon)
            ana = oom*np.dot(hess, v.reshape(3*N))
            error = np.linalg.norm(num-ana)/np.linalg.norm(num)
            if error>threshold:
                print '    IC TEST : test nr. %i for %s FAILED: epsilon=%.6f  norm(num)=%.6e  error=%.6e' %(i,self.name,epsilon,np.linalg.norm(num),error)
