#! /usr/bin/env python

import numpy as np
from tools import global_translation, global_rotation

__all__ = ['ReferenceData']

class ReferenceData(object):
    def __init__(self):
        self.coords = None
        self.grad = None
        self.hess = None

    def _get_natoms(self):
        return len(self.coords)

    natoms = property(_get_natoms)

    def update(self, coords=None, grad=None, hess=None):
        if coords is not None:
            if self.coords is not None:
                assert self.coords.shape==coords.shape
            self.coords = coords
        if grad is not None:
            if self.grad is not None:
                assert self.grad.shape==grad.shape
            self.grad = grad
        if hess is not None:
            if self.hess is not None:
                assert self.hess.shape==hess.shape
            self.hess = hess

    def check(self):
        assert self.coords is not None
        assert self.grad is not None
        assert self.hess is not None
        assert len(self.coords.shape)==2
        assert self.coords.shape[1]==3
        assert self.coords.shape==self.grad.shape
        assert self.hess.shape==(self.coords.shape[0],3,self.coords.shape[0],3)

    def _get_phess(self):
        '''
            Constuct a hessian from which the translational and rotational
            degrees of freedom have been projected out.
        '''
        hess = self.hess.copy().reshape([3*self.natoms, 3*self.natoms])
        VTx, VTy, VTz = global_translation(self.coords)
        VRx, VRy, VRz = global_rotation(self.coords)
        U, S, Vt = np.linalg.svd( np.array([VTx, VTy, VTz, VRx, VRy, VRz]).transpose() )
        P = np.dot(U[:,6:], U[:,6:].T)
        return np.dot(P, np.dot(hess, P)).reshape([self.natoms, 3, self.natoms, 3])

    phess = property(_get_phess)
