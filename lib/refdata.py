import numpy as np
from quickff.tools import global_translation, global_rotation

__all__ = ['ReferenceData']

class ReferenceData(object):
    'A class to store the ab initio reference data'
    def __init__(self):
        self.coords = None
        self.grad = None
        self.hess = None

    def _get_natoms(self):
        'Get the number of atoms in the reference data'
        return len(self.coords)

    natoms = property(_get_natoms)

    def update(self, coords=None, grad=None, hess=None):
        'Method to update one of the attributes'
        if coords is not None:
            if self.coords is not None:
                assert self.coords.shape == coords.shape
            self.coords = coords
        if grad is not None:
            if self.grad is not None:
                assert self.grad.shape == grad.shape
            self.grad = grad
        if hess is not None:
            if self.hess is not None:
                assert self.hess.shape == hess.shape
            self.hess = hess

    def check(self):
        'Internal consistency check'
        assert self.coords is not None
        assert self.grad is not None
        assert self.hess is not None
        assert len(self.coords.shape) == 2
        assert self.coords.shape[1] == 3
        assert self.coords.shape == self.grad.shape
        assert self.hess.shape[0] == self.coords.shape[0]
        assert self.hess.shape[1] == 3
        assert self.hess.shape[2] == self.coords.shape[0]
        assert self.hess.shape[3] == 3

    def _get_phess(self):
        '''
            Constuct a hessian from which the translational and rotational
            degrees of freedom have been projected out.
        '''
        hess = self.hess.copy().reshape([3*self.natoms, 3*self.natoms])
        VTx, VTy, VTz = global_translation(self.coords)
        VRx, VRy, VRz = global_rotation(self.coords)
        U = np.linalg.svd(
            np.array([VTx, VTy, VTz, VRx, VRy, VRz]).transpose()
        )[0]
        proj = np.dot(U[:, 6:], U[:, 6:].T)
        PHP = np.dot(proj, np.dot(hess, proj))
        return PHP.reshape([self.natoms, 3, self.natoms, 3])

    phess = property(_get_phess)
