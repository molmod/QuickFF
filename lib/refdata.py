import numpy as np
from quickff.tools import global_translation, global_rotation

__all__ = ['ReferenceData']

class ReferenceData(object):
    '''
        A class to store the ab initio reference data from which the force field
        will be derived.
    '''
    def __init__(self, coords=None, grad=None, hess=None):
        '''
        **Optional Arguments**
        
        coords
            A (N,3) numpy array containing the cartesian coordinates of all
            atoms in equilibrium.
        
        grad
            A (N,3) numpy array containing the cartesian gradient of all
            atoms in equilibrium. Ideally, all elements of grad are zero.
        
        hess
            A (N,3,N,3) numpy array containing the cartesian Hessian in
            equilibrium.
        '''
        self.coords = coords
        self.grad = grad
        self.hess = hess
        if hess is not None:
            self.phess = self._get_phess()
        else:
            self.phess = None

    def _get_natoms(self):
        'Get the number of atoms in the reference data'
        return len(self.coords)

    _natoms = property(_get_natoms)

    def update(self, coords=None, grad=None, hess=None):
        '''
            Method to update one or more of the attributes. The meaning of the
            optional arguments is the same as with the constructor 
            :meth:`quickff.refdata.ReferenceData`
        '''
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
            self.phess = self._get_phess()

    def _check(self, natoms=None):
        'Internal consistency check'
        assert self.coords is not None
        assert self.grad is not None
        assert self.hess is not None
        assert len(self.coords.shape) == 2
        if natoms is not None:
            assert self.coords.shape[0] == natoms
        assert self.coords.shape[1] == 3
        assert self.grad.shape == self.coords.shape
        assert len(self.hess.shape) == 4
        assert self.hess.shape[0] == self.coords.shape[0]
        assert self.hess.shape[1] == 3
        assert self.hess.shape[2] == self.coords.shape[0]
        assert self.hess.shape[3] == 3

    def _get_phess(self):
        '''
            Constuct a hessian from which the translational and rotational
            degrees of freedom have been projected out.
        '''
        hess = self.hess.copy().reshape([3*self._natoms, 3*self._natoms])
        VTx, VTy, VTz = global_translation(self.coords)
        VRx, VRy, VRz = global_rotation(self.coords)
        U = np.linalg.svd(
            np.array([VTx, VTy, VTz, VRx, VRy, VRz]).transpose()
        )[0]
        proj = np.dot(U[:, 6:], U[:, 6:].T)
        PHP = np.dot(proj, np.dot(hess, proj))
        return PHP.reshape([self._natoms, 3, self._natoms, 3])
