#! /usr/bin/env python

from molmod.units import *

import numpy as np

from tools import get_atoms_within_3bonds, matrix_squared_sum
from scipy.optimize import minimize

__all__ = ['HessianFCCost']

class HessianFCCost(object):
    '''
        A least squares cost function measuring the deviation of the force
        field hessian from the ab initio hessian. Only the force constants
        are regarded as variable parameters, the rest values are kept fixed.
        The initial force constants and fixed rest values should be stored
        in attributes of the members of model.vterms.

        The cost function measures half the sum of the squares of the
        difference between ab initio and force field hessian elements:

            chi^2 = 0.5*sum( (H^ai_ij - H^ff_ij)^2 , ij=...)

        This cost function can be rewritten as:

            chi^2(k) = 0.5*np.dot(k.T, np.dot(A, k)) - np.dot(B.T, k) + 0.5*C

        in which k is the vector of unknown force constants.
    '''
    def __init__(self, system, model):
        '''
            **Arguments*

            system
                An instance of the System class containing all system info

            model
                An instance of the Model class defining the total energy,
                the electrostatic contribution and the valence terms.
        '''
        self.system = system
        self.model = model

    def _update_lstsq_matrices(self):
        '''
            Calculate the matrix A, vector B and scalar C of the cost function
        '''
        ref  = self.system.ref.hess.reshape([3*self.system.natoms, 3*self.system.natoms])
        ref -= self.model.ei.calc_hessian(self.system.ref.coords).reshape([3*self.system.natoms, 3*self.system.natoms])
        self.A = np.zeros([self.model.val.nterms, self.model.val.nterms], float)
        self.B = np.zeros([self.model.val.nterms], float)
        h = np.zeros([self.model.val.nterms, 3*self.system.natoms, 3*self.system.natoms], float)
        for i, icname in enumerate(sorted(self.system.ics.keys())):
            for vterm in self.model.val.vterms[icname]:
                h[i] += vterm.calc_hessian(coords=self.system.ref.coords, k=1.0)
            self.B[i] = matrix_squared_sum(h[i], ref)
            for j in xrange(i+1):
                 self.A[i,j] = matrix_squared_sum(h[i], h[j])
                 self.A[j,i] = self.A[i,j]
        del h
        self.C = matrix_squared_sum(ref, ref)

    def _define_positive_constraints(self):
        constraints = []
        def jac(i):
            result = np.zeros(self.model.val.nterms, float)
            result[i] = 1.0
            return result
        for i in xrange(self.model.val.nterms):
            constraint = {
                'type': 'ineq',
                'fun': lambda k: k[i],
                'jac': lambda k: jac(i),
            }
            constraints.append(constraint)
        return tuple(constraints)

    def fun(self, k, do_grad=False):
        chi2 = 0.5*np.dot(k.T, np.dot(self.A, k)) - np.dot(self.B.T, k) +0.5*self.C
        if do_grad:
            gchi2 = np.dot(self.A, k) - self.B
            return chi2, gchi2
        else:
            return chi2

    def estimate(self, tol=1e-9):
        '''
            Estimate the force constants by minimizing the cost function
        '''
        self._update_lstsq_matrices()
        kinit = self.model.val.get_fcs()
        constraints = self._define_positive_constraints()
        result = minimize(self.fun, kinit, method='SLSQP', constraints=constraints, tol=tol)
        return result.x
