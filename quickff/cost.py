# -*- coding: utf-8 -*-
# QuickFF is a code to quickly derive accurate force fields from ab initio input.
# Copyright (C) 2012 - 2015 Louis Vanduyfhuys <Louis.Vanduyfhuys@UGent.be>
# Steven Vandenbrande <Steven.Vandenbrande@UGent.be>,
# Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center for Molecular Modeling
# (CMM), Ghent University, Ghent, Belgium; all rights reserved unless otherwise
# stated.
#
# This file is part of QuickFF.
#
# QuickFF is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# QuickFF is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
#--

'''Quadratic cost function measuring deviation between ab initio and force-field hessian.

   The force constants are determined by minimizing this cost function while
   keeping the variables within certain boundaries.
'''

import numpy as np

from molmod.units import kjmol, angstrom
from yaff import ForcePartValence, Cosine, Harmonic, ForceField, estimate_cart_hessian

from quickff.tools import get_vterm, boxqp

__all__ = [
    'HessianFCCost',
]

class HessianFCCost(object):
    '''
        A least squares cost function measuring the deviation of the force
        field hessian from the ab initio hessian. Only the force constants
        are regarded as variable parameters, the rest values are kept fixed.
        The initial force constants and fixed rest values should be stored
        in a FFTable object

        The cost function measures half the sum of the squares of the
        difference between ab initio and force field hessian elements:

            chi^2 = 0.5*sum( (H^ai_ij - H^ff_ij)^2 , ij=...)

        This cost function can be rewritten as:

            chi^2(k) = 0.5*np.dot(k.T, np.dot(A, k)) - np.dot(B.T, k) + 0.5*C

        in which k is the vector of unknown force constants.
    '''
    def __init__(self, system, refdata, iclist):
        '''
            **Arguments**

            system
                An instance of the Yaff System class containing all system info

            refdata
                ReferenceData object containing reference coordinates, gradient
                and hessian
        '''
        self.system = system
        self.refdata = refdata
        self.iclist = iclist

    def _update_lstsq_matrices(self, fftab):
        '''
            Calculate the matrix A, vector B and scalar C of the cost function
        '''
        ndofs = 3*self.system.natom
        #TODO Figure out difference between phess and hess on final force field
        ref  = self.refdata.phess.reshape([ndofs, ndofs]).copy()
        if self.refdata.ncff is not None:
            ref -= self.refdata.ff_hessian.reshape([ndofs, ndofs]).copy()
        nterms = len(fftab.pars.keys())
        self._A = np.zeros([nterms, nterms], float)
        self._B = np.zeros([nterms], float)
        h = np.zeros([nterms, ndofs, ndofs], float)
        for i, icname in enumerate(sorted(fftab.pars.keys())):
            num_hessian = False
            val = ForcePartValence(self.system)
            ic_id = np.where(self.iclist.icnames==icname)[0][0]
            for iic in np.where(self.iclist.icname_ids==ic_id)[0]:
                if icname.startswith('dihed'):
                    pars = [fftab.pars[icname]['m'].mean,1.0,fftab.pars[icname]['q0'].mean]
                    indexes = self.iclist.ics[iic].index_pairs
                    term = get_vterm(pars, [indexes[0][1],indexes[0][0],indexes[2][0],indexes[2][1]])
                    if term is None:
                        # Requested a dihedral term that can not be expressed as
                        # a polynomial in cos(phi), but Yaff hessian is not stable
                        # for such terms. Resort to numerical hessian
                        num_hessian = True
                        term = Cosine(pars[0], pars[1], pars[2], self.iclist.ics[iic])
                    val.add_term(term)
                else:
                    val.add_term(Harmonic(1.0, fftab.pars[icname]['q0'].mean, self.iclist.ics[iic]))
            ff = ForceField(self.system,[val])
            ff.update_pos(self.refdata.coords)
            if num_hessian:
                h[i] = estimate_cart_hessian(ff)
            else:
                ff.compute(hess=h[i])
            self._B[i] = np.sum(h[i]*ref)
            for j in xrange(i+1):
                self._A[i, j] = np.sum(h[i]*h[j])
                self._A[j, i] = self._A[i, j]
        self._C = np.sum(ref*ref)
        return h

    def fun(self, k, do_grad=False):
        '''
            Calculate the actual cost

            **Optional Arguments**

            do_grad
                also calculate and return the analytical gradient of the cost
                function towards the parameters
        '''
        chi2 = 0.5*np.dot(k.T, np.dot(self._A, k)) \
             - np.dot(self._B.T, k) + 0.5*self._C
        if do_grad:
            gchi2 = np.dot(self._A, k) - self._B
            return chi2, gchi2
        else:
            return chi2

    def estimate(self, fftab):
        '''
            Estimate the force constants by minimizing the cost function

        '''
        h = self._update_lstsq_matrices(fftab)
        kinit = []
        for icname in sorted(fftab.pars.keys()):
            if icname.startswith('dihed'):
                kinit.append(10.0*kjmol)
            else:
                kinit.append(fftab.pars[icname]['k'].mean)
        kinit = np.asarray(kinit)
        bndl = np.zeros((len(fftab.pars.keys())))
        bndu = np.zeros((len(fftab.pars.keys()))) + np.inf
        for i, icname in enumerate(sorted(fftab.pars.keys())):
            if icname.startswith('dihed'):
                print i, icname
                bndu[i] = 200.0*kjmol
        x = boxqp(self._A, self._B, bndl, bndu, kinit)
        for i, icname in enumerate(sorted(fftab.pars.keys())):
            fftab.pars[icname]['k'].data[:] = x[i]
            fftab.pars[icname]['k'].update_statistics()
        return fftab


class BaseConstraint(object):
    def __init__(self, ctype, npars):
        '''
            A class for defining constraints in the minimalization of the cost.

            **Arguments**

            ctype
                the type of the constraint, can be 'eq' (equality, zero) or
                'ineq' (inequality, non-negative)

            npars
                the number of parameters of the cost function
        '''
        self.ctype = ctype
        self.npars = npars

    def __call__(self):
        'return the constraint in scipy.optimize.minimize format'
        return {
            'type': self.ctype,
            'fun' : self._fun(),
            'jac' : self._jac()
        }

    def _fun(self):
        '''
            Method to return the function defining the constraint. The arguments
            of the returned function should be the force constants.
        '''
        raise NotImplementedError

    def _jac(self):
        '''
            Method to return the jacobian of the function. The argument
            of the returned function should be the force constants.
        '''
        raise NotImplementedError


class FixedValueConstraint(BaseConstraint):
    def __init__(self, index, value, npars):
        '''
            A fixed value constraint

            **Arguments**

            index
                the index of the constrained fc

            value
                the fixed value of the constrained fc

            npars
                the number of parameters of the cost function
        '''
        self.index = index
        self.value = value
        BaseConstraint.__init__(self, 'eq', npars)

    def _fun(self):
        '''
            Method to return the function defining the constraint. The arguments
            of the returned function should be the force constants.
        '''
        return lambda k: k[self.index] - self.value

    def _jac(self):
        '''
            Method to return the jacobian of the function. The argument
            of the returned function should be the force constants.
        '''
        jac = np.zeros(self.npars, float)
        jac[self.index] = 1.0
        return lambda k: jac


class LowerLimitConstraint(BaseConstraint):
    def __init__(self, index, lower, npars):
        '''
            An upper limit constraint

            **Arguments**

            index
                the index of the constrained fc

            lower
                the upper limit of the constrained fc

            npars
                the number of parameters of the cost function
        '''
        self.index = index
        self.lower = lower
        BaseConstraint.__init__(self, 'ineq', npars)

    def _fun(self):
        '''
            Method to return the function defining the constraint. The arguments
            of the returned function should be the force constants.
        '''
        return lambda k: k[self.index] - self.lower

    def _jac(self):
        '''
            Method to return the jacobian of the function. The argument
            of the returned function should be the force constants.
        '''
        jac = np.zeros(self.npars, float)
        jac[self.index] = 1.0
        return lambda k: jac


class UpperLimitConstraint(BaseConstraint):
    def __init__(self, index, upper, npars):
        '''
            An upper limit constraint

            **Arguments**

            index
                the index of the constrained fc

            upper
                the upper limit of the constrained fc

            npars
                the number of parameters of the cost function
        '''
        self.index = index
        self.upper = upper
        BaseConstraint.__init__(self, 'ineq', npars)

    def _fun(self):
        '''
            Method to return the function defining the constraint. The arguments
            of the returned function should be the force constants.
        '''
        return lambda k: self.upper - k[self.index]

    def _jac(self):
        '''
            Method to return the jacobian of the function. The argument
            of the returned function should be the force constants.
        '''
        jac = np.zeros(self.npars, float)
        jac[self.index] = -1.0
        return lambda k: jac
