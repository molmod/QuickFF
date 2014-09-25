# -*- coding: utf-8 -*-
#QuickFF is a code to quickly derive accurate force fields from ab initio input.
#Copyright (C) 2012 - 2014 Louis Vanduyfhuys <Louis.Vanduyfhuys@UGent.be>
#Steven Vandenbrande <Steven.Vandenbrande@UGent.be>,
#Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center for Molecular Modeling
#(CMM), Ghent University, Ghent, Belgium; all rights reserved unless otherwise
#stated.
#
#This file is part of QuickFF.
#
#QuickFF is free software; you can redistribute it and/or
#modify it under the terms of the GNU General Public License
#as published by the Free Software Foundation; either version 3
#of the License, or (at your option) any later version.
#
#QuickFF is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, see <http://www.gnu.org/licenses/>
#
#--

from molmod.units import angstrom, deg, parse_unit
from molmod.periodic import periodic as pt
from molmod.io.xyz import XYZWriter

from scipy.optimize import minimize
import numpy as np

from quickff.tools import fitpar
from quickff.evaluators import *

__all__ = ['BasePertTheory', 'RelaxedGeoPertTheory']


class BasePertTheory(object):
    '''
       Base class to generate and analyze the perturbation trajectories of an
       ic.
    '''
    def __init__(self, system, model):
        '''
            **Arguments**

            system
                an instance of the System class

            model
                an instance of the Model class
        '''
        self.system = system
        self.model = model

    def generate(self, ic, start=None, end=None, steps=11):
        '''
            This method should be implemented in derived classes.
        '''
        raise NotImplementedError

    def analyze(self, trajectory, evaluators):
        '''
            A method to analyze the perturbation trajectory, i.e. apply all
            evaluators to the trajectory and return an array of values along
            the trajectory for each evaluator given.

            **Arguments**

            trajectory
                A (F, N, 3) numpy array defining the perturbation trajectory.
                It contains F frames of (N,3)-dimensional geometry arrays.

            evaluators
                A list of evaluators that is to be applied to the trajectory.
        '''
        values = [[] for i in xrange(len(evaluators))]
        for dx in trajectory:
            coords = self.system.ref.coords + dx
            for i, evaluator in enumerate(evaluators):
                values[i].append(evaluator(self.model, coords))
        return np.array(values)

    def write(self, trajectory, filename):
        '''
            Method to write the given trajectory to a file

            **Arguments**

            trajectory
                a (F,N,3) numpy array defining the perturbation trajectory.
                It contains F frames of (N,3)-dimensional geometry arrays.

            filename
                a string defining the name of the output file
        '''
        _f = open(filename, 'w')
        xyz = XYZWriter(_f, [pt[Z].symbol for Z in self.system.numbers])
        for idx, dx in enumerate(trajectory):
            coords = self.system.ref.coords + dx
            xyz.dump('frame %i' %idx, coords)
        _f.close()

    def plot(self, ic, trajectory, filename, eunit='kjmol'):
        '''
            Method to plot the energy contributions along a perturbation
            trajectory associated to a given ic.

            **Arguments**

            ic
                an instance of the class :class:`quickff.ic.IC` defining for
                which ic the plot will be made.

            trajectory
                a (F,N,3) numpy array defining the perturbation trajectory
                associated to the given ic. It contains F frames of
                (N,3)-dimensional geometry arrays.

            filename
                a string defining the name of the figure

            **Optional Arguments**

            eunit
                a string describing the conversion of the unit of energy. More
                info regarding possible strings can be found in the
                `MolMod documentation <http://molmod.github.io/molmod/reference/const.html#module-molmod.units>`_.
        '''
        import matplotlib.pyplot as pp
        evaluators = [eval_ic(ic), eval_energy('ai'), eval_energy('ei'), eval_energy('vdw')]
        qs, tot, ei, vdw = self.analyze(trajectory, evaluators)
        label = {}
        for (name, obs) in zip(['ai', 'ei', 'vdw', 'cov'], [tot, ei, vdw, tot-ei-vdw]):
            pars = fitpar(qs, obs, rcond=1e-6)
            k = 2*pars[0]
            if k != 0.0:
                q0 = -0.5*pars[1]/pars[0]
                label[name] = '(K=%.0f q0=%.3f)' %(
                    k/parse_unit(ic.kunit),
                    q0/parse_unit(ic.qunit)
                )
            else:
                label[name] = '(K=0.0 q0=None)'
            if name=='cov':
                cov = (pars[0]*qs**2+pars[1]*qs+pars[2])
        fig, ax = pp.subplots()
        ax.plot(
            qs/parse_unit(ic.qunit), tot/parse_unit(eunit),
            'k--', linewidth=4, label='AI total %s' %label['ai']
        )
        ax.plot(
            qs/parse_unit(ic.qunit), ei/parse_unit(eunit),
            'b--', linewidth=2, label='FF elec  %s' %label['ei']
        )
        ax.plot(
            qs/parse_unit(ic.qunit), vdw/parse_unit(eunit),
            'g--', linewidth=2, label='FF vdW   %s' %label['vdw']
        )
        ax.plot(
            qs/parse_unit(ic.qunit), cov/parse_unit(eunit),
            'r-', linewidth=2, label='FF cov   %s' %label['cov']
        )
        ax.set_title(ic.name)
        ax.set_xlabel('%s [%s]' % (ic.name.split('/')[0], ic.qunit), fontsize=16)
        ax.set_ylabel('Energy [%s]' %eunit, fontsize=16)
        ax.grid()
        ax.legend(loc='best', fontsize=16)
        fig.set_size_inches([8, 8])
        fig.savefig(filename)

    def estimate(self, ic, trajectory):
        '''
            Method to estimate the FF parameters for the given ic from the given
            perturbation trajectory by fitting a harmonic potential to the
            covalent energy along the trajectory.

            **Arguments**

            ic
                an instance of the class :class:`quickff.ic.IC` defining for
                which ic the FF parameters will be estimated

            trajectory
                a (F,N,3) numpy array defining the perturbation trajectory
                associated to the given ic. It contains F frames of
                (N,3)-dimensional geometry arrays.
        '''
        evaluators = [eval_ic(ic), eval_energy('ai'), eval_energy('ei'), eval_energy('vdw')]
        qs, tot, ei, vdw = self.analyze(trajectory, evaluators)
        pars = fitpar(qs, tot-ei-vdw, rcond=1e-6)
        return 2*pars[0], -pars[1]/(2*pars[0])


class RelaxedGeoPertTheory(BasePertTheory):
    '''
        Class for generating relaxed perturbation trajectories. These
        trajectories are constructed as an array of geometries that are
        perturbed in the direction of a given ic and relaxed in all other
        directions. The relaxation is implemented as the minimization of the
        strain due to all other ics.
    '''
    def __init__(self, system, model):
        '''
            **Arguments**

            system
                an instance of the System class

            model
                an instance of the Model class
        '''
        BasePertTheory.__init__(self, system, model)

    def get_strain_matrix(self, ic):
        '''
            Method to calculate the strain matrix.

            If sandwiched between a geometry perturbation vector, this
            represents the weighted sum of the deviations of the ics
            from their equilibrium values, except for the ic given in args.

            **Arguments**

            ic
                an instance of the class :class:`quickff.ic.IC` defining which
                ic is to be excluded from the strain cost.
        '''
        ndofs = 3*self.system.natoms
        strain = np.zeros([ndofs, ndofs], float)
        for icname, ics in sorted(self.system.ics.iteritems()):
            if ic is None:
                Gq = np.array([ic0.grad(self.system.ref.coords) for ic0 in ics])
            else:
                Gq = np.array([ic0.grad(self.system.ref.coords) for ic0 in ics if ic0.name!=ic.name])
            if len(Gq) > 0:
                U, S, Vt = np.linalg.svd(Gq)
                svals = np.array([1.0 for s in S if s > 1e-6])
                rank = len(svals)
                if rank > 0:
                    S2 = np.diag(svals**2)
                    V  = Vt.T[:, :rank]
                    Vo = Vt.T[:, rank:]
                    strain += np.dot(V, np.dot(S2, V.T)) \
                            + np.dot(Vo, Vo.T)/(100*ndofs)
                else:
                    #if the gradient of the current ic is zero in equilibrium,
                    #use the second order contribution to the Taylor expansion
                    #of the ic (without the square)
                    print '    WARNING: %s '  % icname + 'gradient is zero,' +\
                          'using second order Taylor to estimate strain'
                    Hq = np.zeros([ndofs, ndofs], float)
                    for ic0 in ics:
                        if ic is None or ic0.name != ic.name:
                            Hq += ic0.hess(self.system.ref.coords)
                    S, V = np.linalg.eigh(Hq)
                    S = np.diag([1.0 for s in S if abs(s) > 1e-6])
                    rank = len(S)
                    V = V[:, :rank]
                    strain += np.dot(V, np.dot(S, V.T))
        return strain

    def generate(self, ic, start=None, end=None, steps=11):
        '''
            Method to calculate the perturbation trajectory, i.e. the trajectory
            that arises when the geometry is perturbed in the direction of ic
            and relaxed in all other directions.

            **Arguments**

            ic
                an instance of the class :class:`quickff.ic.IC` that defines
                for which ic the perturbation trajectory will be calculated.

            **Optional Arguments**

            start
                a float defining the lower limit of the perturbation value of
                the given ic. If not given, a standard value is choosen
                according to the kind of internal coordinate.

            end
                a float defining the upper limit of the perturbations value of
                the given ic. If not given, a standard value is choosen
                according to the kind of internal coordinate.

            steps
                an integer defining the number of steps in the perturbation
                trajectory. The default value is 11 steps.
        '''
        ndofs = 3*self.system.natoms
        #initialization
        q0 = ic.value(self.system.ref.coords)
        if ic.name.startswith('bond'):
            if start is None: start = q0 - 0.05*angstrom
            if end is None:   end   = q0 + 0.05*angstrom
        elif ic.name.startswith('angle'):
            if start is None:   start = q0 - 5*deg
            if end is None:       end = q0 + 5*deg
            if end > 180.0*deg:   end = 179.0*deg
            if q0 == 180.0*deg:   end = 180.0*deg
        elif ic.name.startswith('dihed'):
            if start is None:      start = q0 - 5*deg
            if end is None:          end = q0 + 5*deg
            if start < -180.0*deg: start = -180.0*deg
            if end > 180.0*deg:      end = 180.0*deg
        elif ic.name.startswith('opdist'):
            if start is None: start = q0 - 0.05*angstrom
            if end is None:     end = q0 + 0.05*angstrom
        qarray = start + (end-start)/(steps-1)*np.array(range(steps), float)
        trajectory = np.zeros([steps, self.system.natoms, 3], float)
        #Define cost function that needs to be minimized
        S = self.get_strain_matrix(ic)
        def chi(dx):
            return 0.5*np.dot(dx.T, np.dot(S, dx))
        #Guess delta_x first time
        guess = np.zeros(ndofs, float)
        for iq, q in enumerate(qarray):
            #Only use if the minimum is located in 180*deg
            if ic.name.startswith('angle') and q == 180*deg and q0 == 180*deg:
                trajectory[iq] = np.zeros([self.system.natoms, 3], float)
                continue
            #Define the constraint under which the cost function needs
            #to be minimized
            constraints = ({
                'type': 'eq',
                'fun' : lambda dx: ic.value(self.system.ref.coords + dx.reshape((-1, 3))) - q,
                'jac' : lambda dx: ic.grad(self.system.ref.coords + dx.reshape((-1, 3))),
            },)
            result = minimize(
                chi, guess, method='SLSQP',
                constraints=constraints,
                tol=1e-9
            )
            trajectory[iq] = result.x.reshape([-1, 3])
            #Use the result just found as the new guess to ensure continuity
            guess = result.x
        return trajectory
