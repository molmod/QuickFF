#! /usr/bin/env python

from molmod.units import *
from molmod.periodic import periodic as pt
from molmod.io.xyz import XYZWriter

from scipy.optimize import minimize
import numpy as np, matplotlib.pyplot as pp

from tools import *
from evaluators import *

__all__ = ['BasePertTheory', 'RelaxedGeoPertTheory']


class BasePertTheory(object):
    def __init__(self, system, model):
        '''
            Base class to generate and analyze the perturbation trajectories
            of an ic.

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

        **Arguments**

            ic
                An instance of the IC class

        **Optional Arguments**

            start
                a float defining the start value of the ic in the
                perturbation trajectory.

            end
                a float defining the end value of the ic in the
                perturbation trajectory.

            steps
                an integere defining the number of steps in the
                perturbation trajectory.
        '''
        raise NotImplementedError

    def analyze(self, trajectory, evaluators):
        values = [[] for i in xrange(len(evaluators))]
        for idx, dx in enumerate(trajectory):
            coords = self.system.ref.coords + dx
            for i, evaluator in enumerate(evaluators):
                values[i].append(evaluator(self.model, coords))
        return np.array(values)

    def write(self, trajectory, fn):
        f = open(fn, 'w')
        xyz = XYZWriter(f, [pt[number].symbol for number in self.system.numbers])
        for idx, dx in enumerate(trajectory):
            coords = self.system.ref.coords + dx
            xyz.dump('frame %i' %idx, coords)
        f.close()

    def plot(self, ic, trajectory, fn, eunit='kjmol'):
        evaluators = [eval_ic(ic), eval_energy('total'), eval_energy('ei')]
        qs, tot, ei = self.analyze(trajectory, evaluators)
        pp.clf()
        pp.plot(qs/parse_unit(ic.qunit), tot/parse_unit(eunit), 'k', linewidth=2)
        pp.plot(qs/parse_unit(ic.qunit), ei/parse_unit(eunit), 'b')
        pp.title(ic.name)
        pp.xlabel('%s [%s]' %(ic.name.split('/')[0], ic.qunit))
        pp.ylabel('Energy [%s]' %eunit)
        fig = pp.gcf()
        fig.set_size_inches([8,8])
        pp.savefig(fn)

    def estimate(self, ic, start=None, end=None, steps=11):
        trajectory = self.generate(ic, start, end, steps)
        evaluators = [eval_ic(ic), eval_energy('total'), eval_energy('ei')]
        qs, tot, ei = self.analyze(trajectory, evaluators)
        pars = fitpar(qs, tot-ei, rcond=1e-6)
        k = 2*pars[0]
        q0 = -pars[1]/k
        return k, q0


class RelaxedGeoPertTheory(BasePertTheory):
    def __init__(self, system, model, weight_strain=1.0, weight_energy=0.0):
        self.weight_strain = weight_strain
        self.weight_energy = weight_energy
        BasePertTheory.__init__(self, system, model)

    def get_strain_matrix(self, ic):
        '''
            Method to calculate the strain matrix.

            If sandwiched between a geometry perturbation vector, this
            represents the weighted sum of the deviations of the ics
            from their equilibrium values, except for the ic given in args.
        '''
        strain = np.zeros([3*self.system.natoms, 3*self.system.natoms], float)
        for icname, ics in sorted(self.system.ics.iteritems()):
            if ic is None:
                Gq = np.array([ic0.grad(self.system.ref.coords) for ic0 in ics])
            else:
                Gq = np.array([ic0.grad(self.system.ref.coords) for ic0 in ics if ic0.name!=ic.name])
            if len(Gq)>0:
                U, S, Vt = np.linalg.svd(Gq)
                svals = np.array([1.0 for s in S if s > 1e-6])
                rank = len(svals)
                if rank>0:
                    S2 = np.diag(svals**2)
                    V  = Vt.T[:,:rank]
                    Vo = Vt.T[:,rank:]
                    strain += np.dot(V, np.dot(S2, V.T)) + np.dot(Vo, Vo.T)/(3*self.system.natoms)
                else:
                    #if the gradient of the current ic is zero in equilibrium (e.g.: the cosine of an
                    #out-of-plane bend of a planar opbend-quadruple), use the second order contribution
                    #to the ic-Taylor expansion (without the square)
                    Hq = np.zeros([3*self.system.natoms, 3*self.system.natoms], float)
                    for ic0 in ics:
                        if ic is None or ic0.name!=ic.name:
                            Hq += ic0.hess(self.system.ref.coords)
                    S, V = np.linalg.eigh(Hq)
                    S = np.diag([1.0 for s in S if abs(s) > 1e-6])
                    rank = len(S)
                    V = V[:,:rank]
                    strain += np.dot(V, np.dot(S, V.T))
        return strain

    def generate(self, ic, start=None, end=None, steps=11):
        '''
            Calculate the perturbation trajectory, i.e. the trajectory that
            arises when the geometry is perturbed in the direction of ic
            and relaxed in all other directions.

            The relaxation is implemented as the minimization of weighted
            sum of the strain and the energy.
        '''
        #initialization
        q0 = ic.value(self.system.ref.coords)
        if ic.name.startswith('bond'):
            if start is None: start = q0 - 0.05*angstrom
            if end is None:   end   = q0 + 0.05*angstrom
        elif ic.name.startswith('angle'):
            if start is None: start = q0 - 5*deg
            if end is None:     end = q0 + 5*deg
            if end>180.0*deg:   end = 179.0*deg
            if q0==180.0*deg:   end = 180.0*deg
        elif ic.name.startswith('dihed'):
            if start is None:    start = q0 - 5*deg
            if end is None:        end = q0 + 5*deg
            if start<-180.0*deg: start = -180.0*deg
            if end>180.0*deg:      end = 180.0*deg
        elif ic.name.startswith('opbend'):
            if start is None: start = q0 - 0.05*angstrom
            if end is None:     end = q0 + 0.05*angstrom
        qarray = start + (end-start)/(steps-1)*np.array(range(steps),float)
        trajectory = np.zeros([steps, self.system.natoms, 3], float)
        #Define cost function that needs to be minimized
        H = self.model.total.hessian.reshape([3*self.system.natoms, 3*self.system.natoms])
        S = self.get_strain_matrix(ic)
        def chi(dx):
            strain = 0.5*np.dot(dx.T, np.dot(S, dx))
            energy = 0.5*np.dot(dx.T, np.dot(H, dx))
            return self.weight_strain*strain + self.weight_energy*energy
        #Guess delta_x first time
        if ic.name.startswith('opbend'):
            guess = np.random.normal(loc=0.0,scale=0.01,size=3*self.system.natoms)
        else:
            guess = np.random.normal(loc=0.0,scale=0.000001,size=3*self.system.natoms)
        for iq, q in enumerate(qarray):
            #Only use if the minimum is located in 180*deg
            if ic.name.startswith('angle') and q==180*deg and q0==180*deg:
                trajectory[iq] = np.zeros([self.system.natoms, 3], float)
                continue
            #Define the constraint under which the cost function needs to be minimized
            constraints = ({
                'type': 'eq',
                'fun' : lambda dx: ic.value(self.system.ref.coords + dx.reshape((-1, 3))) - q,
                'jac' : lambda dx: ic.grad(self.system.ref.coords + dx.reshape((-1, 3))),
            },)
            result = minimize(chi, guess, method='SLSQP', constraints=constraints, tol=1e-9)
            trajectory[iq] = result.x.reshape([-1, 3])
            #Use the result just found as the new guess to ensure continuity
            guess = result.x
        return trajectory
