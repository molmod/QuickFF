from molmod.units import angstrom, deg, parse_unit
from molmod.periodic import periodic as pt
from molmod.io.xyz import XYZWriter

from scipy.optimize import minimize
import numpy as np

from quickff.tools import fitpar
from quickff.evaluators import *

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
        for dx in trajectory:
            coords = self.system.ref.coords + dx
            for i, evaluator in enumerate(evaluators):
                values[i].append(evaluator(self.model, coords))
        return np.array(values)

    def write(self, trajectory, filename):
        _f = open(filename, 'w')
        xyz = XYZWriter(_f, [pt[Z].symbol for Z in self.system.numbers])
        for idx, dx in enumerate(trajectory):
            coords = self.system.ref.coords + dx
            xyz.dump('frame %i' %idx, coords)
        _f.close()

    def plot(self, ic, trajectory, filename, eunit='kjmol'):
        import matplotlib.pyplot as pp
        evaluators = [eval_ic(ic), eval_energy('ai'), eval_energy('ei')]
        qs, tot, ei = self.analyze(trajectory, evaluators)
        pars = fitpar(qs, tot-ei, rcond=1e-6)
        pp.clf()
        pp.plot(
            qs/parse_unit(ic.qunit),
            tot/parse_unit(eunit),
            'k--', linewidth=2, label='AI total'
        )
        pp.plot(
            qs/parse_unit(ic.qunit),
            ei/parse_unit(eunit),
            'b--', linewidth=1, label='FF electrostatic'
        )
        pp.plot(
            qs/parse_unit(ic.qunit),
            (pars[0]*qs**2+pars[1]*qs+pars[2])/parse_unit(eunit),
            'b-', linewidth=1, label='FF covalent fitted'
        )
        pp.title(ic.name)
        pp.xlabel('%s [%s]' % (ic.name.split('/')[0], ic.qunit), fontsize=16)
        pp.ylabel('Energy [%s]' %eunit, fontsize=16)
        pp.legend(loc='best', fontsize=16)
        fig = pp.gcf()
        fig.set_size_inches([8, 8])
        pp.savefig(filename)

    def estimate(self, ic, trajectory):
        evaluators = [eval_ic(ic), eval_energy('ai'), eval_energy('ei')]
        qs, tot, ei = self.analyze(trajectory, evaluators)
        pars = fitpar(qs, tot-ei, rcond=1e-6)
        return 2*pars[0], -pars[1]/(2*pars[0])


class RelaxedGeoPertTheory(BasePertTheory):
    def __init__(self, system, model):
        BasePertTheory.__init__(self, system, model)

    def get_strain_matrix(self, ic):
        '''
            Method to calculate the strain matrix.

            If sandwiched between a geometry perturbation vector, this
            represents the weighted sum of the deviations of the ics
            from their equilibrium values, except for the ic given in args.
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
                            + np.dot(Vo, Vo.T)/ndofs
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
            Calculate the perturbation trajectory, i.e. the trajectory that
            arises when the geometry is perturbed in the direction of ic
            and relaxed in all other directions.

            The relaxation is implemented as the minimization of weighted
            sum of the strain and the energy.
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
        #if ic.name.startswith('opdist'):
        #    guess = np.random.normal(loc=0.0, scale=0.01, size=ndofs)
        #else:
        #    guess = np.random.normal(loc=0.0, scale=0.000001, size=ndofs)
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
