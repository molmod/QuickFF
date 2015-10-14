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

'''Tools to construact perturbation trajectories for internal coordinates.
'''


import warnings
warnings.filterwarnings('ignore', 'The iteration is not making good progress')
from scipy.optimize import minimize
from scipy.optimize import fsolve
import numpy as np
from datetime import datetime

from yaff import ForcePartValence, Harmonic, Chebychev1, ForceField
from molmod.units import angstrom, deg, parse_unit, kjmol
from molmod.periodic import periodic as pt
from molmod.io.xyz import XYZWriter
from molmod.ic import *
from molmod import parse_unit

from quickff.tools import fitpar

__all__ = ['BasePertTheory', 'RelaxedGeoPertTheory', 'Strain', 'StrainTaylor']


class BasePertTheory(object):
    '''
       Base class to generate and analyze the perturbation trajectories of an
       ic.
    '''
    def __init__(self, system, refdata, iclist):
        '''
            **Arguments**

            system
                an instance of the System class

            refdata
                an instance of the RefData class

            iclist
                an instance of the ICList class
        '''
        self.system = system
        self.refdata = refdata
        self.iclist = iclist

    def generate(self, ic, start=None, end=None, steps=11):
        '''
            This method should be implemented in derived classes.
        '''
        raise NotImplementedError

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
        for itraj, traj in enumerate(trajectory):
            xyz.dump('frame %i' %itraj, traj)
        _f.close()

    def plot(self, iic, trajectory, filename, eunit='kjmol'):
        '''
            Method to plot the energy contributions along a perturbation
            trajectory associated to a given ic.

            **Arguments**

            ic
                Index of the ic in the associated iclist for
                which the plot will be made.

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
        qs = np.zeros(len(trajectory))
        tot = np.zeros(len(trajectory))
        noncov = np.zeros(len(trajectory))
        for istep, pos in enumerate(trajectory):
            tot[istep] = self.refdata.calc_energy(pos, ff=False)
            noncov[istep] = self.refdata.calc_energy(pos, ai=False)
            self.system.pos[:] = pos
            self.iclist.dlist.forward()
            self.iclist.forward()
            qs[istep] = self.iclist.ictab[iic]['value']
        tot[:] -= tot[len(trajectory)/2]
        noncov[:] -= noncov[len(trajectory)/2]
        label = {}
        for (name, obs) in zip(['ai', 'noncov', 'cov'], [tot, noncov, tot-noncov]):
            pars = fitpar(qs, obs, rcond=1e-6)
            k = 2*pars[0]
            if k != 0.0:
                q0 = -0.5*pars[1]/pars[0]
                label[name] = '(K=%.0f q0=%.3f)' %(
                    k/parse_unit(self.iclist.kunits[self.iclist.icname_ids[iic]]),
                    q0/parse_unit(self.iclist.qunits[self.iclist.icname_ids[iic]])
                )
            else:
                label[name] = '(K=0.0 q0=None)'
            if name=='cov':
                cov = (pars[0]*qs**2+pars[1]*qs+pars[2])
        fig, ax = pp.subplots()
        ax.plot(
            qs/parse_unit(self.iclist.qunits[self.iclist.icname_ids[iic]]), tot/parse_unit(eunit),
            'k--', linewidth=4, label='AI total %s' %label['ai']
        )
        ax.plot(
            qs/parse_unit(self.iclist.qunits[self.iclist.icname_ids[iic]]), noncov/parse_unit(eunit),
            'b--', linewidth=2, label='FF noncov  %s' %label['noncov']
        )
        ax.plot(
            qs/parse_unit(self.iclist.qunits[self.iclist.icname_ids[iic]]), cov/parse_unit(eunit),
            'r-', linewidth=2, label='FF cov   %s' %label['cov']
        )
        ax.set_title(self.iclist.icnames[self.iclist.icname_ids[iic]])
        ax.set_xlabel('%s [%s]' % (self.iclist.icnames[self.iclist.icname_ids[iic]].split('/')[0], self.iclist.qunits[self.iclist.icname_ids[iic]]), fontsize=16)
        ax.set_ylabel('Energy [%s]' %eunit, fontsize=16)
        ax.grid()
        ax.legend(loc='best', fontsize=16)
        fig.set_size_inches([8, 8])
        fig.savefig(filename)

    def estimate(self, iic, trajectory):
        '''
            Method to estimate the FF parameters for the given ic from the given
            perturbation trajectory by fitting a harmonic potential to the
            covalent energy along the trajectory.

            **Arguments**

            iic
                The index of the ic in the associated iclist

            trajectory
                a (F,N,3) numpy array defining the perturbation trajectory
                associated to the given ic. It contains F frames of
                (N,3)-dimensional geometry arrays.
        '''
        Es = np.zeros(len(trajectory))
        qs = np.zeros(len(trajectory))
        for istep, pos in enumerate(trajectory):
            Es[istep] = self.refdata.calc_energy(pos)
            self.system.pos[:] = pos
            self.iclist.dlist.forward()
            self.iclist.forward()
            qs[istep] = self.iclist.ictab[iic]['value']
        pars = fitpar(qs, Es, rcond=1e-6)
        if pars[0]!=0.0:
            return 2*pars[0], -pars[1]/(2*pars[0])
        else: return 0.0, qs[len(qs)/2]

    def refineq(self, iic, trajectory, old, val_ff):
        '''
            Method to refine the rest values of the covalent terms assuming the
            force constants are known.

            **Arguments**

            iic
                The index of the ic in the associated iclist
            trajectory
                a (F,N,3) numpy array defining the perturbation trajectory
                associated to the given ic. It contains F frames of
                (N,3)-dimensional geometry arrays.
            old
                A 2-typle representing the previous force constant and rest
                value for the given ic.
        '''
        kold, q0old = old
        if kold<1.0*kjmol:
            return kold, q0old
        else:
            Es = np.zeros(len(trajectory))
            qs = np.zeros(len(trajectory))
            for istep, pos in enumerate(trajectory):
                val_ff.update_pos(pos)
                Es[istep] = self.refdata.calc_energy(pos) - val_ff.compute()
                self.system.pos[:] = pos
                self.iclist.dlist.forward()
                self.iclist.forward()
                qs[istep] = self.iclist.ictab[iic]['value']
            pars = fitpar(qs, Es, rcond=1e-6)
            q0new = -0.5*pars[1]/pars[0]
            if q0new<0.0:
                q0new = 0.0
            return kold, q0new


class RelaxedGeoPertTheory(BasePertTheory):
    '''
        Class for generating relaxed perturbation trajectories. These
        trajectories are constructed as an array of geometries that are
        perturbed in the direction of a given ic and relaxed in all other
        directions. The relaxation is implemented as the minimization of the
        strain due to all other ics.
    '''
    def __init__(self, system, refdata, iclist, taylor=False):
        '''
            **Arguments**

            system
                an instance of the System class

            model
                an instance of the Model class
        '''
        BasePertTheory.__init__(self, system, refdata, iclist)
        self.taylor = taylor

    def generate(self, iic, start=None, end=None, steps=11):
        '''
            Method to calculate the perturbation trajectory, i.e. the trajectory
            that arises when the geometry is perturbed in the direction of ic
            and relaxed in all other directions.

            **Arguments**

            iic
                The index of the ic in the associated iclist for which the
                perturbation trajectory will be generated

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
        if self.taylor: strain = StrainTaylor(self.system, self.iclist, [iic], self.refdata)
        else: strain = Strain(self.system, self.iclist, [iic], self.refdata)
        self.system.pos[:] = self.refdata.coords
        self.iclist.dlist.forward()
        self.iclist.forward()
        q0 = self.iclist.ictab[iic]['value']
        if self.iclist.icnames[self.iclist.icname_ids[iic]].startswith('bond'):
            if start is None: start = q0 - 0.05*angstrom
            if end is None:   end   = q0 + 0.05*angstrom
        elif self.iclist.icnames[self.iclist.icname_ids[iic]].startswith('angle'):
            if start is None:   start = q0 - 5*deg
            if end is None:       end = q0 + 5*deg
            if end > 180.0*deg:   end = 179.0*deg
            if q0 == 180.0*deg:   end = 180.0*deg
        elif self.iclist.icnames[self.iclist.icname_ids[iic]].startswith('dihed'):
            if start is None:      start = q0 - 5*deg
            if end is None:          end = q0 + 5*deg
            if start < -180.0*deg: start = -180.0*deg
            if end > 180.0*deg:      end = 180.0*deg
        elif self.iclist.icnames[self.iclist.icname_ids[iic]].startswith('opdist'):
            if start is None: start = q0 - 0.05*angstrom
            if end is None:     end = q0 + 0.05*angstrom
        else: raise NotImplementedError
        qarray = start + (end-start)/(steps-1)*np.array(range(steps), float)
        trajectory = np.zeros([steps, self.system.natom, 3], float)
        for iq, q in enumerate(qarray):
            strain.constrain_values[0] = q
            # Minimize strain function
            guess = np.zeros((self.system.natom*3+1,))
            sol = fsolve(strain.gradient, guess, xtol=1e-9)
            trajectory[iq,:,:] = self.refdata.coords + sol[:self.system.natom*3].reshape((-1,3))
        k,q = self.estimate(iic,trajectory)
        return trajectory


class Strain(object):
    '''
        The strain is the sum of quadratic deviations of the internal coordinates
        from their values in the reference structure. Typically one or more ics
        are not included in the strain but rather constrained to a fixed value.
    '''
    def __init__(self, system, iclist, iics, refdata):
        '''
            **Arguments**
                system: a Yaff System object
                iclist: a list of all the relevant internal coordinates
                iics: list of indexes of internal coordinates that will not be
                    included in the strain. Typically the values of these
                    internal coordinates will be constrained
                refdata: ReferenceData object, necessary here to get the
                    reference coordinates.
        '''
        self.refdata = refdata
        self.system = system
        self.iclist = iclist
        self.iics = iics
        self.constrain_values = np.zeros((len(self.iics),))
        self.system.pos[:] = self.refdata.coords
        self.iclist.dlist.forward()
        self.iclist.forward()
        strainpart = ForcePartValence(self.system)
        for jic in xrange(self.iclist.nic):
            if jic in iics: continue
            icname = self.iclist.icnames[self.iclist.icname_ids[jic]]
            #TODO Figure out why it is best to throw away dihedrals
            if self.iclist.icnames[self.iclist.icname_ids[jic]].startswith('dihed'):
                continue
            strainpart.add_term(Harmonic(1.0,self.iclist.ictab[jic]['value'],self.iclist.ics[jic]))
        constraintpart = ForcePartValence(self.system)
        for iic in iics:
            # Abuse the Chebychev1 polynomial to simply get the value of q-1
            constraintpart.add_term(Chebychev1(-2.0,self.iclist.ics[iic]))
        self.strain = ForceField(self.system, [strainpart])
        self.constraint = ForceField(self.system, [constraintpart])

    def value(self, X):
        self.strain.update_pos(self.refdata.coords + X[:3*self.system.natom].reshape((-1,3)))
        return self.strain.compute()

    def gradient(self, X):
        '''
        Compute the gradient of the strain wrt Cartesian coordinates of the
        system. For every ic that needs to be constrained, a Lagrange multiplier
        is included.
        '''
        assert X.shape[0] == 3*self.system.natom + len(self.iics)
        func = np.zeros((len(X),))
        gstrain = np.zeros(self.system.pos.shape)
        gconstraint = np.zeros(self.system.pos.shape)
        self.strain.update_pos(self.refdata.coords + X[:3*self.system.natom].reshape((-1,3)))
        self.iclist.dlist.forward()
        self.iclist.forward()
        self.strain.compute(gpos=gstrain)
        self.constraint.update_pos(self.refdata.coords + X[:3*self.system.natom].reshape((-1,3)))
        self.constraint.compute(gpos=gconstraint)
        for i in xrange(self.constraint.parts[0].vlist.nv):
            func[-1-i] =  self.constraint.parts[0].vlist.vtab[0]['energy'] + 1.0 - self.constrain_values[i]
        func[:3*self.system.natom] = gstrain.reshape((-1,)) + X[-1]*gconstraint.reshape((-1,)) + 0.01*X[:3*self.system.natom]
        return func


class StrainTaylor(Strain):
    def __init__(self, system, iclist, iics, refdata):
        super(StrainTaylor, self).__init__(system, iclist, iics, refdata)
        self.S = self.get_strain_matrix()

    def get_strain_matrix(self):
        def fill_gpos(jacobian, iic):
            gpos = np.zeros(np.prod(self.system.pos.shape))
            ic_jac = jacobian[iic]
            kind = self.iclist.ictab[iic]['kind']
            # Loop over all relative vectors making up this internal coordinate
            index = self.iclist.ictab[iic]['i0']
            i,j = self.iclist.dlist.deltas[index]['i'],self.iclist.dlist.deltas[index]['j']
            d = self.iclist.ictab[iic]['i0']
            gpos[3*i:3*(i+1)] -= ic_jac[3*d:3*(d+1)]
            gpos[3*j:3*(j+1)] += ic_jac[3*d:3*(d+1)]
            if not kind in [0,5]:
                index = self.iclist.ictab[iic]['i1']
                i,j = self.iclist.dlist.deltas[index]['i'],self.iclist.dlist.deltas[index]['j']
                d = self.iclist.ictab[iic]['i1']
                gpos[3*i:3*(i+1)] -= ic_jac[3*d:3*(d+1)]
                gpos[3*j:3*(j+1)] += ic_jac[3*d:3*(d+1)]
            if not kind in [0,1,2,5]:
                index = self.iclist.ictab[iic]['i2']
                i,j = self.iclist.dlist.deltas[index]['i'],self.iclist.dlist.deltas[index]['j']
                d = self.iclist.ictab[iic]['i2']
                gpos[3*i:3*(i+1)] -= ic_jac[3*d:3*(d+1)]
                gpos[3*j:3*(j+1)] += ic_jac[3*d:3*(d+1)]
            return gpos
        ndofs = 3*self.system.natom
        strain = np.zeros([ndofs, ndofs], float)
        self.system.pos[:] = self.refdata.coords
        self.iclist.dlist.forward()
        self.iclist.forward()
        jacobian = self.iclist.jacobian()
        for icname in self.iclist.icnames:
            ic_id = np.where(self.iclist.icnames==icname)[0][0]
            Gq = np.array([fill_gpos(jacobian, iic) for iic in np.where(self.iclist.icname_ids==ic_id)[0] if not iic in self.iics])
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
                    raise NotImplementedError
                    #if the gradient of the current ic is zero in equilibrium,
                    #use the second order contribution to the Taylor expansion
                    #of the ic (without the square)
        return strain

    def value(self, X):
        return 0.5*np.dot(X.T, np.dot(self.S, X))

    def gradient(self, X):
        assert X.shape[0] == 3*self.system.natom + len(self.iics)
        func = np.zeros((len(X),))
        gconstraint = np.zeros(self.system.pos.shape)
        self.iclist.dlist.forward()
        self.iclist.forward()
        #TODO Currently the evaluation of this gradient is slow.
        #Typically the constraint force field will only contain one term, so
        #gradient can be calculated more efficiently than what is done now.
        self.constraint.update_pos(self.refdata.coords + X[:3*self.system.natom].reshape((-1,3)))
        self.constraint.compute(gpos=gconstraint)
        for i in xrange(self.constraint.parts[0].vlist.nv):
            func[-1-i] =  self.constraint.parts[0].vlist.vtab[0]['energy'] + 1.0 - self.constrain_values[i]
        func[:3*self.system.natom] = np.dot(self.S, X[:3*self.system.natom]) + X[-1]*gconstraint.reshape((-1,))
        return func
