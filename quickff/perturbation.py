# -*- coding: utf-8 -*-
# QuickFF is a code to quickly derive accurate force fields from ab initio input.
# Copyright (C) 2012 - 2016 Louis Vanduyfhuys <Louis.Vanduyfhuys@UGent.be>
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

from molmod.io.xyz import XYZWriter
from molmod.units import *

from yaff.pes.ff import ForceField, ForcePartValence
from yaff.pes.vlist import Chebychev1

from quickff.tools import fitpar

import numpy as np, scipy.optimize, warnings
warnings.filterwarnings('ignore', 'The iteration is not making good progress')

__all__ = ['Trajectory', 'RelaxedStrain']

class Trajectory(object):
    '''
        A class to store a perturbation trajectory
    '''
    def __init__(self, index, name, numbers, start, end, steps=11, qunit='au', kunit='au'):
        '''
            **Arguments**
            
            index
                the index in the valence force field vlist of the valence term 
                associated with the current trajectory 
            
            name
                the name of the valence term associated with the current
                trajectory
            
            numbers
                a numpy array containing the atomic numbers of the system
            
            start
                a float defining the lower limit of the perturbation value of
                the given ic. If not given, a standard value is choosen
                according to the kind of internal coordinate.

            end
                a float defining the upper limit of the perturbations value of
                the given ic. If not given, a standard value is choosen
                according to the kind of internal coordinate.

            **Optional Arguments**
            
            steps
                an integer defining the number of steps in the perturbation
                trajectory. The default value is 11 steps.
            
            qunit
                the unit of the internal coordinate value
            
            kunit
                the unit of the force constant associated with the IC
        '''
        self.index = index
        self.name = name
        self.numbers = numbers
        self.qunit = qunit
        self.kunit = kunit
        self.qvals = start + (end-start)/(steps-1)*np.array(range(steps), float)
        self.coords =  np.zeros([steps, len(numbers), 3], float)
        self.fc = None
        self.rv = None
    
    def plot(self, fn, refs, eunit='kjmol'):
        '''
            Method to plot the energy contributions along a perturbation
            trajectory associated to a given ic. This method assumes that the
            coords, fc and rv attributes have been computed and assigned.

            **Arguments**
            
            fn
                a string defining the name of the figure
            
            refs
                a list of Reference instances. The value of each of these
                contributions will be plotted.

            **Optional Arguments**

            eunit
                a string describing the conversion of the unit of energy. More
                info regarding possible strings can be found in the
                `MolMod documentation <http://molmod.github.io/molmod/reference/const.html#module-molmod.units>`_.
        '''
        import matplotlib.pyplot as pp
        fig, ax = pp.subplots()
        #plot reference data
        for iref, ref in enumerate(refs):
            #collect data
            data = np.zeros([len(self.coords)], float)
            for istep, pos in enumerate(self.coords):
                data[istep] = ref.energy(pos)
            data[:] -= data[len(data)/2]
            #configure labels, linestyle, color
            pars = fitpar(self.qvals, data, rcond=1e-6)
            k = 2*pars[0]
            if k != 0.0:
                q0 = -0.5*pars[1]/pars[0]
                label = '(K=%.0f q0=%.3f)' %(
                    k/parse_unit(self.kunit),
                    q0/parse_unit(self.qunit)
                )
            else:
                label = '(K=0.0 q0=None)'
            if 'ai' in ref.name.lower():
                linestyle = 'k-'
                linewidth = 4.0
            else:
                linestyle = ':'
                linewidth = 2.0
            #plot data
            ax.plot(
                self.qvals/parse_unit(self.qunit), data/parse_unit(eunit), linestyle, 
                linewidth=linewidth, label='%8s %s' %(ref.name, label)
            )
        #plot current term 
        cov = 0.5*self.fc*(self.qvals - self.rv)**2
        cov -= cov[len(cov)/2]
        label = '(K=%.0f q0=%.3f)' %(
            k/parse_unit(self.kunit),
            q0/parse_unit(self.qunit)
        )
        ax.plot(
            self.qvals/parse_unit(self.qunit), cov/parse_unit(eunit),
            'r-', linewidth=2, label='FF cov   %s' %covlabel
        )
        #decorate plot with title, labels, ...
        ax.set_title(self.name)
        ax.set_xlabel('%s [%s]' % (self.name.split('/')[0], self.qunit), fontsize=16)
        ax.set_ylabel('Energy [%s]' %eunit, fontsize=16)
        ax.grid()
        ax.legend(loc='best', fontsize=16)
        fig.set_size_inches([8, 8])
        fig.savefig(fn)
    
    def to_xyz(self, fn):
        '''
            Method to write the given trajectory to an XYZ file. This method
            assumes that the coords attribute has been assigned.

            **Arguments**

            fn
                a string defining the name of the output file
        '''
        assert self.coords is None, 'Coords attribute not assigned'
        f = open(filename, 'w')
        xyz = XYZWriter(f, [pt[Z].symbol for Z in self.numbers])
        for frame, coord in enumerate(self.coords):
            xyz.dump('frame %i' %frame, coords)
        f.close()




class BasePerturbationTheory(object):
    def __init__(self, system, refs, valence):
        self.system = system
        self.ai = None
        self.refs = []
        for ref in refs:
            if 'ai' in ref.name.lower():
                self.ai = ref
            else:
                self.refs.append(ref)
        if self.ai is None:
            raise IOError("No Ab Initio reference found. Be sure to add the string 'ai' to its name.")
        self.valence = valence

    def prepare(self):
        '''
            Method to initialize trajectories and configure everything required
            for the generate method.
        '''
        #set system coords to ab initio reference and compute the values of
        #all ics
        self.system.pos = self.ai.coords.copy()
        self.valence.dlist.forward()
        self.valence.iclist.forward()
        #list of diagonal terms for which the perturbation trajectories will be
        #computed also set fc to 1.0 and rv to ai equilibrium of each harmonic
        #term
        trajectories = []
        for index in xrange(self.valence.vlist.nv):
            term = self.valence.vlist.vtab[index]
            ic = self.valence.iclist.ictab[term['ic0']]
            self.valence.set_params(index,fc=1.0,rv0=ic['value'])
            if 'PT_ALL' in self.valence.data[index].tasks:
                if 'EQ_RV' in self.valence.data[index].tasks:
                    raise ValueError, 'Both PT_ALL and EQ_RV are in the tasks of %s' \
                        % self.valence.data[index].name
                if term['kind']==0:#harmonic
                    kunit, qunit = self.valence.data[index].units
                else:
                    raise ValueError('Perturbation trajectory not supported for terms of kind %i' %term['kind'])
                if ic['kind'] in [0, 10]:
                    start=ic['value']-0.05*angstrom
                    end=ic['value']+0.05*angstrom
                elif ic['kind'] in [2, 4]:
                    start=ic['value']-5*deg
                    end=ic['value']+5*deg
                    if end>180*deg: end=180*deg
                    if start<0*deg:
                        raise ValueError(
                            'Starting angle for %s is smaller than zero' %(
                                self.valence.data[index].name
                            )
                        )
                else:
                    raise NotImplementedError
                traj = Trajectory(
                    index, self.valence.data[index].name+'('+str(index)+')',
                    self.system.numbers, start, end, steps=11, qunit=qunit, 
                    kunit=kunit
                )
                trajectories.append(traj)
        return trajectories
    
    def generate(self, trajectory):
        '''
            Method to calculate the perturbation trajectory, i.e. the trajectory
            that scans the geometry along the direction of the ic figuring in
            the term with the given index (should be a diagonal term). This
            method should be implemented in the derived classes.

            **Arguments**

            trajectory
                a Trajectory instance representing the perturbation trajectory
        '''
        raise NotImplementeError



    def estimate(self, trajectory):
        '''
            Method to estimate the FF parameters for the given ic from the given
            perturbation trajectory by fitting a harmonic potential to the
            covalent energy along the trajectory.

            **Arguments**

            trajectory
                a Trajectory instance representing the perturbation trajectory
        '''
        qs = trajectory.qvals.copy()
        Es = np.zeros(len(trajectory.coords))
        for istep, pos in enumerate(trajectory.coords):
            Es[istep] += self.ai.energy(pos)
            for ref in self.refs:
                Es[istep] -= ref.energy(pos)
        pars = fitpar(qs, Es, rcond=1e-6)
        if pars[0]!=0.0:
            trajectory.fc = 2*pars[0]
            trajectory.rv = -pars[1]/(2*pars[0])
        else:
            trajectory.fc = 0
            trajectory.rv = qs[len(qs)/2]

    def refinervs(self, trajectory):
        '''
            Method to refine the rest values of the covalent terms assuming the
            force constants are known.

            **Arguments**

            trajectory
                a Trajectory instance representing the perturbation trajectory
        '''
        term = self.valence.vlist.vtab[trajectory.index]
        kold, q0old = trajectory.fc, trajectory.rv
        if kold<1.0*kjmol:
            return kold, q0old
        else:
            Es = np.zeros(len(trajectory.coords))
            qs = trajectory.qvals.copy()
            for istep, pos in enumerate(trajectory.coords):
                Es[istep] += self.ai.energy(pos)
                for ref in self.refs:
                    Es[istep] -= ref.energy(pos)
                #substract valence energy minus the current term energy
                Es[istep] -= (self.valence.energy(pos) - term['energy'])
            pars = fitpar(qs, Es, rcond=1e-6)
            q0new = -0.5*pars[1]/pars[0]
            if q0new<0.0:
                q0new = 0.0
            trajectory.rv = q0new


class RelaxedStrain(BasePerturbationTheory):
    def __init__(self, system, refs, valence, do_taylor=False):
        self.do_taylor = do_taylor
        BasePerturbationTheory.__init__(self, system, refs, valence)
    
    def generate(self, trajectory):
        term = self.valence.vlist.vtab[trajectory.index]
        strain = Strain(self.system, self.valence, trajectory.index)
        if self.do_taylor:
            #TODO: extend for taylor strain
            raise NotImplementedError
        for iq, q in enumerate(trajectory.qvals):
            strain.constrain_value = q
            guess = np.zeros((self.system.natom*3+1,))
            sol = scipy.optimize.fsolve(strain.gradient, guess, xtol=1e-9)
            trajectory.coords[iq,:,:] = self.ai.coords + sol[:self.system.natom*3].reshape((-1,3))


class Strain(object):
    def __init__(self, system, valence, index):
        self.coords0 = system.pos.copy()
        self.ndof = np.prod(self.coords0.shape)
        self.strain = ForceField(system, [valence])
        self.constrain_value = None
        constraintpart = ForcePartValence(system)
        # Abuse the Chebychev1 polynomial to simply get the value of q-1
        assert len(valence.data[index].ics)==1#should not be cross term
        constraintpart.add_term(Chebychev1(-2.0,valence.data[index].ics[0]))
        self.constraint = ForceField(system, [constraintpart])

    def value(self, X):
        self.strain.update_pos(self.refdata.coords + X[:self.ndof].reshape((-1,3)))
        return self.strain.compute()

    def gradient(self, X):
        '''
            Compute the gradient of the strain wrt Cartesian coordinates of the
            system. For every ic that needs to be constrained, a Lagrange multiplier
            is included.
        '''
        #small check
        assert X.shape[0] ==  self.ndof + 1
        #initialize return value
        grad = np.zeros((len(X),))
        #compute strain gradient
        gstrain = np.zeros(self.coords0.shape)
        self.strain.update_pos(self.coords0 + X[:self.ndof].reshape((-1,3)))
        self.strain.compute(gpos=gstrain)
        #compute constraint gradient
        gconstraint = np.zeros(self.coords0.shape)
        self.constraint.update_pos(self.coords0 + X[:self.ndof].reshape((-1,3)))
        self.constraint.compute(gpos=gconstraint)
        #construct gradient
        grad[:self.ndof] = gstrain.reshape((-1,)) + X[self.ndof]*gconstraint.reshape((-1,)) #+ 0.01*X[:self.ndof]
        grad[self.ndof] = self.constraint.parts[0].vlist.vtab[0]['energy'] + 1.0 - self.constrain_value
        return grad
