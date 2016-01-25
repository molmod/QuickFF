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
from molmod.periodic import periodic as pt

from yaff.pes.ff import ForceField, ForcePartValence
from yaff.pes.vlist import Chebychev1, Harmonic
from yaff.pes.iclist import Bond, BendAngle, DihedAngle, OopDist

from quickff.tools import fitpar

import numpy as np, scipy.optimize, warnings
warnings.filterwarnings('ignore', 'The iteration is not making good progress')

__all__ = ['Trajectory', 'RelaxedStrain']

class Trajectory(object):
    '''
        A class to store a perturbation trajectory
    '''
    def __init__(self, term, start, end, numbers, steps=11):
        '''
            **Arguments**
            
            term
                an instance of the Term class for which the perturbation
                trajectory will be computed
            
            numbers
                list of atom numbers required for dumping xyz coords
            
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
        '''
        self.term = term
        self.numbers = numbers
        self.qunit = term.units[1]
        self.kunit = term.units[0]
        self.qvals = start + (end-start)/(steps-1)*np.array(range(steps), float)
        self.coords = np.zeros([steps, len(numbers), 3], float)
        self.fc = None
        self.rv = None
    
    def plot(self, refs, fn=None, eunit='kjmol', valence=None):
        '''
            Method to plot the energy contributions along a perturbation
            trajectory associated to a given ic. This method assumes that the
            coords, fc and rv attributes have been computed and assigned.

            **Arguments**
            
            refs
                a list of Reference instances. The value of each of these
                contributions will be plotted.

            **Optional Arguments**
            
            fn
                a string defining the name of the figure

            eunit
                a string describing the conversion of the unit of energy. More
                info regarding possible strings can be found in the
                `MolMod documentation <http://molmod.github.io/molmod/reference/const.html#module-molmod.units>`_.
           
            valence
                an instance of ValenceFF which will be used to plot the covalent
                contribution. If not given, only the contribution of the IC
                corresponding to the trajectory will be plotted using the values
                of fc and rv
        '''
        import matplotlib.pyplot as pp
        fig, ax = pp.subplots()
        #plot reference data
        ffrefs = np.zeros([len(self.coords)], float)
        for iref, ref in enumerate(refs):
            #collect data
            data = np.zeros([len(self.coords)], float)
            for istep, pos in enumerate(self.coords):
                data[istep] = ref.energy(pos)
            data[:] -= data[len(data)/2]
            if not 'ai' in ref.name.lower():
                ffrefs += data
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
                linewidth=linewidth, label='%8s %s' %(ref.name+' '*(8-len(ref.name)), label)
            )
        #plot covalent part
        if valence is not None:
            cov = np.zeros(len(self.qvals), float)
            for istep, pos in enumerate(self.coords):
                cov[istep] = valence.calc_energy(pos)
            pars = fitpar(self.qvals, cov, rcond=1e-6)
            k = 2*pars[0]
            if k != 0.0:
                q0 = -0.5*pars[1]/pars[0]
                label = '(K=%.0f q0=%.3f)' %(
                    k/parse_unit(self.kunit),
                    q0/parse_unit(self.qunit)
                )
            else:
                label = '(K=0.0 q0=None)'
        else:
            if self.fc is None or self.rv is None:
                raise ValueError('The values of fc and rv for %s(%i) are not set' %(self.term.basename, self.term.index))
            cov = 0.5*self.fc*(self.qvals - self.rv)**2
            label = '(K=%.0f q0=%.3f)' %(
                self.fc/parse_unit(self.kunit),
                self.rv/parse_unit(self.qunit)
            )
        cov -= cov[len(cov)/2]
        ax.plot(
            self.qvals/parse_unit(self.qunit), cov/parse_unit(eunit),
            'r-', linewidth=2, label='FF cov %s' %label
        )
        #plot ff total if valence was given
        if valence is not None:
            pars = fitpar(self.qvals, cov+ffrefs, rcond=1e-6)
            k = 2*pars[0]
            if k != 0.0:
                q0 = -0.5*pars[1]/pars[0]
                label = '(K=%.0f q0=%.3f)' %(
                    k/parse_unit(self.kunit),
                    q0/parse_unit(self.qunit)
                )
            else:
                label = '(K=0.0 q0=None)'
            ax.plot(
                self.qvals/parse_unit(self.qunit), (cov+ffrefs)/parse_unit(eunit),
                'k--', linewidth=2, label='FF tot %s' %label
            )
        #decorate plot with title, labels, ...
        ax.set_title(self.term.basename)
        ax.set_xlabel('%s [%s]' % (self.term.basename.split('/')[0], self.qunit), fontsize=16)
        ax.set_ylabel('Energy [%s]' %eunit, fontsize=16)
        ax.grid()
        ax.legend(loc='best', fontsize=16)
        fig.set_size_inches([8, 8])
        if fn is None: fn = self.term.basename.replace('/', '-')+'-%i.png'%self.term.index
        fig.savefig(fn)
    
    def to_xyz(self, fn=None):
        '''
            Method to write the given trajectory to an XYZ file. This method
            assumes that the coords attribute has been assigned.

            **Optional Arguments**

            fn
                a string defining the name of the output file
        '''
        if fn is None: fn = self.term.basename.replace('/', '-')+'-%i.xyz'%self.term.index
        f = open(fn, 'w')
        xyz = XYZWriter(f, [pt[Z].symbol for Z in self.numbers])
        for frame, coord in enumerate(self.coords):
            xyz.dump('frame %i' %frame, coord)
        f.close()



class RelaxedStrain(object):
    def __init__(self, system, refs, do_taylor=False):
        self.system = system
        if do_taylor:
            #TODO: extend for taylor strain
            raise NotImplementedError
        self.ai = None
        self.refs = []
        for ref in refs:
            if 'ai' in ref.name.lower():
                self.ai = ref
            else:
                self.refs.append(ref)
        if self.ai is None:
            raise IOError("No Ab Initio reference found. Be sure to add the string 'ai' to its name.")
        self.strains = []

    def prepare(self, valence, do_terms):
        '''
            Method to initialize trajectories and configure everything required
            for the generate method.
        '''
        trajectories = []
        self.strains = [None,]*len(valence.terms)
        for term in do_terms:
            assert term.kind==0, 'Only harmonic terms supported for pert traj'
            ic = valence.iclist.ictab[valence.vlist.vtab[term.index]['ic0']]
            kunit, qunit = term.units
            if ic['kind'] in [0, 10, 11]:
                start=ic['value']-0.05*angstrom
                end=ic['value']+0.05*angstrom
            elif ic['kind'] in [2, 4]:
                start=ic['value']-5*deg
                end=ic['value']+5*deg
                if end>180*deg: end=180*deg
                if start<0*deg:
                    raise ValueError('Starting angle for %s is smaller than zero' %(term.basename))
            else:
                raise NotImplementedError
            #collect all other ics
            ics = []
            for term2 in valence.terms:
                if term2.index == term.index: continue
                if term2.kind == 3: continue
                ics.append(term2.ics[0])
            self.strains[term.index] = Strain(term.index, self.system, term.ics[0], ics)
            traj = Trajectory(term, start, end, self.system.numbers, steps=11)
            trajectories.append(traj)
        return trajectories
    
    def generate(self, trajectory, remove_com=True):
        '''
            Method to calculate the perturbation trajectory, i.e. the trajectory
            that scans the geometry along the direction of the ic figuring in
            the term with the given index (should be a diagonal term). This
            method should be implemented in the derived classes.

            **Arguments**

            trajectory
                a Trajectory instance representing the perturbation trajectory
        '''
        strain = self.strains[trajectory.term.index]
        assert strain.term_index == trajectory.term.index
        for iq, q in enumerate(trajectory.qvals):
            strain.constrain_value = q
            guess = np.zeros([strain.ndof+1], float)
            sol = scipy.optimize.fsolve(strain.gradient, guess, xtol=1e-9)
            x = strain.coords0 + sol[:self.system.natom*3].reshape((-1,3))
            if remove_com:
                com = (x.T*self.system.masses).sum(axis=1)/self.system.masses.sum()
                for i in xrange(self.system.natom):
                    x[i,:] -= com
            trajectory.coords[iq,:,:] = x

    def estimate(self, trajectory):
        '''
            Method to estimate the FF parameters for the relevant ic from the
            given perturbation trajectory by fitting a harmonic potential to the
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
            trajectory.fc = 0.0
            trajectory.rv = qs[len(qs)/2]



class Strain(ForceField):
    def __init__(self, term_index, system, cons_ic, ics):
        '''
            **Arguments**
            
            term_index
                Integer defining the index of the term in valence.terms for
                which the current strain is designed.
            
            system
            
            cons_ic
                An instance of Yaff Internal Coordinate representing the
                constrained term in the strain
            
            ics
                A list of Yaff Internal Coordinate instances for which the
                strain needs to be minimized
        '''
        self.term_index = term_index
        self.coords0 = system.pos.copy()
        self.ndof = np.prod(self.coords0.shape)
        part = ForcePartValence(system)
        for ic in ics:
            part.add_term(Harmonic(1.0, None, ic))
        part.dlist.forward()
        part.iclist.forward()
        for iterm in xrange(part.vlist.nv):
            term = part.vlist.vtab[iterm]
            ic = part.iclist.ictab[term['ic0']]
            term['par1'] = ic['value']
        ForceField.__init__(self, system, [part])
        #Abuse the Chebychev1 polynomial to simply get the value of q-1 and
        #implement the contraint
        part = ForcePartValence(system)
        part.add_term(Chebychev1(-2.0,cons_ic))
        self.constraint = ForceField(system, [part])
        self.constraint_value = None

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
        self.update_pos(self.coords0 + X[:self.ndof].reshape((-1,3)))
        self.compute(gpos=gstrain)
        #compute constraint gradient
        gconstraint = np.zeros(self.coords0.shape)
        self.constraint.update_pos(self.coords0 + X[:self.ndof].reshape((-1,3)))
        self.constraint.compute(gpos=gconstraint)
        #construct gradient
        grad[:self.ndof] = gstrain.reshape((-1,)) + X[self.ndof]*gconstraint.reshape((-1,)) #+ 0.01*X[:self.ndof]
        grad[self.ndof] = self.constraint.parts[0].vlist.vtab[0]['energy'] + 1.0 - self.constrain_value
        return grad
