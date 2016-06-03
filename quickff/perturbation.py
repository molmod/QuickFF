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
from yaff.pes.iclist import Bond, BendAngle, BendCos, DihedAngle, OopDist

from quickff.tools import fitpar
from quickff.log import log

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
        if not term.kind in [0]:
            raise NotImplementedError('Perturbation trajectory only implemented for Harmonic terms')
        self.term = term
        self.numbers = numbers
        self.qunit = term.units[1]
        self.kunit = term.units[0]
        self.step = (end-start)/(steps-1)
        self.targets = start + (end-start)/(steps-1)*np.array(range(steps), float)
        self.values = np.zeros(steps, float)
        self.coords = np.zeros([steps, len(numbers), 3], float)
        self.active = True
        self.fc = None
        self.rv = None

    def plot(self, ai, ffrefs=[], valence=None, fn=None, eunit='kjmol'):
        '''
            Method to plot the energy contributions along a perturbation
            trajectory associated to a given ic. This method assumes that the
            coords, fc and rv attributes have been computed and assigned.

            **Arguments**

            ai
                an instance of the Reference representing the ab initio input

            **Optional Arguments**

            ffrefs
                a list of Reference instances representing possible a priori
                determined contributions to the force field (such as eg.
                electrostatics and van der Waals)

            valence
                an instance of ValenceFF which will be used to plot the covalent
                contribution. If not given, only the contribution of the IC
                corresponding to the trajectory will be plotted using the values
                of fc and rv

            fn
                a string defining the name of the figure

            eunit
                a string describing the conversion of the unit of energy. More
                info regarding possible strings can be found in the
                `MolMod documentation <http://molmod.github.io/molmod/reference/const.html#module-molmod.units>`_.
        '''
        import matplotlib.pyplot as pp
        if 'active' in self.__dict__.keys() and not self.active: return
        fig, ax = pp.subplots()
        def add_plot(xs, ys, prefix, kwargs):
            pars = fitpar(xs, ys, rcond=1e-6)
            k = 2*pars[0]
            if k==0: q0 = np.nan
            else: q0 = -pars[1]/k
            label = '%s (K=%.0f q0=%.3f)' %(prefix, k/parse_unit(self.kunit), q0/parse_unit(self.qunit))
            kwargs['label'] = label
            ax.plot(xs/parse_unit(self.qunit), ys/parse_unit(eunit), **kwargs)
        #ai
        data = np.array([ai.energy(pos) for pos in self.coords])
        e_max = max(np.ceil(max(data-min(data))/parse_unit(eunit)), 1.0)
        add_plot(self.values, data-min(data), 'AI ref', {'linestyle': 'none', 'marker': 'o', 'markerfacecolor': 'k', 'markersize': 12})
        #ffrefs
        totff = np.zeros([len(self.coords)], float)
        colors = ['b', 'g', 'm', 'y', 'c']
        for i, ffref in enumerate(ffrefs):
            data = np.array([ffref.energy(pos) for pos in self.coords])
            totff += data
            add_plot(self.values, data-min(data), ffref.name, {'linestyle': ':', 'color': colors[i], 'linewidth': 2.0})
        #residual valence model if given
        if valence is not None:
            for term in valence.iter_terms():
                valence.check_params(term, ['all'])
            fc = valence.get_params(self.term.index, only='fc')
            rv = valence.get_params(self.term.index, only='rv')
            data = np.array([valence.calc_energy(pos) for pos in self.coords]) - 0.5*fc*(self.values-rv)**2
            totff += data
            add_plot(self.values, data-min(data), 'Residual Valence', {'linestyle': '--', 'color': 'r', 'linewidth':2.0})
        else:
            fc = self.fc
            rv = self.rv
        #plot contribution of current term
        data = 0.5*fc*(self.values-rv)**2
        totff += data
        add_plot(self.values, data-min(data), 'PT Term', {'linestyle': '-', 'color': 'r', 'linewidth':2.0})
        add_plot(self.values, totff-min(totff), 'Total FF', {'linestyle': '-', 'color': [0.4,0.4,0.4], 'linewidth':3.0})
        #decorate plot
        ax.set_xlim([(min(self.values)-self.step)/parse_unit(self.qunit), (max(self.values)+self.step)/parse_unit(self.qunit)])
        ax.set_title('%s-%i' %(self.term.basename, self.term.index))
        ax.set_xlabel('%s [%s]' % (self.term.basename.split('/')[0], self.qunit), fontsize=16)
        ax.set_ylabel('Energy [%s]' %eunit, fontsize=16)
        ax.set_ylim([-0.2, e_max])
        ax.grid()
        ax.legend(loc='upper center', fontsize=16)
        fig.set_size_inches([8, 8])
        if fn is None:
            fn = 'trajectory-%s-%i.png' %(self.term.basename.replace('/', '-'),self.term.index)
        fig.savefig(fn)

    def to_xyz(self, fn=None):
        '''
            Method to write the given trajectory to an XYZ file. This method
            assumes that the coords attribute has been assigned.

            **Optional Arguments**

            fn
                a string defining the name of the output file
        '''
        if 'active' in self.__dict__.keys() and not self.active: return
        if fn is None:
            fn = 'trajectory-%s-%i.xyz' %(self.term.basename.replace('/', '-'),self.term.index)
        f = open(fn, 'w')
        xyz = XYZWriter(f, [pt[Z].symbol for Z in self.numbers])
        for frame, coord in enumerate(self.coords):
            xyz.dump('frame %i' %frame, coord)
        f.close()



class RelaxedStrain(object):
    def __init__(self, system, valence):
        '''
            **Arguments**

            system
                an instance of the Yaff System class defining the system

            valence
                an instance of ValenceFF defining the valence force field
        '''
        self.system = system
        self.valence = valence
        self.strains = None

    def prepare(self, do_terms):
        '''
            Method to initialize trajectories and configure everything required
            for the generate method.
        '''
        self.strains = [None,]*len(self.valence.terms)
        trajectories = [None,]*len(self.valence.terms)
        for term in do_terms:
            assert term.kind==0, 'Only harmonic terms supported for pert traj'
            ic = self.valence.iclist.ictab[self.valence.vlist.vtab[term.index]['ic0']]
            kunit, qunit = term.units
            if ic['kind'] in [0]:
                start=ic['value']-0.02*angstrom
                end=ic['value']+0.02*angstrom
                if start<0.0: start = 0.0
            elif ic['kind'] in [2]:
                start=ic['value']-2*deg
                end=ic['value']+2*deg
                if start<0*deg: start = 0.0
                if end>180*deg: end=180*deg
            elif ic['kind']  in [10]:
                start=-0.01*angstrom
                end=0.01*angstrom
            elif ic['kind'] in [11]:
                start=ic['value']-0.01*angstrom**2
                end=ic['value']+0.01*angstrom**2
                if start<0.0: start=0.0
            else:
                raise NotImplementedError
            #collect all other ics
            ics = []
            for term2 in self.valence.terms:
                #if term2.index == term.index: continue #not the current term ic
                if term2.kind == 3: continue #not a cross term ic
                ics.append(term2.ics[0])
            self.strains[term.index] = Strain(self.system, term.ics[0], ics)
            trajectories[term.index] = Trajectory(term, start, end, self.system.numbers, steps=7)
        return trajectories

    def generate(self, trajectory, remove_com=True):
        '''
            Method to calculate the perturbation trajectory, i.e. the trajectory
            that scans the geometry along the direction of the ic figuring in
            the term with the given index (should be a diagonal term). This
            method should be implemented in the derived classes.

            **Arguments**

            trajectory
                instance of Trajectory class representing the perturbation
                trajectory

            **Optional Arguments**

            remove_com
                if set to True, removes the center of mass translation from the
                resulting perturbation trajectories [default=True].
        '''
        index = trajectory.term.index
        strain = self.strains[index]
        natom = self.system.natom
        if strain is None:
            log.dump('Strain for term %i (%s) is not initialized, skipping.' %(index, self.valence.terms[index].basename))
            return
        q0 = self.valence.iclist.ictab[self.valence.vlist.vtab[index]['ic0']]['value']
        diag = np.ones([strain.ndof+1], float)
        diag[:strain.ndof] *= 0.1*angstrom
        diag[strain.ndof] *= abs(q0-trajectory.targets[0])
        for iq, target in enumerate(trajectory.targets):
            strain.constrain_target = target
            init = np.zeros([strain.ndof+1], float)
            sol, infodict, ier, mesg = scipy.optimize.fsolve(strain.gradient, init, xtol=1e-3, full_output=True, diag=diag)
            if ier==5:
                #fsolve did not converge, flag this frame for deletion
                log.dump('Trajectory of frame %i (target=%.3f) for %s(atoms=%s) did not converge.' %(iq, target, self.valence.terms[index].basename, trajectory.term.get_atoms()))
                trajectory.targets[iq] = np.nan
                continue
            x = strain.coords0 + sol[:3*natom].reshape((-1,3))
            trajectory.values[iq] = strain.constrain_value
            if remove_com:
                com = (x.T*self.system.masses).sum(axis=1)/self.system.masses.sum()
                for i in xrange(natom):
                    x[i,:] -= com
            trajectory.coords[iq,:,:] = x
        #delete flagged frames
        targets = []
        values = []
        coords = []
        for target, value, coord in zip(trajectory.targets, trajectory.values, trajectory.coords):
            if not np.isnan(target):
                targets.append(target)
                values.append(value)
                coords.append(coord)
        trajectory.targets = np.array(targets)
        trajectory.values = np.array(values)
        trajectory.coords = np.array(coords)
        return trajectory

    def estimate(self, trajectory, ai, ffrefs=[], do_valence=False):
        '''
            Method to estimate the FF parameters for the relevant ic from the
            given perturbation trajectory by fitting a harmonic potential to the
            covalent energy along the trajectory.

            **Arguments**

            trajectory
                a Trajectory instance representing the perturbation trajectory

            ai
                an instance of the Reference representing the ab initio input

            **Optional Arguments**

            ffrefs
                a list of Reference instances representing possible a priori
                determined contributions to the force field (such as eg.
                electrostatics and van der Waals)

            do_valence
                If set to True, the current valence force field (stored in
                self.valence) will be used to compute the valence contribution
        '''
        term = trajectory.term
        index = term.index
        basename = term.basename
        if 'active' in trajectory.__dict__.keys() and not trajectory.active:
            log.dump('Trajectory of %s was deactivated: skipping' %(basename))
            return
        qs = trajectory.values.copy()
        AIs = np.zeros(len(trajectory.coords))
        FFs = np.zeros(len(trajectory.coords))
        RESs = np.zeros(len(trajectory.coords))
        for istep, pos in enumerate(trajectory.coords):
            AIs[istep] = ai.energy(pos)
            for ref in ffrefs:
                FFs[istep] += ref.energy(pos)
        if do_valence:
            fc = self.valence.get_params(index, only='fc')
            rv = self.valence.get_params(index, only='rv')
            for istep, pos in enumerate(trajectory.coords):
                RESs[istep] += self.valence.calc_energy(pos) - 0.5*fc*(qs[istep]-rv)**2
        pars = fitpar(qs, AIs-FFs-RESs-min(AIs-FFs-RESs), rcond=-1)
        if pars[0]!=0.0:
            trajectory.fc = 2.0*pars[0]
            trajectory.rv = -pars[1]/(2.0*pars[0])
        else:
            trajectory.fc = 0.0
            trajectory.rv = qs[len(qs)/2]
            log.dump('force constant of %s is zero: rest value set to middle value' %basename)
        #no negative rest values for all ics except dihedrals and bendcos
        if term.ics[0].kind not in [1,3,4,11]:
            if trajectory.rv<0:
                trajectory.rv = 0.0
                log.dump('rest value of %s was negative: set to zero' %basename)




class Strain(ForceField):
    def __init__(self, system, cons_ic, ics, cart_penalty=1.0*angstrom):
        '''
            A class deriving from the Yaff ForceField class to implement the
            strain of a molecular geometry associated with the term defined by
            term_index.

            **Arguments**

            term_index
                Integer defining the index of the term in valence.terms for
                which the current strain is designed.

            system
                A Yaff System instance containing all system information.

            cons_ic
                An instance of Yaff Internal Coordinate representing the
                constrained term in the strain.

            ics
                A list of Yaff Internal Coordinate instances for which the
                strain needs to be minimized.
            
            cart_penalty
                Magnitude of an extra term added to the strain that penalises
                a deviation of the cartesian coordinates of each atom with
                respect to the equilibrium coordinates. This penalty is equal
                to norm(R-R_eq)**2/(2.0*3*Natoms*cart_penalty**2) and prevents
                global translations, global rotations as well as rotations of
                molecular fragments far from the IC under consideration.
        '''
        self.coords0 = system.pos.copy()
        self.ndof = np.prod(self.coords0.shape)
        self.cart_penalty = cart_penalty
        part = ForcePartValence(system)
        for ic in ics:
            part.add_term(Harmonic(1.0, None, ic))
        #set the rest values to the equilibrium values
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
        self.constrain_target = None
        self.constrain_value = None
        self.value = None

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
        self.value = self.compute(gpos=gstrain)
        #compute constraint gradient
        gconstraint = np.zeros(self.coords0.shape)
        self.constraint.update_pos(self.coords0 + X[:self.ndof].reshape((-1,3)))
        self.constrain_value = self.constraint.compute(gpos=gconstraint) + 1.0
        #construct gradient
        grad[:self.ndof] = gstrain.reshape((-1,)) + X[self.ndof]*gconstraint.reshape((-1,))
        grad[self.ndof] = self.constrain_value - self.constrain_target
        #cartesian penalty, i.e. extra penalty for deviation w.r.t. cartesian equilibrium coords
        grad[:self.ndof] += X[:self.ndof]/(self.ndof*self.cart_penalty**2)
        return grad
