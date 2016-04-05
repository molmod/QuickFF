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
        self.qvals = start + (end-start)/(steps-1)*np.array(range(steps), float)
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
        index_min = np.where(data==min(data))[0][0]
        e_max = max(np.ceil(max(data-min(data))/parse_unit(eunit)), 1.0)
        add_plot(self.qvals, data-data[index_min], 'AI ref', {'linestyle': 'none', 'marker': 'o', 'markerfacecolor': 'k', 'markersize': 12})
        #ffrefs
        totff = np.zeros([len(self.qvals)], float)
        colors = ['b', 'g', 'm', 'y', 'c']
        for i, ffref in enumerate(ffrefs):
            data = np.array([ffref.energy(pos) for pos in self.coords])
            totff += data
            add_plot(self.qvals, data-data[index_min], ffref.name, {'linestyle': ':', 'color': colors[i], 'linewidth': 2.0})
        #complete fitted valence model if given
        if valence is not None:
            for term in valence.iter_terms():
                valence.check_params(term, ['all'])
            data = np.array([valence.calc_energy(pos) for pos in self.coords])
            add_plot(self.qvals, data-data[index_min], 'Fitted Valence', {'linestyle': '-', 'color': 'r', 'linewidth':2.0})
            add_plot(self.qvals, totff+data-totff[index_min]-data[index_min], 'Total FF', {'linestyle': '-', 'color': [0.4,0.4,0.4], 'linewidth':3.0})
        #complete fitted term if valence is not given
        else:
            assert self.fc is not None and self.rv is not None
            start = self.qvals[0]
            end = self.qvals[-1]
            step = float((end-start)/100)
            if step==0:
                start -= 0.05*start
                end += 0.05*end
                step = float((end-start)/100)
            qs = np.arange(start, end+step, step)
            data = 0.5*self.fc*(qs - self.rv)**2
            add_plot(qs, data-min(data), 'Fitted Term', {'linestyle': '-', 'color': 'r', 'linewidth':2.0})
        #decorate plot
        ax.set_title('%s-%i' %(self.term.basename, self.term.index))
        ax.set_xlabel('%s [%s]' % (self.term.basename.split('/')[0], self.qunit), fontsize=16)
        ax.set_ylabel('Energy [%s]' %eunit, fontsize=16)
        ax.set_ylim([-0.5, e_max])
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
                start=ic['value']-0.05*angstrom
                end=ic['value']+0.05*angstrom
                if start<0.0: start = 0.0
            elif ic['kind'] in [2]:
                start=ic['value']-5*deg
                end=ic['value']+5*deg
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
                if term2.index == term.index: continue #not the current term ic
                if term2.kind == 3: continue #not a cross term ic
                ics.append(term2.ics[0])
            self.strains[term.index] = Strain(self.system, term.ics[0], ics)
            trajectories[term.index] = Trajectory(term, start, end, self.system.numbers, steps=9)
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
        if strain is None:
            log.dump('Strain for term %i (%s) is not initialized, skipping.' %(index, self.valence.terms[index].basename))
            return
        for iq, q in enumerate(trajectory.qvals):
            strain.constrain_value = q
            if iq==0:
                init = np.zeros([strain.ndof+1], float)
            else:
                init = sol.copy()
            sol, infodict, ier, mesg = scipy.optimize.fsolve(strain.gradient, init, xtol=1e-9, full_output=True)
            if False:#ier==5:
                #fsolve did not converge, flag this frame for deletion
                trajectory.qvals[iq] = np.nan
                continue
            x = strain.coords0 + sol[:self.system.natom*3].reshape((-1,3))
            if remove_com:
                com = (x.T*self.system.masses).sum(axis=1)/self.system.masses.sum()
                for i in xrange(self.system.natom):
                    x[i,:] -= com
            trajectory.coords[iq,:,:] = x
        #delete flagged frames
        qvals = []
        coords = []
        for q, coord in zip(trajectory.qvals, trajectory.coords):
            if not np.isnan(q):
                qvals.append(q)
                coords.append(coord)
        trajectory.qvals = np.array(qvals)
        trajectory.coords = np.array(coords)
        return trajectory

    def estimate(self, trajectory, trajectories, ai, ffrefs=[], do_valence=False):
        '''
            Method to estimate the FF parameters for the relevant ic from the
            given perturbation trajectory by fitting a harmonic potential to the
            covalent energy along the trajectory.

            **Arguments**

            trajectory
                a Trajectory instance representing the perturbation trajectory

            trajectories
                list of Trajectory instances representing all perturbation
                trajectories

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
            log.dump('WARNING: Trajectory of term %i (%s) was deactivated: skipping' %(index, basename))
            return
        qs = trajectory.qvals.copy()
        AIs = np.zeros(len(trajectory.coords))
        FFs = np.zeros(len(trajectory.coords))
        for istep, pos in enumerate(trajectory.coords):
            eai = ai.energy(pos)
            AIs[istep] += eai
            for ref in ffrefs:
                FFs[istep] += ref.energy(pos)
        if do_valence:
            k = self.valence.get_params(index, only='fc')
            self.valence.set_params(index, fc=0.0)
            for istep, pos in enumerate(trajectory.coords):
                FFs[istep] += self.valence.calc_energy(pos)
            self.valence.set_params(index, fc=k)
        pars = fitpar(qs, AIs-FFs, rcond=1e-6)
        if pars[0]!=0.0:
            trajectory.fc = 2*pars[0]
            trajectory.rv = -pars[1]/(2*pars[0])
        else:
            trajectory.fc = 0.0
            vterm = self.valence.vlist.vtab[index]
            ic = self.valence.iclist.ictab[vterm['ic0']]
            trajectory.rv = ic['value']
            log.dump('WARNING: force constant of %s is zero: rest value set to AI equilibrium' %basename)
        #no negative rest values for all ics except dihedrals and bendcos
        if term.ics[0].kind not in [1,3,4]:
            if trajectory.rv<0:
                trajectory.rv = 0
                log.dump('WARNING: rest value of %s was negative: set to zero' %basename)
        #no bending angles larger than 180*deg
        if term.ics[0].kind in [2]:
            if trajectory.rv>180*deg:
                log.dump('WARNING: rest value of %s exceeds 180 deg: term set to BendCHarm with cos(phi0)=-1, deactivated perturbation trajectory' %basename)
                for other_traj in trajectories:
                    if other_traj is None or not other_traj.term.basename==basename: continue
                    other_traj.rv = None
                    other_traj.fc = None
                    other_traj.active = False
                    self.valence.modify_term(
                        other_traj.term.index,
                        Harmonic, [BendCos(*other_traj.term.get_atoms())],
                        basename.replace('BendAHarm', 'BendCHarm'),
                        ['HC_FC_DIAG'], ['kjmol', 'au']
                    )
                    self.valence.set_params(other_traj.term.index, fc=0.0, rv0=-1.0)



class Strain(ForceField):
    def __init__(self, system, cons_ic, ics):
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
        '''
        self.coords0 = system.pos.copy()
        self.ndof = np.prod(self.coords0.shape)
        part = ForcePartValence(system)
        #construct ordered list of atoms in the constrained ic
        if cons_ic.kind==0:
            cons_ic_atoms = cons_ic.index_pairs[0]
        elif cons_ic.kind in [1,2]:
            cons_ic_atoms = [
                cons_ic.index_pairs[0][1],
                cons_ic.index_pairs[1][0],
                cons_ic.index_pairs[1][1],
            ]
        elif cons_ic.kind in [3,4]:
            cons_ic_atoms = [
                cons_ic.index_pairs[0][1],
                cons_ic.index_pairs[1][0],
                cons_ic.index_pairs[1][1],
                cons_ic.index_pairs[2][1],
            ]
        elif cons_ic.kind in [10,11]:
            cons_ic_atoms = [
                cons_ic.index_pairs[0][0],
                cons_ic.index_pairs[0][1],
                cons_ic.index_pairs[2][0],
                cons_ic.index_pairs[2][1],
            ]
        else:
            raise ValueError('IC of kind %i not supported' %cons_ic.kind)
        #add all bonds, bends and diheds that are not constrained
        for bond in system.iter_bonds():
            if not (list(bond)==list(cons_ic_atoms) or list(bond)==list(cons_ic_atoms[::-1])):
                part.add_term(Harmonic(1.0, None, Bond(*bond)))
        for angle in system.iter_angles():
            if not (list(angle)==list(cons_ic_atoms) or list(angle)==list(cons_ic_atoms[::-1])):
                part.add_term(Harmonic(1.0, None, BendAngle(*angle)))
        for dihed in system.iter_dihedrals():
            if not (list(dihed)==list(cons_ic_atoms) or list(dihed)==list(cons_ic_atoms[::-1])):
                part.add_term(Harmonic(1.0, None, DihedAngle(*dihed)))
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
        econstraint = self.constraint.compute(gpos=gconstraint)
        #construct gradient
        grad[:self.ndof] = gstrain.reshape((-1,)) + X[self.ndof]*gconstraint.reshape((-1,)) #+ 0.01*X[:self.ndof]
        grad[self.ndof] = econstraint + 1.0 - self.constrain_value #self.constraint.parts[0].vlist.vtab[0]['energy'] + 1.0 - self.constrain_value
        return grad
