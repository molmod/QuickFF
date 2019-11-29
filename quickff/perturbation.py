# -*- coding: utf-8 -*-
# QuickFF is a code to quickly derive accurate force fields from ab initio input.
# Copyright (C) 2012 - 2019 Louis Vanduyfhuys <Louis.Vanduyfhuys@UGent.be>
# Steven Vandenbrande <Steven.Vandenbrande@UGent.be>,
# Jelle Wieme <Jelle.Wieme@UGent.be>,
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

from __future__ import absolute_import

from molmod.io.xyz import XYZWriter
from molmod.units import *
from molmod.periodic import periodic as pt

from yaff.system import System
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
    def __init__(self, term, start, end, numbers, nsteps=11):
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

            nsteps
                an integer defining the number of steps in the perturbation
                trajectory. The default value is 11 steps.
        '''
        if not term.kind in [0,2,11,12]:
            raise NotImplementedError('Perturbation trajectory only implemented for Harmonic, Fues, MM3Quartic or MM3Bend terms')
        self.term = term
        self.numbers = numbers
        self.qunit = term.units[1]
        self.kunit = term.units[0]
        self.step = (end-start)/(nsteps-1)
        self.targets = start + (end-start)/(nsteps-1)*np.array(list(range(nsteps)), float)
        self.values = np.zeros(nsteps, float)
        self.coords = np.zeros([nsteps, len(numbers), 3], float)
        self.active = True
        self.fc = None
        self.rv = None

    def plot(self, ai, ffrefs=[], valence=None, fn='default', eunit='kjmol', suffix=''):
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

            suffix
                a string to be added to the filename at the end. Is overwritten
                when fn is specified.
        '''
        import matplotlib.pyplot as pp
        if 'active' in list(self.__dict__.keys()) and not self.active: return
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
        add_plot(self.values, data-min(data), 'AI ref', {'linestyle': 'none', 'marker': 'o', 'markerfacecolor': 'k', 'markersize': 12, 'markeredgecolor': 'k'})
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
            valence.set_params(self.term.index, fc=0.0)
            valence.set_params(self.term.index, rv0=0.0)
            data = np.array([valence.calc_energy(pos) for pos in self.coords]) #- 0.5*fc*(self.values-rv)**2
            valence.set_params(self.term.index, fc=fc)
            valence.set_params(self.term.index, rv0=rv)
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
        if fn is 'show':
            pp.show()
        else:
            if fn is 'default':
                fn = 'trajectory-%s-%i%s.png' %(self.term.basename.replace('/', '-'),self.term.index,suffix)
            fig.savefig(fn)
        pp.close()

    def to_xyz(self, fn=None):
        '''
            Method to write the given trajectory to an XYZ file. This method
            assumes that the coords attribute has been assigned.

            **Optional Arguments**

            fn
                a string defining the name of the output file
        '''
        if 'active' in list(self.__dict__.keys()) and not self.active: return
        if fn is None:
            fn = 'trajectory-%s-%i.xyz' %(self.term.basename.replace('/', '-'),self.term.index)
        f = open(fn, 'w')
        xyz = XYZWriter(f, [pt[Z].symbol for Z in self.numbers])
        for frame, coord in enumerate(self.coords):
            xyz.dump('frame %i' %frame, coord)
        f.close()


class RelaxedStrain(object):
    def __init__(self, system, valence, settings):
        '''
            **Arguments**

            system
                a Yaff `System` object defining the system

            valence
                an instance of ValenceFF defining the valence force field

            settings
                a `Settings` instance defining all QuickFF settings
        '''
        self.system0 = system
        self.system_rvecs = system.cell.rvecs.copy()
        self.valence = valence
        self.settings = settings

    def _get_system_copy(self):
        'Routine to get a copy of the equilibrium system'
        numbers = self.system0.numbers.copy()
        coords = self.system0.pos.copy()
        rvecs = self.system_rvecs.copy()
        masses = self.system0.masses.copy()
        bonds = self.system0.bonds.copy()
        ffatypes = self.system0.ffatypes.copy()
        ffatype_ids = self.system0.ffatype_ids.copy()
        return System(
            numbers, coords, rvecs=rvecs,
            ffatypes=ffatypes, ffatype_ids=ffatype_ids,
            masses=masses, bonds=bonds
        )

    system = property(_get_system_copy)

    def prepare(self, do_terms):
        '''
            Method to initialize trajectories and configure everything required
            for the generate method.
        '''
        trajectories = []
        #TODO: make settings options out of the range of the ics in the trajectories as well as the number of steps
        for term in do_terms:
            assert term.kind in [0,2,11,12], 'Only Harmonic, Fues, MM3Quartic or MM3Bend terms supported for pert traj, got term.kind=%i' %term.kind
            ic = self.valence.iclist.ictab[self.valence.vlist.vtab[term.index]['ic0']]
            kunit, qunit = term.units
            if ic['kind'] in [0]:
                start=ic['value']-0.02*angstrom
                end=ic['value']+0.02*angstrom
                if start<0.0: start = 0.0
            elif ic['kind'] in [2,4]:
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
            trajectories.append(Trajectory(term, start, end, self.system0.numbers.copy(), nsteps=7))
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
        #TODO: find out why system.cell is not parsed correctly when using scoop
        #force correct rvecs
        self.system0.cell.update_rvecs(self.system_rvecs)
        with log.section('PTGEN', 4, timer='PT Generate'):
            log.dump('  Generating %s(atoms=%s)' %(trajectory.term.basename, trajectory.term.get_atoms()))
            strain = Strain(self.system, trajectory.term, self.valence.terms)
            natom = self.system0.natom
            q0 = self.valence.iclist.ictab[self.valence.vlist.vtab[trajectory.term.index]['ic0']]['value']
            diag = np.array([0.1*angstrom,]*3*natom+[abs(q0-trajectory.targets[0])])
            sol = None
            for iq, target in enumerate(trajectory.targets):
                log.dump('    Frame %i (target=%.3f)' %(iq, target))
                strain.constrain_target = target
                if abs(target-q0)<1e-6:
                    sol = np.zeros([strain.ndof+1],float)
                    #call strain.gradient once to compute/store/log relevant information
                    strain.gradient(sol)
                else:
                    if sol is not None:
                        init = sol.copy()
                    else:
                        init = np.zeros([3*natom+1], float)
                    init[-1] = np.sign(q0-target)
                    sol, infodict, ier, mesg = scipy.optimize.fsolve(strain.gradient, init, xtol=self.settings.pert_traj_tol, full_output=True, diag=diag)
                    if ier!=1:
                        #fsolve did not converge, try again after adding small random noise
                        log.dump('      %s' %mesg.replace('\n', ' '))
                        log.dump('    Frame %i (target=%.3f) %s(%s) did not converge. Trying again with slightly perturbed initial conditions.' %(
                            iq, target, trajectory.term.basename, trajectory.term.get_atoms()
                        ))
                        #try one more time
                        init = sol.copy()
                        init[:3*natom] += np.random.normal(0.0, 0.01, [3*natom])*angstrom
                        sol, infodict, ier, mesg = scipy.optimize.fsolve(strain.gradient, init, xtol=self.settings.pert_traj_tol, full_output=True, diag=diag)
                        #fsolve did STILL not converge, flag this frame for deletion
                        if ier!=1:
                            log.dump('      %s' %mesg.replace('\n', ' '))
                            log.dump('    Frame %i (target=%.3f) %s(%s) STILL did not converge.' %(
                                iq, target, trajectory.term.basename, trajectory.term.get_atoms()
                            ))
                            trajectory.targets[iq] = np.nan
                            continue
                x = self.system0.pos.copy() + sol[:3*natom].reshape((-1,3))
                trajectory.values[iq] = strain.constrain_value
                log.dump('    Converged (value=%.3f, lagmult=%.3e)' %(strain.constrain_value,sol[3*natom]))
                if remove_com:
                    com = (x.T*self.system0.masses.copy()).sum(axis=1)/self.system0.masses.sum()
                    for i in range(natom):
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

    def estimate(self, trajectory, ai, ffrefs=[], do_valence=False, energy_noise=None, Nerrorsteps=100):
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

            energy_noise
                If set to a float, the parabolic fitting will be repeated
                Nerrorsteps times including normal noise on top of the reference
                value. The mean of the noise is 0, while the std equals the
                number given by energy_noise. The resulting fits give a
                distribution of force  constants and rest values instead of
                single value, the std is used to identify bad estimates, the
                mean is used for the actual FF parametrs. If set to nan, the
                parabolic fit is performed only once without any noise.
        '''
        with log.section('PTEST', 3, timer='PT Estimate'):
            term = trajectory.term
            index = term.index
            basename = term.basename
            if 'active' in list(trajectory.__dict__.keys()) and not trajectory.active:
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
                self.valence.set_params(index, fc=0.0)
                self.valence.set_params(index, rv0=0.0)
                for istep, pos in enumerate(trajectory.coords):
                    RESs[istep] += self.valence.calc_energy(pos) #- 0.5*fc*(qs[istep]-rv)**2
                self.valence.set_params(index, fc=fc)
                self.valence.set_params(index, rv0=rv)
            pars = fitpar(qs, AIs-FFs-RESs-min(AIs-FFs-RESs), rcond=-1)
            if energy_noise is None:
                if pars[0]!=0.0:
                    trajectory.fc = 2.0*pars[0]
                    trajectory.rv = -pars[1]/(2.0*pars[0])
                else:
                    trajectory.fc = 0.0
                    trajectory.rv = qs[len(qs)//2]
                    log.dump('force constant of %s is zero: rest value set to middle value' %basename)
            else:
                with log.section('PTEST', 4, timer='PT Estimate'):
                    log.dump('Performing noise analysis for trajectory of %s' %basename)
                    As = [pars[0]]
                    Bs = [pars[1]]
                    for i in range(Nerrorsteps):
                        pars = fitpar(qs, AIs-FFs-RESs-min(AIs-FFs-RESs)+np.random.normal(0.0, energy_noise, size=AIs.shape), rcond=-1)
                        As.append(pars[0])
                        Bs.append(pars[1])
                    if 0.0 in As:
                        log.dump('  force constant of zero detected, removing the relevant runs from analysis')
                    Bs = np.array([b for a,b in zip(As,Bs) if a!=0.0])
                    As = np.array([a for a in As if a!=0.0])
                    ks = As*2.0
                    q0s = -Bs/(2.0*As)
                    kunit = trajectory.term.units[0]
                    qunit = trajectory.term.units[1]
                    log.dump('    k  = %8.3f +- %6.3f (noisefree: %8.3f) %s' %(ks.mean()/parse_unit(kunit), ks.std()/parse_unit(kunit), ks[0]/parse_unit(kunit), kunit))
                    log.dump('    q0 = %8.3f +- %6.3f (noisefree: %8.3f) %s' %(q0s.mean()/parse_unit(qunit), q0s.std()/parse_unit(qunit), q0s[0]/parse_unit(qunit), qunit))
                    if q0s.std()/q0s.mean()>0.01:
                        with log.section('PTEST', 3, timer='PT Estimate'):
                            fc, rv = self.valence.get_params(trajectory.term.index)
                            if rv is None:
                                log.dump('Noise on rest value of %s to high, using ab initio rest value' %basename)
                                pars = fitpar(qs, AIs-FFs-RESs-min(AIs-FFs-RESs)+np.random.normal(0.0, energy_noise, size=AIs.shape), rcond=-1)
                                if pars[0]!=0.0:
                                    trajectory.fc = 2.0*pars[0]
                                    trajectory.rv = -pars[1]/(2.0*pars[0])
                                else:
                                    trajectory.fc = 0.0
                                    trajectory.rv = qs[len(qs)//2]
                                    log.dump('AI force constant of %s is zero: rest value set to middle value' %basename)
                            else:
                                log.dump('Noise on rest value of %s to high, using previous value' %basename)
                                trajectory.fc = fc
                                trajectory.rv = rv
                    else:
                        trajectory.fc = ks.mean()
                        trajectory.rv = q0s.mean()
            #no negative rest values for all ics except dihedrals and bendcos
            if term.ics[0].kind not in [1,3,4,11]:
                if trajectory.rv<0:
                    trajectory.rv = 0.0
                    log.dump('rest value of %s was negative: set to zero' %basename)




class Strain(ForceField):
    def __init__(self, system, term, other_terms, cart_penalty=1e-3*angstrom):
        '''
            A class deriving from the Yaff ForceField class to implement the
            strain of a molecular geometry associated with the term defined by
            term_index.

            **Arguments**

            system
                A Yaff System instance containing all system information.

            term
                a Term instance representing the term of the perturbation
                trajectory of the current strain

            other_terms
                a list of Term instances representing all other terms for ICs
                for which a strain contribution should be added

            **Keyword Arguments**

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
        self.cons_ic_atindexes = term.get_atoms()
        #construct main strain
        strain = ForcePartValence(system)
        for other in other_terms:
            if other.kind == 3: continue #no cross terms
            strain.add_term(Harmonic(1.0, None, other.ics[0]))
        #set the rest values to the equilibrium values
        strain.dlist.forward()
        strain.iclist.forward()
        for iterm in range(strain.vlist.nv):
            vterm = strain.vlist.vtab[iterm]
            ic = strain.iclist.ictab[vterm['ic0']]
            vterm['par1'] = ic['value']
        ForceField.__init__(self, system, [strain])
        #Abuse the Chebychev1 polynomial to simply get the value of q-1 and
        #implement the contraint
        constraint = ForcePartValence(system)
        constraint.add_term(Chebychev1(-2.0,term.ics[0]))
        self.constraint = ForceField(system, [constraint])
        self.constrain_target = None
        self.constrain_value = None
        self.value = None

    def gradient(self, X):
        '''
            Compute the gradient of the strain w.r.t. Cartesian coordinates of
            the system. For the ic that needs to be constrained, a Lagrange
            multiplier is included.
        '''
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
        indices = np.array([[3*i,3*i+1,3*i+2] for i in range(self.ndof//3) if i not in self.cons_ic_atindexes]).ravel()
        if len(indices)>0:
            grad[indices] += X[indices]/(self.ndof*self.cart_penalty**2)
        with log.section('PTGEN', 4, timer='PT Generate'):
            log.dump('      Gradient:  rms = %.3e  max = %.3e  cnstr = %.3e' %(np.sqrt((grad[:self.ndof]**2).mean()), max(grad[:self.ndof]), grad[self.ndof]))
        return grad
