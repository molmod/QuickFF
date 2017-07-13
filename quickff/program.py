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
from molmod.units import *

from quickff.valence import ValenceFF
from quickff.perturbation import RelaxedStrain
from quickff.cost import HessianFCCost
from quickff.paracontext import paracontext
from quickff.io import dump_charmm22_prm, dump_charmm22_psf, dump_yaff
from quickff.log import log
from quickff.tools import chebychev

from yaff.system import System
from yaff.pes.vlist import Cosine, Harmonic, Chebychev1, Chebychev4
from yaff.pes.iclist import BendAngle, BendCos, OopDist

import os, cPickle, numpy as np, datetime

__all__ = [
    'BaseProgram', 'MakeTrajectories', 'PlotTrajectories', 'DeriveFF',
]

class BaseProgram(object):
    '''
        Base program which implements all possible steps of a force field
        fitting program. The actual sequence of the steps are defined in the
        deriving classes.
    '''
    def __init__(self, system, ai, settings, ffrefs=[]):
        '''
            **Arguments**

            system
                a Yaff `System` instance defining the system

            ai
                a `Reference` instance corresponding to the ab initio input data

            settings
                a `Settings` instance defining all QuickFF settings

            **Optional Arguments**

            ffrefs
                a list of `Reference` instances defining the a-priori force
                field contributions.
        '''
        with log.section('PROG', 2, timer='Initializing'):
            log.dump('Initializing program')
            self.settings = settings
            self.system = system
            self.ai = ai
            self.ffrefs = ffrefs
            self.valence = ValenceFF(system, settings)
            self.perturbation = RelaxedStrain(system, self.valence, settings)
            self.trajectories = None

    def reset_system(self):
        '''
            routine to reset the system coords to the ai equilbrium
        '''
        log.dump('Resetting system coordinates to ab initio ref')
        self.system.pos = self.ai.coords0.copy()
        self.valence.dlist.forward()
        self.valence.iclist.forward()


    def update_trajectory_terms(self):
        '''
            Routine to make ``self.valence.terms`` and the term attribute of each
            trajectory in ``self.trajectories`` consistent again. This is usefull
            if the trajectory were read from a file and the ``valenceFF`` instance
            was modified.
        '''
        log.dump('Updating terms of trajectories to current valenceFF terms')
        with log.section('PTUPD', 3):
            #update the terms in the trajectories to match the terms in
            #self.valence
            for traj in self.trajectories:
                found = False
                for term in self.valence.iter_terms():
                    if traj.term.get_atoms()==term.get_atoms():
                        if found: raise ValueError('Found two terms for trajectory %s with atom indices %s' %(traj.term.basename, str(traj.term.get_atoms())))
                        traj.term = term
                        if 'PT_ALL' not in term.tasks:
                            log.dump('PT_ALL not in tasks of %s-%i, deactivated PT' %(term.basename, term.index))
                            traj.active = False
                        found = True
                if not found:
                    log.warning('No term found for trajectory %s with atom indices %s, deactivating trajectory' %(traj.term.basename, str(traj.term.get_atoms())))
                    traj.active = False
            #check if every term with task PT_ALL has a trajectory associated
            #with it. It a trajectory is missing, generate it.
            for term in self.valence.iter_terms():
                if 'PT_ALL' not in term.tasks: continue
                found = False
                for traj in self.trajectories:
                    if term.get_atoms()==traj.term.get_atoms():
                        if found: raise ValueError('Found two trajectories for term %s with atom indices %s' %(term.basename, str(term.get_atoms())))
                        found =True
                if not found:
                    log.warning('No trajectory found for term %s with atom indices %s. Generating it now.' %(term.basename, str(term.get_atoms())))
                    trajectory = self.perturbation.prepare([term])[term.index]
                    self.perturbation.generate(trajectory)
                    self.trajectories.append(trajectory)

    def average_pars(self):
        '''
            Average force field parameters over master and slaves.
        '''
        log.dump('Averaging force field parameters over master and slaves')
        for master in self.valence.iter_masters():
            npars = len(self.valence.get_params(master.index))
            pars = np.zeros([len(master.slaves)+1, npars], float)
            pars[0,:] = np.array(self.valence.get_params(master.index))
            for i, islave in enumerate(master.slaves):
                pars[1+i,:] = np.array(self.valence.get_params(islave))
            if master.kind in [0,2,11,12]:#harmonic,fues,MM3Quartic,MM3Bend
                fc, rv = pars.mean(axis=0)
                self.valence.set_params(master.index, fc=fc, rv0=rv)
                for islave in master.slaves:
                    self.valence.set_params(islave, fc=fc, rv0=rv)
            elif master.kind==1:
                a0, a1, a2, a3 = pars.mean(axis=0)
                self.valence.set_params(master.index, a0=a0, a1=a1, a2=a2, a3=a3)
                for islave in master.slaves:
                    self.valence.set_params(islave, a0=a0, a1=a1, a2=a2, a3=a3)
            elif master.kind==3:#cross
                fc, rv0, rv1 = pars.mean(axis=0)
                self.valence.set_params(master.index, fc=fc, rv0=rv0, rv1=rv1)
                for islave in master.slaves:
                    self.valence.set_params(islave, fc=fc, rv0=rv0, rv1=rv1)
            elif master.kind==4:#cosine
                assert pars[:,0].std()<1e-6, 'dihedral multiplicity not unique'
                m, fc, rv = pars.mean(axis=0)
                self.valence.set_params(master.index, fc=fc, rv0=rv, m=m)
                for islave in master.slaves:
                    self.valence.set_params(islave, fc=fc, rv0=rv, m=m)
            elif master.kind in [5, 6, 7, 8, 9]:#chebychev
                assert pars.shape[1]==2
                fc = pars[:,0].mean()
                self.valence.set_params(master.index, fc=fc)
                for islave in master.slaves:
                    self.valence.set_params(islave, fc=fc)
            else:
                raise NotImplementedError

    def make_output(self):
        '''
            Dump Yaff parameters, Yaff system, plot energy contributions along
            perturbation trajectories and dump perturbation trajectories to XYZ
            files.
        '''
        if self.settings.fn_yaff is not None:
            dump_yaff(self.valence, self.settings.fn_yaff)
        if self.settings.fn_charmm22_prm is not None:
            dump_charmm22_prm(self.valence, self.settings.fn_charmm22_prm)
        if self.settings.fn_charmm22_psf is not None:
            dump_charmm22_psf(self.system, self.valence, self.settings.fn_charmm22_psf)
        if self.settings.fn_sys is not None:
            self.system.to_file(self.settings.fn_sys)
        self.plot_trajectories(do_valence=True)

    def plot_trajectories(self, do_valence=False):
        '''
            Plot energy contributions along perturbation trajectories and dump
            perturbation trajectories to XYZ files.
        '''
        only = self.settings.only_traj
        if not isinstance(only, list): only = [only]
        with log.section('PLOT', 3, timer='PT plot energy'):
            if self.settings.plot_traj:
                valence = None
                if do_valence: valence=self.valence
                for trajectory in self.trajectories:
                    if trajectory is None: continue
                    for pattern in only:
                        if pattern=='PT_ALL' or pattern in trajectory.term.basename:
                            trajectory.plot(self.ai, ffrefs=self.ffrefs, valence=valence)
        with log.section('XYZ', 3, timer='PT dump XYZ'):
            if self.settings.xyz_traj:
                for trajectory in self.trajectories:
                    if trajectory is None: continue
                    for pattern in only:
                        if pattern=='PT_ALL' or pattern in trajectory.term.basename:
                            trajectory.to_xyz()

    def do_pt_generate(self):
        'Generate perturbation trajectories.'
        with log.section('PTGEN', 2, timer='PT Generate'):
            #read if an existing file was specified through fn_traj
            fn_traj = self.settings.fn_traj
            if fn_traj is not None and os.path.isfile(fn_traj):
                self.trajectories = cPickle.load(open(fn_traj, 'r'))
                log.dump('Trajectories read from file %s' %fn_traj)
                self.update_trajectory_terms()
                newname = 'updated_'+fn_traj.split('/')[-1]
                cPickle.dump(self.trajectories, open(newname, 'w'))
                return
            #configure
            self.reset_system()
            only = self.settings.only_traj
            if isinstance(only, str):
                do_terms = [term for term in self.valence.terms if only in term.tasks]
            else:
                do_terms = []
                for pattern in only:
                    for term in self.valence.iter_terms(pattern):
                        do_terms.append(term)
            trajectories = self.perturbation.prepare(do_terms)
            #compute
            log.dump('Constructing trajectories')
            self.trajectories = paracontext.map(self.perturbation.generate, [traj for traj in trajectories if traj.active])
            #write the trajectories to the non-existing file fn_traj
            if fn_traj is not None:
                assert not os.path.isfile(fn_traj)
                cPickle.dump(self.trajectories, open(fn_traj, 'w'))
                log.dump('Trajectories stored to file %s' %fn_traj)

    def do_pt_estimate(self, do_valence=False, logger_level=3):
        '''
            Estimate force constants and rest values from the perturbation
            trajectories

            **Optional Arguments**

            do_valence
                if set to True, the current valence force field will be used to
                estimate the contribution of all other valence terms.
        '''
        with log.section('PTEST', 2, timer='PT Estimate'):
            self.reset_system()
            message = 'Estimating FF parameters from perturbation trajectories'
            if do_valence: message += ' with valence reference'
            log.dump(message)
            #compute fc and rv from trajectory
            only = self.settings.only_traj
            for traj in self.trajectories:
                if traj is None: continue
                if only is not None:
                    basename = self.valence.terms[traj.term.master].basename
                    if isinstance(only, list) and basename not in only:
                        continue
                    elif only!=basename:
                        continue
                self.perturbation.estimate(traj, self.ai, ffrefs=self.ffrefs, do_valence=do_valence)
            #set force field parameters to computed fc and rv
            for traj in self.trajectories:
                if traj is None: continue
                if only is not None:
                    basename = self.valence.terms[traj.term.master].basename
                    if isinstance(only, list) and basename not in only:
                        continue
                    elif only!=basename:
                        continue
                self.valence.set_params(traj.term.index, fc=traj.fc, rv0=traj.rv)
            #output
            self.valence.dump_logger(print_level=logger_level)
            #do not add average here since the fluctuation on the parameters is
            #required for do_pt_postprocess. Average will be done at the end of
            #do_pt_postprocess

    def do_pt_postprocess(self):
        '''
            Do some first post processing of the ff parameters estimated from
            the perturbation trajectories including:

                * detecting bend patterns with rest values of 90 and 180 deg
                * detecting bend patterns with rest values only close to 180 deg
                * transforming SqOopDist with rv=0.0 to OopDist
                * averaging parameters
        '''
        with log.section('PTPOST', 2, timer='PT Post process'):
            if self.settings.do_squarebend:
                self.do_squarebend()
            if self.settings.do_bendclin:
                self.do_bendclin()
            if self.settings.do_sqoopdist_to_oopdist:
                self.do_sqoopdist_to_oopdist()
            self.average_pars()

    def do_eq_setrv(self, tasks, logger_level=3):
        '''
            Set the rest values to their respective AI equilibrium values.
        '''
        with log.section('EQSET', 2, timer='Equil Set RV'):
            self.reset_system()
            log.dump('Setting rest values to AI equilibrium values for tasks %s' %' '.join(tasks))
            for term in self.valence.terms:
                vterm = self.valence.vlist.vtab[term.index]
                if np.array([task in term.tasks for task in tasks]).any():
                    if term.kind==3:#cross term
                        ic0 = self.valence.iclist.ictab[vterm['ic0']]
                        ic1 = self.valence.iclist.ictab[vterm['ic1']]
                        self.valence.set_params(term.index, rv0=ic0['value'], rv1=ic1['value'])
                    elif term.kind==4 and term.ics[0].kind==4:#Cosine of DihedAngle
                        ic = self.valence.iclist.ictab[vterm['ic0']]
                        m = self.valence.get_params(term.index, only='m')
                        rv = ic['value']%(360.0*deg/m)
                        with log.section('EQSET', 4, timer='Equil Set RV'):
                            log.dump('Set rest value of %s(%s) (eq=%.3f deg) to %.3f deg' %(
                                term.basename,
                                '.'.join([str(at) for at in term.get_atoms()]),
                                ic['value']/deg, rv/deg
                            ))
                        self.valence.set_params(term.index, rv0=rv)
                    else:
                        rv = self.valence.iclist.ictab[vterm['ic0']]['value']
                        self.valence.set_params(term.index, rv0=rv)
            self.valence.dump_logger(print_level=logger_level)
            self.average_pars()

    def do_hc_estimatefc(self, tasks, logger_level=3, do_svd=False, do_mass_weighting=True):
        '''
            Refine force constants using Hessian Cost function.

            **Arguments**

            tasks
                A list of strings identifying which terms should have their
                force constant estimated from the hessian cost function. Using
                such a flag, one can distinguish between for example force
                constant refinement (flag=HC_FC_DIAG) of the diagonal terms and
                force constant estimation of the cross terms (flag=HC_FC_CROSS).
                If the string 'all' is present in tasks, all fc's will be
                estimated.

            **Optional Arguments**

            logger_level
                print level at which the resulting parameters should be dumped to
                the logger. By default, the parameters will only be dumped at
                the highest log level.

            do_svd
                whether or not to do an SVD decomposition before solving the
                set of equations and explicitly throw out the degrees of
                freedom that correspond to the lowest singular values.

            do_mass_weighting
                whether or not to apply mass weighing to the ab initio hessian
                and the force field contributions before doing the fitting.
        '''
        with log.section('HCEST', 2, timer='HC Estimate FC'):
            self.reset_system()
            log.dump('Estimating force constants from Hessian cost for tasks %s' %' '.join(tasks))
            term_indices = []
            for index in xrange(self.valence.vlist.nv):
                term = self.valence.terms[index]
                flagged = False
                for flag in tasks:
                    if flag in term.tasks:
                        flagged = True
                        break
                if flagged:
                    #first check if all rest values and multiplicities have been defined
                    if term.kind==0: self.valence.check_params(term, ['rv'])
                    if term.kind==1: self.valence.check_params(term, ['a0', 'a1', 'a2', 'a3'])
                    if term.kind==3: self.valence.check_params(term, ['rv0','rv1'])
                    if term.kind==4: self.valence.check_params(term, ['rv', 'm'])
                    if term.is_master():
                        term_indices.append(index)
                else:
                    #first check if all pars have been defined
                    if term.kind==0: self.valence.check_params(term, ['fc', 'rv'])
                    if term.kind==1: self.valence.check_params(term, ['a0', 'a1', 'a2', 'a3'])
                    if term.kind==3: self.valence.check_params(term, ['fc', 'rv0','rv1'])
                    if term.kind==4: self.valence.check_params(term, ['fc', 'rv', 'm'])
            if len(term_indices)==0:
                log.dump('No terms (with task in %s) found to estimate FC from HC' %(str(tasks)))
                return
            cost = HessianFCCost(self.system, self.ai, self.valence, term_indices, ffrefs=self.ffrefs, do_mass_weighting=do_mass_weighting)
            fcs = cost.estimate(do_svd=do_svd)
            for index, fc in zip(term_indices, fcs):
                master = self.valence.terms[index]
                assert master.is_master()
                self.valence.set_params(index, fc=fc)
                for islave in master.slaves:
                    self.valence.set_params(islave, fc=fc)
            self.valence.dump_logger(print_level=logger_level)

    def do_cross_init(self):
        '''
            Set the rest values of cross terms to the rest values of the
            corresponding diagonal terms. The force constants are initialized
            to zero.
        '''
        with log.section('VAL', 2, 'Initializing'):
            self.reset_system()
            self.valence.init_cross_angle_terms()

            #function to find rest value
            def find_rest_value(label):
                candidates = [cand for cand in self.valence.iter_masters(label=label, use_re=True)]
                assert len(candidates)<2, 'Multiple masters found for %s: %s' %(label, ','.join([cand.basename for cand in candidates]))
                if len(candidates)==0:
                    if label.startswith('^Bond') or label.startswith('^Bend') or label.startswith('^Tors'):
                        sublabels = label.split('|')
                        prefix0, types0 = sublabels[0].lstrip('^').rstrip('$').split('/')
                        label  = '^'+prefix0+'/'+'\.'.join(types0.split('.')[::-1])+'$'
                        if len(sublabels)>1:
                            prefix1, types1 = sublabels[1].lstrip('^').rstrip('$').split('/')
                            label += '|'
                            label += '^'+prefix1+'/'+'\.'.join(types1.split('.')[::-1])+'$'
                        candidates = [cand for cand in self.valence.iter_masters(label=label, use_re=True)]
                        assert len(candidates)<2, 'Multiple masters found for %s: %s' %(label, ','.join([cand.basename for cand in candidates]))
                        if len(candidates)==0:
                            return None
                can = candidates[0]
                if can.basename.startswith('TorsCheby') or can.basename.startswith('BendCheby'):
                    return -self.valence.get_params(can.index, only='sign')
                else:
                    return self.valence.get_params(can.index, only='rv')

            #set rest values and initialize fc for bond-bond cross
            for term in self.valence.iter_masters('^Cross/.*/bb$', use_re=True):
                types = term.basename.split('/')[1].split('.')
                assert len(types)==3, 'Found angle cross terms with more/less than 3 atom types'
                rv0 = find_rest_value('^Bond.*/%s$' %('.'.join(types[:2])))
                rv1 = find_rest_value('^Bond.*/%s$' %('.'.join(types[1:])))
                assert rv0 is not None, 'Rest value of BondHarm/%s not found' %('.'.join(types[:2]))
                assert rv1 is not None, 'Rest value of BondHarm/%s not found' %('.'.join(types[1:]))
                self.valence.set_params(term.index, fc=0.0, rv0=rv0, rv1=rv1)
                for index in term.slaves: self.valence.set_params(index, fc=0.0, rv0=rv0, rv1=rv1)

            #set rest values and initialize fc for bond-angle cross
            for term in self.valence.iter_masters('^Cross/.*/b0a$', use_re=True):
                types = term.basename.split('/')[1].split('.')
                assert len(types)==3, 'Found angle cross terms with more/less than 3 atom types'
                rv0 = find_rest_value('^Bond.*/%s$' %('.'.join(types[:2])))
                rv1 = find_rest_value('^BendAHarm/%s$|^BendMM3/%s$' %('.'.join(types),'.'.join(types)))
                assert rv0 is not None, 'Rest value of BondHarm/%s not found' %('.'.join(types[:2]))
                assert rv1 is not None, 'Rest value of BendAHarm|BendMM3/%s not found' %('.'.join(types))
                self.valence.set_params(term.index, fc=0.0, rv0=rv0, rv1=rv1)
                for index in term.slaves: self.valence.set_params(index, fc=0.0, rv0=rv0, rv1=rv1)

            for term in self.valence.iter_masters('^Cross/.*/b1a$', use_re=True):
                types = term.basename.split('/')[1].split('.')
                assert len(types)==3, 'Found angle cross terms with more/less than 3 atom types'
                rv0 = find_rest_value('^Bond.*/%s$' %('.'.join(types[1:])))
                rv1 = find_rest_value('^BendAHarm/%s$|^BendMM3/%s$' %('.'.join(types),'.'.join(types)))
                assert rv0 is not None, 'Rest value of BondHarm/%s not found' %('.'.join(types[1:]))
                assert rv1 is not None, 'Rest value of BendAHarm|BendMM3/%s not found' %('.'.join(types))
                self.valence.set_params(term.index, fc=0.0, rv0=rv0, rv1=rv1)
                for index in term.slaves: self.valence.set_params(index, fc=0.0, rv0=rv0, rv1=rv1)

            if self.settings.do_cross_DSS or self.settings.do_cross_DSD or self.settings.do_cross_DAD or self.settings.do_cross_DAA:
                self.valence.init_cross_dihed_terms()
                #set rest values and initialize fc for bond-bond cross
                for m in [1,2,3,4,6]:
                    for term in self.valence.iter_masters('^CrossBondDih%i/.*/bb$' %m, use_re=True):
                        types = term.basename.split('/')[1].split('.')
                        assert len(types)==4, 'Found angle cross terms with more/less than 4 atom types'
                        rv0 = find_rest_value('^Bond.*/%s$' %('.'.join(types[:2])))
                        rv1 = find_rest_value('^Bond.*/%s$' %('.'.join(types[2:])))
                        self.valence.set_params(term.index, fc=0.0, rv0=rv0, rv1=rv1)
                        for index in term.slaves: self.valence.set_params(index, fc=0.0, rv0=rv0, rv1=rv1)

                #set rest values and initialize fc for bond-dihed cross
                for m in [1,2,3,4,6]:
                    for term in self.valence.iter_masters('^CrossBondDih%i/.*/b0d$' %m, use_re=True):
                        types = term.basename.split('/')[1].split('.')
                        assert len(types)==4, 'Found angle cross terms with more/less than 4 atom types'
                        rv0 = find_rest_value('^Bond.*/%s$' %('.'.join(types[:2])))
                        rv1 = find_rest_value('^TorsCheby%i/%s$' %(m,'.'.join(types)))
                        self.valence.set_params(term.index, fc=0.0, rv0=rv0, rv1=rv1)
                        for index in term.slaves: self.valence.set_params(index, fc=0.0, rv0=rv0, rv1=rv1)

                    for term in self.valence.iter_masters('^CrossBondDih%i/.*/b1d$' %m, use_re=True):
                        types = term.basename.split('/')[1].split('.')
                        assert len(types)==4, 'Found angle cross terms with more/less than 4 atom types'
                        rv0 = find_rest_value('^Bond.*/%s$' %('.'.join(types[1:3])))
                        rv1 = find_rest_value('^TorsCheby%i/%s$' %(m,'.'.join(types)))
                        self.valence.set_params(term.index, fc=0.0, rv0=rv0, rv1=rv1)
                        for index in term.slaves: self.valence.set_params(index, fc=0.0, rv0=rv0, rv1=rv1)

                    for term in self.valence.iter_masters('^CrossBondDih%i/.*/b2d$' %m, use_re=True):
                        types = term.basename.split('/')[1].split('.')
                        assert len(types)==4, 'Found angle cross terms with more/less than 4 atom types'
                        rv0 = find_rest_value('^Bond.*/%s$' %('.'.join(types[2:])))
                        rv1 = find_rest_value('^TorsCheby%i/%s$' %(m,'.'.join(types)))
                        self.valence.set_params(term.index, fc=0.0, rv0=rv0, rv1=rv1)
                        for index in term.slaves: self.valence.set_params(index, fc=0.0, rv0=rv0, rv1=rv1)

                #set rest values and initialize fc for angle-angle cross
                for m in [1,2,3,4,6]:
                    for term in self.valence.iter_masters('^CrossBendDih%i/.*/aa$' %m, use_re=True):
                        types = term.basename.split('/')[1].split('.')
                        assert len(types)==4, 'Found angle cross terms with more/less than 4 atom types'
                        rv0 = find_rest_value('^BendAHarm/%s$|^BendMM3/%s$' %('.'.join(types[:3]),'.'.join(types[:3])))
                        rv1 = find_rest_value('^BendAHarm/%s$|^BendMM3/%s$' %('.'.join(types[1:]),'.'.join(types[1:])))
                        self.valence.set_params(term.index, fc=0.0, rv0=rv0, rv1=rv1)
                        for index in term.slaves: self.valence.set_params(index, fc=0.0, rv0=rv0, rv1=rv1)

                #set rest values and initialize fc for angle-dihed cross
                for m in [1,2,3,4,6]:
                    for term in self.valence.iter_masters('^CrossBendDih%i/.*/a0d$' %m, use_re=True):
                        types = term.basename.split('/')[1].split('.')
                        assert len(types)==4, 'Found angle cross terms with more/less than 4 atom types'
                        rv0 = find_rest_value('^BendAHarm/%s$|^BendMM3/%s$' %('.'.join(types[:3]),'.'.join(types[:3])))
                        rv1 = find_rest_value('^TorsCheby%i/%s$' %(m,'.'.join(types)))
                        self.valence.set_params(term.index, fc=0.0, rv0=rv0, rv1=rv1)
                        for index in term.slaves: self.valence.set_params(index, fc=0.0, rv0=rv0, rv1=rv1)

                    for term in self.valence.iter_masters('^CrossBendDih%i/.*/a1d$' %m, use_re=True):
                        types = term.basename.split('/')[1].split('.')
                        assert len(types)==4, 'Found angle cross terms with more/less than 4 atom types'
                        rv0 = find_rest_value('^BendAHarm/%s$|^BendMM3/%s$' %('.'.join(types[1:]),'.'.join(types[1:])))
                        rv1 = find_rest_value('^TorsCheby%i/%s$' %(m,'.'.join(types)))
                        self.valence.set_params(term.index, fc=0.0, rv0=rv0, rv1=rv1)
                        for index in term.slaves: self.valence.set_params(index, fc=0.0, rv0=rv0, rv1=rv1)

                    for term in self.valence.iter_masters('^CrossCBendDih%i/.*/a0d$' %m, use_re=True):
                        types = term.basename.split('/')[1].split('.')
                        assert len(types)==4, 'Found angle cross terms with more/less than 4 atom types'
                        rv0 = find_rest_value('^BendCLin/%s$' %('.'.join(types[:3])))
                        rv1 = find_rest_value('^TorsCheby%i/%s$' %(m,'.'.join(types)))
                        self.valence.set_params(term.index, fc=0.0, rv0=rv0, rv1=rv1)
                        for index in term.slaves: self.valence.set_params(index, fc=0.0, rv0=rv0, rv1=rv1)

                    for term in self.valence.iter_masters('^CrossCBendDih%i/.*/a1d$' %m, use_re=True):
                        types = term.basename.split('/')[1].split('.')
                        assert len(types)==4, 'Found angle cross terms with more/less than 4 atom types'
                        rv0 = find_rest_value('^BendCLin/%s$' %('.'.join(types[1:])))
                        rv1 = find_rest_value('^TorsCheby%i/%s$' %(m,'.'.join(types)))
                        self.valence.set_params(term.index, fc=0.0, rv0=rv0, rv1=rv1)
                        for index in term.slaves: self.valence.set_params(index, fc=0.0, rv0=rv0, rv1=rv1)

    def do_squarebend(self, thresshold=10*deg):
        '''
            Identify bend patterns in which 4 atoms of type A surround a central
            atom of type B with A-B-A angles of 90/180 degrees. A simple
            harmonic pattern will not be adequate since a rest value of 90 and
            180 degrees is possible for the same A-B-A term. Therefore, a
            cosine term with multiplicity of 4 is used (which corresponds to a
            chebychev4 potential with sign=-1):

                  V = K/2*[1-cos(4*theta)]

            To identify the patterns, it is assumed that the rest values have
            already been estimated from the perturbation trajectories. For each
            master and slave of a BENDAHARM term, its rest value is computed and
            checked if it lies either the interval [90-thresshold,90+thresshold]
            or [180-thresshold,180]. If this is the case, the new cosine term
            is used.

            **Optional arguments**

            thresshold
                the (half) the width of the interval around 180 deg (90 degrees)
                to check if a square BA4
        '''
        for master in self.valence.iter_masters(label='BendAHarm'):
            rvs = np.zeros([len(master.slaves)+1], float)
            rvs[0] = self.valence.get_params(master.index, only='rv')
            for i, islave in enumerate(master.slaves):
                rvs[1+i] = self.valence.get_params(islave, only='rv')
            n90 = 0
            n180 = 0
            nother = 0
            for i, rv in enumerate(rvs):
                if 90*deg-thresshold<=rv and rv<=90*deg+thresshold: n90 += 1
                elif 180*deg-thresshold<=rv and rv<=180*deg+thresshold: n180 += 1
                else: nother += 1
            if n90>0 and n180>0:
                log.dump('%s has rest values around 90 deg and 180 deg, converted to BendCheby4' %master.basename)
                #modify master and slaves
                indices = [master.index]
                for slave in master.slaves: indices.append(slave)
                for index in indices:
                    term = self.valence.terms[index]
                    self.valence.modify_term(
                        index,
                        Chebychev4, [BendCos(*term.get_atoms())],
                        term.basename.replace('BendAHarm', 'BendCheby4'),
                        ['HC_FC_DIAG'], ['kjmol', 'au']
                    )
                    self.valence.set_params(index, sign=-1)
                    for traj in self.trajectories:
                        if traj.term.index==index:
                            traj.active = False
                            traj.fc = None
                            traj.rv = None

    def do_bendclin(self, thresshold=5*deg):
        '''
            No Harmonic bend can have a rest value equal to are large than
            180 deg - thresshold. If a master (or its slaves) has such a rest
            value, convert master and all slaves to BendCLin (which corresponds
            to a chebychev1 potential with sign=+1):

                0.5*K*[1+cos(theta)]
        '''
        for master in self.valence.iter_masters(label='BendAHarm'):
            indices = [master.index]
            for slave in master.slaves: indices.append(slave)
            found = False
            for index in indices:
                rv = self.valence.get_params(index, only='rv')
                if rv>=180.0*deg-thresshold:
                    found = True
                    break
            if found:
                log.dump('%s has rest value > 180-%.0f deg, converted to BendCheby1' %(master.basename, thresshold/deg))
                for index in indices:
                    term = self.valence.terms[index]
                    self.valence.modify_term(
                        index,
                        Chebychev1, [BendCos(*term.get_atoms())],
                        term.basename.replace('BendAHarm', 'BendCheby1'),
                        ['HC_FC_DIAG'], ['kjmol', 'au']
                    )
                    self.valence.set_params(index, fc=0.0, sign=1.0)
                    for traj in self.trajectories:
                        if traj.term.index==index:
                            traj.rv = None
                            traj.fc = None
                            traj.active = False

    def do_sqoopdist_to_oopdist(self, thresshold=1e-4*angstrom):
        '''
            Transform a SqOopdist term with a rest value that has been set to
            zero, to a term Oopdist (harmonic in Oopdist instead of square of
            Oopdist) with a rest value of 0.0 A.
        '''
        for master in self.valence.iter_masters(label='SqOopdist'):
            indices = [master.index]
            for slave in master.slaves: indices.append(slave)
            found = False
            for index in indices:
                rv = self.valence.get_params(index, only='rv')
                if rv<=thresshold:
                    found = True
                    break
            if found:
                log.dump('%s has rest value <= %.0f A^2, converted to Oopdist with d0=0' %(master.basename, thresshold/angstrom))
                for index in indices:
                    term = self.valence.terms[index]
                    self.valence.modify_term(
                        index,
                        Harmonic, [OopDist(*term.get_atoms())],
                        term.basename.replace('SqOopdist', 'Oopdist'),
                        ['HC_FC_DIAG'], ['kjmol/A**2', 'A']
                    )
                    self.valence.set_params(index, fc=0.0, rv0=0.0)

    def run(self):
        '''
            Sequence of instructions, should be implemented in the inheriting
            classes. The various inheriting classes distinguish themselves by
            means of the instructions implemented in this routine.
        '''
        raise NotImplementedError


class MakeTrajectories(BaseProgram):
    '''
        Construct the perturbation trajectories and store them. This program
        does not derive the force field.
    '''
    def run(self):
        with log.section('PROGRAM', 2):
            assert self.settings.fn_traj is not None, 'It is useless to run the MakeTrajectories program without specifying a trajectory filename fn_traj!'
            assert not os.path.isfile(self.settings.fn_traj), 'Given file %s to store trajectories to already exists!' %fn_traj
            self.do_pt_generate()

class PlotTrajectories(BaseProgram):
    '''
        Read the perturbation trajectories, dump to XYZ files and plot the
        energy contributions. This program does not derive the force field.
    '''
    def run(self):
        with log.section('PROGRAM', 2):
            assert self.settings.fn_traj is not None, 'The PlotTrajectories program requires a trajectory filename fn_traj!'
            assert os.path.isfile(self.settings.fn_traj), 'Given file %s to read trajectories does not exists!' %fn_traj
            self.settings.set('xyz_traj', True)
            self.settings.set('plot_traj', True)
            self.do_pt_generate()
            self.do_pt_estimate()
            self.plot_trajectories()

class DeriveFF(BaseProgram):
    '''
        Derive a force field for the given system. After the hessian fit of the
        force constants, the rest values are refined by revisiting the
        perturbation trajectories with an extra a priori term representing the
        current valence contribution. Finally, the force constants are refined
        by means of a final hessian fit.
    '''
    def run(self):
        with log.section('PROGRAM', 2):
            self.do_eq_setrv(['EQ_RV'])
            self.do_pt_generate()
            self.do_pt_estimate()
            self.do_pt_postprocess()
            self.do_cross_init()
            self.do_hc_estimatefc(['HC_FC_DIAG', 'HC_FC_CROSS_ASS', 'HC_FC_CROSS_ASA'], do_mass_weighting=self.settings.do_hess_mass_weighting)
            self.do_pt_estimate(do_valence=True)
            self.do_pt_postprocess()
            self.do_hc_estimatefc(['HC_FC_DIAG', 'HC_FC_CROSS_ASS', 'HC_FC_CROSS_ASA'], do_mass_weighting=self.settings.do_hess_mass_weighting)
            self.do_hc_estimatefc([
                'HC_FC_CROSS_ASS', 'HC_FC_CROSS_ASA', 'HC_FC_CROSS_DSS',
                'HC_FC_CROSS_DSD', 'HC_FC_CROSS_DAA', 'HC_FC_CROSS_DAD'
            ], logger_level=1, do_mass_weighting=self.settings.do_hess_mass_weighting, do_svd=self.settings.do_cross_svd)
            self.make_output()
