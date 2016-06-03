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
from quickff.log import log

from yaff.pes.vlist import Cosine, Harmonic
from yaff.pes.iclist import BendAngle, BendCos, OopDist

import os, cPickle, numpy as np, datetime

__all__ = [
    'BaseProgram', 'MakeTrajectories', 'PlotTrajectories',
    'DeriveDiagFF', 'DeriveNonDiagFF'
]

class BaseProgram(object):
    '''
        Base program which implements all possible steps of a force field
        fitting program. The actual sequence of the steps are defined in the
        deriving classes.
    '''
    def __init__(self, system, ai, **kwargs):
        '''
            **Arguments**

            system
                a Yaff `System` object defining the system

            ai
                a `Reference` instance corresponding to the ab initio input data

            **Keyword Arguments**

            ffrefs
                a list of `Reference` objects corresponding to a priori determined
                contributions to the force field (such as eg. electrostatics
                or van der Waals contributions)

            fn_yaff
                the name of the file to write the final parameters to in Yaff
                format. The default is `pars.txt`.

            fn_sys
                the name of the file to write the system to. The default is
                `system.chk`.

            fn_traj
                a cPickle filename to read/write the perturbation trajectories
                from/to. If the file exists, the trajectories are read from the
                file. If the file does not exist, the trajectories are written
                to the file.

            only_traj
                specifier to determine for which terms a perturbation trajectory
                needs to be constructed. If ONLY_TRAJ is a single string, it is
                interpreted as a task (only terms that have this task in their
                tasks attribute will get a trajectory). If ONLY_TRAJ is a list
                of strings, each string is interpreted as the basename of the
                term for which a trajectory will be constructed.

            plot_traj
                if set to True, all energy contributions along each perturbation
                trajectory will be plotted using the final force field.

            xyz_traj
                if set to True, each perturbation trajectory will be written to
                an XYZ file.
        '''
        with log.section('PROG', 2, timer='Initializing'):
            log.dump('Initializing program')
            self.system = system
            self.ai = ai
            self.kwargs = kwargs
            self.valence = ValenceFF(system)
            self.perturbation = RelaxedStrain(system, self.valence)
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
        with log.section('PTUPD', 4):
            for traj in self.trajectories:
                found = False
                for term in self.valence.iter_terms():
                    if traj.term.get_atoms()==term.get_atoms():
                        if found: raise ValueError('Found two trajectories for term %i (%s) with atom indices %s' %(term.index, term.basename, str(term.get_atoms())))
                        traj.term = term
                        if 'PT_ALL' not in term.tasks: 
                            log.dump('PT_ALL not in tasks of %s-%i, deactivated PT' %(term.basename, term.index))
                            traj.active = False
                        found = True
                        break
                if not found: log.dump('WARNING: No trajectory found for term %i (%s) with atom indices %s' %(term.index, term.basename, str(term.get_atoms())))

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
            if master.kind==0:#harmonic
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
            else:
                raise NotImplementedError

    def make_output(self):
        '''
            Dump Yaff parameters, Yaff system, plot energy contributions along
            perturbation trajectories and dump perturbation trajectories to XYZ
            files.
        '''
        fn_yaff = self.kwargs.get('fn_yaff', None)
        if fn_yaff is None:
            fn_yaff = 'pars_cov%s.txt' %(self.kwargs.get('suffix', ''))
        self.valence.dump_yaff(fn_yaff)
        fn_sys = self.kwargs.get('fn_sys', None)
        if fn_sys is None:
            fn_sys = 'system%s.chk' %(self.kwargs.get('suffix', ''))
        self.system.to_file(fn_sys)
        self.plot_trajectories(do_valence=True)

    def plot_trajectories(self, do_valence=False):
        '''
            Plot energy contributions along perturbation trajectories and dump
            perturbation trajectories to XYZ files.
        '''
        with log.section('PLOT', 4, timer='PT plot energy'):
            if self.kwargs.get('plot_traj', False):
                ffrefs = self.kwargs.get('ffrefs', [])
                valence = None
                if do_valence: valence=self.valence
                for trajectory in self.trajectories:
                    if trajectory is not None:
                        trajectory.plot(self.ai, ffrefs=ffrefs, valence=valence)
        with log.section('XYZ', 4, timer='PT dump XYZ'):
            if self.kwargs.get('xyz_traj', False):
                for trajectory in self.trajectories:
                    if trajectory is not None:
                        trajectory.to_xyz()

    def do_pt_generate(self, do='PT_ALL'):
        '''
            Generate perturbation trajectories.

            **Optional Arguments**

            do
              List of term basenames for which the perturbation trajectories
              should be constructed. Can also be a string specifying a task,
              a trajectory will then be constructed for each term that has this
              task in its tasks attribute.
        '''
        with log.section('PTGEN', 2, timer='PT Generate'):
            #read if an existing file was specified through fn_traj
            fn_traj = self.kwargs.get('fn_traj', None)
            if fn_traj is not None and os.path.isfile(fn_traj):
                self.trajectories = cPickle.load(open(fn_traj, 'r'))
                log.dump('Trajectories read from file %s' %fn_traj)
                self.update_trajectory_terms()
                return
            #configure
            self.reset_system()
            if isinstance(do, str):
                do_terms = [term for term in self.valence.terms if do in term.tasks]
            elif isinstance(do, list):
                do_terms = []
                for master in do:
                    for term in self.valence.iter_terms(master):
                        do_terms.append(term)
            else:
                raise IOError("Invalid value for optional argument 'do', recieved %s" %str(do))
            trajectories = self.perturbation.prepare(do_terms)
            #compute
            log.dump('Constructing trajectories')
            self.trajectories = paracontext.map(self.perturbation.generate, [traj for traj in trajectories if traj is not None])
            #write the trajectories to the non-existing file fn_traj
            if fn_traj is not None:
                assert not os.path.isfile(fn_traj)
                cPickle.dump(self.trajectories, open(fn_traj, 'w'))
                log.dump('Trajectories stored to file %s' %fn_traj)

    def do_pt_estimate(self, do_valence=False, plot=False):
        '''
            Estimate force constants and rest values from the perturbation
            trajectories

            **Optional Arguments**

            do_valence
                if set to True, the current valence force field will be used to
                estimate the contribution of all other valence terms.
            
            plot
                if set to True, plots the energy contributions immediately
                after estimating them (before averaging).
        '''
        with log.section('PTEST', 2, timer='PT Estimate'):
            self.reset_system()
            message = 'Estimating FF parameters from perturbation trajectories'
            if do_valence: message += ' with valence reference'
            log.dump(message)
            ffrefs = self.kwargs.get('ffrefs', [])
            #compute fc and rv from trajectory
            for traj in self.trajectories:
                if traj is None: continue
                self.perturbation.estimate(traj, self.ai, ffrefs=ffrefs, do_valence=do_valence)
            #plot trajectories if needed
            if plot: self.plot_trajectories(do_valence=do_valence)
            #set force field parameters to computed fc and rv
            for traj in self.trajectories:
                if traj is None: continue
                self.valence.set_params(traj.term.index, fc=traj.fc, rv0=traj.rv)
            #output
            self.valence.dump_logger(print_level=4)
            self.do_squarebend()
            self.do_bendcharm()
            #self.do_sqoopdist_to_oopdist()
            self.average_pars()

    def do_eq_setrv(self, tasks):
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
                        with log.section('EQSET', 3, timer='Equil Set RV'):
                            log.dump('Set rest value of %s(%s) (eq=%.3f deg) to %.3f deg' %(
                                term.basename,
                                '.'.join([str(at) for at in term.get_atoms()]),
                                ic['value']/deg, rv/deg
                            ))
                        self.valence.set_params(term.index, rv0=rv)
                    else:
                        rv = self.valence.iclist.ictab[vterm['ic0']]['value']
                        self.valence.set_params(term.index, rv0=rv)
            self.valence.dump_logger(print_level=4)
            self.average_pars()

    def do_hc_estimatefc(self, tasks, logger_level=4):
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
        '''
        with log.section('HCEST', 2, timer='HC Estimate FC'):
            self.reset_system()
            log.dump('Estimating force constants from Hessian cost for tasks %s' %' '.join(tasks))
            ffrefs = self.kwargs.get('ffrefs', [])
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
            cost = HessianFCCost(self.system, self.ai, self.valence, term_indices, ffrefs=ffrefs)
            fcs = cost.estimate()
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
        with log.section('CRINI', 2):
            self.reset_system()
            self.valence.init_cross_terms()
            for index in xrange(self.valence.vlist.nv):
                term = self.valence.vlist.vtab[index]
                if term['kind']!=3: continue
                rv0, rv1 = None, None
                for index2 in xrange(self.valence.vlist.nv):
                    term2 = self.valence.vlist.vtab[index2]
                    if term2['kind']==3: continue
                    if term['ic0']==term2['ic0']:
                        assert rv0 is None
                        rv0 = self.valence.get_params(index2, only='rv')
                    if term['ic1']==term2['ic0']:
                        assert rv1 is None
                        rv1 = self.valence.get_params(index2, only='rv')
                if rv0 is None or rv1 is None:
                    raise ValueError('No rest values found for %s' %self.valence.terms[index].basename)
                self.valence.set_params(index, fc=0.0, rv0=rv0, rv1=rv1)
    
    def do_squarebend(self, thresshold=10*deg):
        '''
            Identify bend patterns in which 4 atoms of type A surround a central
            atom of type B with A-B-A angles of 90/180 degrees. A simple
            harmonic pattern will not be adequate since a rest value of 90 and
            180 degrees is possible for the same A-B-A term. Therefore, a
            cosine term with multiplicity of 4 is used:
           
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
        with log.section('SQBEND', 2):
            for master in self.valence.iter_masters(label='BendAHarm'):
                rvs = np.zeros([len(master.slaves)+1], float)
                rvs[0] = self.valence.get_params(master.index, only='rv')
                for i, islave in enumerate(master.slaves):
                    rvs[1+i] = self.valence.get_params(islave, only='rv')
                n90 = 0
                n180 = 0
                nother = 0
                for rv in rvs:
                    if 90*deg-thresshold<=rv and rv<=90*deg+thresshold: n90 += 1
                    elif 180*deg-thresshold<=rv and rv<=180*deg+thresshold: n180 += 1
                    else: nother += 1
                if nother==0 and n90>0 and n180>0:
                    log.dump('%s has rest values around 90 deg and 180 deg, converted to BendCos with m=4' %master.basename)
                    #modify master and slaves
                    indices = [master.index]
                    for slave in master.slaves: indices.append(slave)
                    for index in indices:
                        term = self.valence.terms[index]
                        self.valence.modify_term(
                            index,
                            Cosine, [BendAngle(*term.get_atoms())],
                            term.basename.replace('BendAHarm', 'BendCos'),
                            ['HC_FC_DIAG'], ['au', 'kjmol', 'deg']
                        )
                        self.valence.set_params(index, rv0=0.0, m=4)
                        for traj in self.trajectories:
                            if traj.term.index==index:
                                traj.active = False
                                traj.fc = None
                                traj.rv = None
    
    def do_bendcharm(self, thresshold=2*deg):
        '''
            No Harmonic bend can have a rest value equal to are large than 
            180 deg - thresshold. If a master (or its slaves) has such a rest 
            value, convert master and all slaves to BendCharm with
            cos(phis0)=-1.
        '''
        with log.section('BNDCHRM', 2):
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
                    log.dump('%s has rest value > 180-%.0f deg, converted to BendCHarm with cos(phi0)=-1' %(master.basename, thresshold/deg))
                    for index in indices:
                        term = self.valence.terms[index]
                        self.valence.modify_term(
                            index,
                            Harmonic, [BendCos(*term.get_atoms())],
                            term.basename.replace('BendAHarm', 'BendCHarm'),
                            ['HC_FC_DIAG'], ['kjmol', 'au']
                        )
                        self.valence.set_params(index, fc=0.0, rv0=-1.0)
                        for traj in self.trajectories:
                            if traj.term.index==index:
                                traj.rv = None
                                traj.fc = None
                                traj.active = False

    def do_sqoopdist_to_oopdist(self,thresshold=1e-4*angstrom):
        '''
            Transform a SqOopdist term with a rest value that has been set to
            zero, to a term Oopdist (harmonic in Oopdist instead of square of
            Oopdist) with a rest value of 0.0 A.
        '''
        with log.section('SQOOP', 2):
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
            fn_traj = self.kwargs.get('fn_traj', None)
            assert fn_traj is not None, 'It is useless to run the MakeTrajectories program without specifying a trajectory filename fn_traj!'
            assert not os.path.isfile(fn_traj), 'Given file %s to store trajectories to already exists!' %fn_traj
            self.do_pt_generate(do=self.kwargs.get('only_traj', 'PT_ALL'))

class PlotTrajectories(BaseProgram):
    '''
        Read the perturbation trajectories, dump to XYZ files and plot the 
        energy contributions.
    '''
    def run(self):
        with log.section('PROGRAM', 2):
            fn_traj = self.kwargs.get('fn_traj', None)
            assert fn_traj is not None, 'The PlotTrajectories program requires a trajectory filename fn_traj!'
            assert os.path.isfile(fn_traj), 'Given file %s to read trajectories does not exists!' %fn_traj
            self.do_pt_generate()
            self.do_pt_estimate()
            self.kwargs['xyz_traj'] = True
            self.kwargs['plot_traj'] = True
            self.plot_trajectories()

class DeriveDiagFF(BaseProgram):
    '''
        Derive a diagonal force field, i.e. does not contain cross terms.
    '''
    def run(self):
        with log.section('PROGRAM', 2):
            self.do_eq_setrv(['EQ_RV'])
            self.do_pt_generate()
            self.do_pt_estimate()
            self.do_hc_estimatefc(['HC_FC_DIAG'])
            self.do_pt_estimate(do_valence=True)
            self.do_hc_estimatefc(['HC_FC_DIAG'], logger_level=1)
            self.make_output()

class DeriveNonDiagFF(BaseProgram):
    '''
        Derive a non-diagonal force field, i.e. contains cross terms.
    '''
    def run(self):
        with log.section('PROGRAM', 2):
            self.do_eq_setrv(['EQ_RV'])
            self.do_pt_generate()
            self.do_pt_estimate()
            self.do_cross_init()
            self.do_hc_estimatefc(['HC_FC_DIAG','HC_FC_CROSS'])
            self.do_pt_estimate(do_valence=True)
            self.do_hc_estimatefc(['HC_FC_DIAG','HC_FC_CROSS'], logger_level=1)
            self.make_output()
