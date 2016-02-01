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

import os, cPickle, numpy as np, datetime

__all__ = ['BaseProgram', 'Program']

class BaseProgram(object):
    def __init__(self, system, ai, **kwargs):
        #TODO: finish documentation
        '''
            **Arguments**
            
            system
                a Yaff system object defining the system
            
            ai
                a Reference instance corresponding to the ab initio input data
            
            **Keyword Arguments**
            
            ffrefs
                a list of Reference objects corresponding to a priori determined
                contributions to the force field (such as eg. electrostatics
                or van der Waals contributions)
            
            verbose
                a flag to increase the verbosity
            
            fn_traj
                a cPickle filename to read/write the perturbation trajectories
                from/to. If the file exists, the trajectories are read from the
                file. If the file does not exist, the trajectories are written 
                to the file.
            
            plot_traj
                if set to True, all energy contributions along each perturbation
                trajectory will be plotted using the final force field.
            
            xyz_traj
                if set to True, each perturbation trajectory will be written to
                an XYZ file.
        '''
        self.system = system
        self.ai = ai
        self.kwargs = kwargs
        self.valence = ValenceFF(system)
        self.perturbation = RelaxedStrain(system)

    def reset_system(self):
        '''
            routine to reset the system coords to the ai equilbrium
        '''
        self.system.pos = self.ai.coords0.copy()
        self.valence.dlist.forward()
        self.valence.iclist.forward()
    
    def do_pt_generate(self):
        '''
            Generate perturbation trajectories.
        '''
        #read if an existing file was specified through fn_traj
        fn_traj = self.kwargs.get('fn_traj', None)
        if fn_traj is not None and os.path.isfile(fn_traj):
            trajectories = cPickle.load(open(fn_traj, 'r'))
            if self.kwargs.get('verbose', False): print 'perturbation trajectories read from file %s' %fn_traj
            return trajectories
        #configure
        self.reset_system()
        do_terms = [term for term in self.valence.terms if 'PT_ALL' in term.tasks]
        trajectories = self.perturbation.prepare(self.valence, do_terms)
        #compute
        trajectories = paracontext.map(self.perturbation.generate, trajectories)
        #write the trajectories to the non-existing file fn_traj
        if fn_traj is not None:
            assert not os.path.isfile(fn_traj)
            cPickle.dump(trajectories, open(fn_traj, 'w'))
        return trajectories
        
    def do_pt_estimate(self, trajectories, do_valence=False):
        '''
            Estimate force constants and rest values from the perturbation
            trajectories
            
            **Optional Arguments**
            
            do_valence
                if set to True, the current valence force field will be used to
                estimate the contribution of all other valence terms.
        '''
        self.reset_system()
        ffrefs = self.kwargs.get('ffrefs', [])
        for traj in trajectories:
            #compute fc and rv from trajectory
            if do_valence:
                self.perturbation.estimate(traj, self.ai, ffrefs=ffrefs, valence=self.valence)
            else:
                self.perturbation.estimate(traj, self.ai, ffrefs=ffrefs)
        #set force field parameters to computed fc and rv
        for traj in trajectories:
            self.valence.set_params(traj.term.index, fc=traj.fc, rv0=traj.rv)
    
    def do_eq_setrv(self, tasks):
        '''
            Set the rest values to their respective AI equilibrium values.
        '''
        self.reset_system()
        for term in self.valence.terms:
            vterm = self.valence.vlist.vtab[term.index]
            flagged = False
            for flag in tasks:
                if flag in term.tasks:
                    flagged = True
                    break
            if flagged:
                if term.kind==3:#cross term
                    ic0 = self.valence.iclist.ictab[vterm['ic0']]
                    ic1 = self.valence.iclist.ictab[vterm['ic1']]
                    self.valence.set_params(term.index, rv0=ic0['value'], rv1=ic1['value'])
                elif term.kind==4 and term.ics[0].kind==4:#Cosine of DihedAngle
                    ic = self.valence.iclist.ictab[vterm['ic0']]
                    m = self.valence.get_params(term.index, only='m')
                    rv = ic['value']%(360.0*deg/m)
                    self.valence.set_params(term.index, rv0=rv)
                else:
                    rv = self.valence.iclist.ictab[vterm['ic0']]['value']
                    self.valence.set_params(term.index, rv0=rv)
    
    def do_hc_estimatefc(self, tasks):
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
        '''
        self.reset_system()
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

    def do_cross_init(self):
        '''
            Set the rest values of cross terms to the rest values of the
            corresponding diagonal terms. The force constants are initialized
            to zero.
        '''
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
    
    def do_average_pars(self):
        '''
            Average force field parameters over master and slaves.
        '''
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
    
    def run(self):
        '''
            Sequence of instructions, should be implemented in the inheriting
            classes. The various inheriting classes distinguish themselves by 
            means of the instructions implemented in this routine.
        '''
        raise NotImplementedError


class Program(BaseProgram):
    def run(self):
        time0 = datetime.datetime.now()
        verbose = self.kwargs.get('verbose', False)
        
        if verbose: print 'generating perturbation trajectories ...'
        trajectories = self.do_pt_generate()
        time1 = datetime.datetime.now()
        
        if verbose: print 'estimating pars from perturbation trajectories ...'
        self.do_pt_estimate(trajectories)
        self.valence.dump_screen()
        self.do_average_pars()
        self.valence.dump_yaff('pars_pert.txt')
        time2 = datetime.datetime.now()
        
        if verbose: print 'estimating diagonal force constants from Hessian ...'
        self.do_hc_estimatefc(['HC_FC_DIAG'])
        self.valence.dump_screen()
        self.valence.dump_yaff('pars_hess_diag.txt')
        time3 = datetime.datetime.now()
        
        if verbose: print 'revisiting perturbation trajectories ...'
        self.do_pt_estimate(trajectories, do_valence=True)
        self.valence.dump_screen()
        self.do_average_pars()
        self.valence.dump_yaff('pars_pert2.txt')
        
        if verbose: print 'reestimate diagonal force constants from Hessian ...'
        self.do_hc_estimatefc(['HC_FC_DIAG'])
        self.valence.dump_screen()
        self.valence.dump_yaff('pars_hess_diag2.txt')
        time4 = datetime.datetime.now()
        
        if verbose: print 'estimating cross force constants from Hessian ...'
        self.do_cross_init()
        self.do_hc_estimatefc(['HC_FC_CROSS'])
        self.valence.dump_screen()
        self.valence.dump_yaff('pars_hess_cross.txt')
        time5 = datetime.datetime.now()
        
        if verbose: print 'estimating all force constants from Hessian ...'
        self.do_hc_estimatefc(['HC_FC_DIAG','HC_FC_CROSS'])
        self.valence.dump_screen()
        self.valence.dump_yaff('pars_hess_all.txt')
        time6 = datetime.datetime.now()
        
        if verbose: print 'dumping additional output ...'
        if self.kwargs.get('plot_traj', False):
            for trajectory in trajectories:
                trajectory.plot(self.ai, ffrefs=self.kwargs.get('ffrefs', []), valence=self.valence)
        if self.kwargs.get('xyz_traj', False):
            for trajectory in trajectories:
                trajectory.to_xyz()
        self.system.to_file('system.chk')
        time7 = datetime.datetime.now()
        
        if verbose:
            timings = {
                'Generate PT':          time1-time0,
                'Estimate PT':          time2-time1,
                'Estimate Cost diag':   time3-time2,
                'Revisit  PT':          time4-time3,
                'Estimate Cost cross':  time5-time4,
                'Estimate Cost all':    time6-time5,
                'Dumping  Output':       time7-time6,
            }
            labels = ['Generate PT', 'Estimate PT', 'Estimate Cost diag', 'Revisit  PT', 'Estimate Cost cross', 'Estimate Cost all', 'Dumping  Output']
            print ''
            print 'TIMING:'
            print '-'*45
            for label in labels:
                line = '%30s  ' %(label+' '*(30-len(label)))
                seconds = timings[label].total_seconds()
                line += '  %10.1f' %seconds
                print line
            print '-'*45
