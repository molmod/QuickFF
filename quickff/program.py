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

import os, cPickle, numpy as np

__all__ = ['Program']

class BaseProgram(object):
    def __init__(self, system, refs, **kwargs):
        #TODO: finish documentation
        '''
            **Arguments**
            
            system
                a Yaff system object defining the system
            
            refs
                a list of Reference objects of which one should be identified
                as the ab initio reference by including the string 'ai' in its
                name.
            
            
            **Keyword Arguments**
            
            term_specs
            
            do_taylor
            
            fn_traj
        '''
        self.system = system
        self.refs = refs
        self.kwargs = kwargs
        self.valence = ValenceFF(system, specs=kwargs.get('term_specs', None))
        self.perturbation = RelaxedStrain(
            system, refs, do_taylor=kwargs.get('do_taylor', None)
        )

    def reset_system(self):
        '''
            routine to reset the system coords to the ai equilbrium
        '''
        for ref in self.refs:
            if 'ai' in ref.name.lower():
                self.system.pos = ref.coords0.copy()
                self.valence.dlist.forward()
                self.valence.iclist.forward()
                return
        raise IOError("No Ab Initio reference found. Be sure to add the string 'ai' to its name.")
    
    def do_pt_generate(self):
        '''
            Generate perturbation trajectories.
        '''
        #configure
        self.reset_system()
        do_terms = [term for term in self.valence.terms if 'PT_ALL' in term.tasks]
        trajectories = self.perturbation.prepare(self.valence, do_terms)
        #compute
        paracontext.map(self.perturbation.generate, trajectories)
        return trajectories
        
    def do_pt_estimate(self, trajectories):
        '''
            Estimate force constants and rest values from the perturbation
            trajectories
        '''
        self.reset_system()
        for traj in trajectories:
            #compute fc and rv from trajectory
            self.perturbation.estimate(traj)
            #set force field parameters to computed fc and rv
            self.valence.set_params(traj.term.index, fc=traj.fc, rv0=traj.rv)
    
    def do_eq_setrv(self, task_flag):
        '''
            Set the rest values to their respective AI equilibrium values.
        '''
        self.reset_system()
        for index in xrange(self.valence.vlist.nv):
            term = self.valence.vlist.vtab[index]
            if task_flag is None or task_flag in self.valence.terms[index].tasks:
                if term['kind']==3:#cross term
                    ic0 = self.valence.iclist.ictab[term['ic0']]
                    ic1 = self.valence.iclist.ictab[term['ic1']]
                    self.valence.set_params(index, rv0=ic0['value'], rv1=ic1['value'])
                else:
                    ic = self.valence.iclist.ictab[term['ic0']]
                    self.valence.set_params(index, rv0=ic['value'])
    
    def do_hc_estimatefc(self, task_flag):
        '''
            Refine force constants using Hessian Cost function.
            
            **Arguments**
            
            task_flag
                A string identifying which terms should have their force
                constant estimated from the hessian cost function. Using such a
                flag, one can distinguish between for example force constant
                refinement (flag=HC_FC_DIAG) of the diagonal terms and force 
                constant estimation of the cross terms (flag=HC_FC_CROSS).
        '''
        self.reset_system()
        term_indices = []
        for index in xrange(self.valence.vlist.nv):
            term = self.valence.terms[index]
            if task_flag in term.tasks:
                #first check if all rest values and multiplicities have been defined
                if term.kind==0: self.valence.check_params(term, ['rv'])
                if term.kind==3: self.valence.check_params(term, ['rv0','rv1'])
                if term.kind==4: self.valence.check_params(term, ['rv', 'm'])
                if term.is_master():
                    term_indices.append(index)
            else:
                #first check if all pars have been defined
                if term.kind==0: self.valence.check_params(term, ['fc', 'rv'])
                if term.kind==3: self.valence.check_params(term, ['fc', 'rv0','rv1'])
                if term.kind==4: self.valence.check_params(term, ['fc', 'rv', 'm'])
        cost = HessianFCCost(self.system, self.refs, self.valence, term_indices)   
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
            Average force field parameters over master and slaves. Small
            negative rest values of out-of-plane distance terms will
            also be set to zero.
        '''
        for master in self.valence.iter_masters():
            npars = len(self.valence.get_params(master.index))
            pars = np.zeros([len(master.slaves)+1, npars], float)
            pars[0,:] = np.array(self.valence.get_params(master.index))
            for i, islave in enumerate(master.slaves):
                pars[1+i,:] = np.array(self.valence.get_params(islave))
            if master.kind==0:#harmonic
                fc = pars[:,0].mean()
                rv = pars[:,1].mean()
                if master.ics[0].kind==10:#opdist
                    if rv<0.0 and abs(rv)<1e-2*angstrom:
                        rv=0.0
                self.valence.set_params(master.index, fc=fc, rv0=rv)
                for islave in master.slaves:
                    self.valence.set_params(islave, fc=fc, rv0=rv)
            elif master.kind==3:#cross
                fc = pars[:,0].mean()
                rv0 = pars[:,1].mean()
                rv1 = pars[:,2].mean()
                self.valence.set_params(master.index, fc=fc, rv0=rv0, rv1=rv1)
                for islave in master.slaves:
                    self.valence.set_params(islave, fc=fc, rv0=rv0, rv1=rv1)
            elif master.kind==4:#cosine
                m = pars[:,0].mean()
                assert pars[:,0].std()<1e-6, 'dihedral multiplicity not unique'
                fc = pars[:,1].mean()
                rv = pars[:,2].mean()
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
        trajectories = self.do_pt_generate()
        self.do_pt_estimate(trajectories)
        self.do_eq_setrv('EQ_RV')
        self.do_average_pars()
        self.do_hc_estimatefc('HC_FC_DIAG')
        self.do_cross_init()
        self.do_hc_estimatefc('HC_FC_CROSS')
