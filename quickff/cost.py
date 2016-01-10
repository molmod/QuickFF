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
from yaff.pes.ff import ForceField, ForcePartValence
from yaff.pes.vlist import Harmonic, Fues, Cross, Cosine
from yaff.sampling.harmonic import estimate_cart_hessian
from molmod.units import *
from quickff.tools import dihed_to_chebychev, boxqp
import numpy as np

__all__ = ['HessianFCCost']


kind_to_term = {0: Harmonic, 2: Fues, 3: Cross, 4: Cosine}

class HessianFCCost(object):
    def __init__(self, system, refs, valence, term_indices):
        '''
            **Arguments**
            
            system
                a Yaff system object
            
            refs
                a list of ReferenceData objecst, one of which should be the
                ab initio reference indicated with the string 'ai' in its
                title attribute.
            
            valence
                A ValenceFF object containing all valence terms.
            
            term_indices
                a list of indices indicating the terms for which the force
                constants should be determined.
        
        '''
        #initialization
        self.init = np.zeros(len(term_indices), float)
        self.upper = np.zeros(len(term_indices), float)+np.inf
        self.lower = np.zeros(len(term_indices), float)
        self.A = np.zeros([len(term_indices), len(term_indices)], float)
        self.B = np.zeros([len(term_indices)], float)
        
        #compute the reference hessian (including valence terms not fitted)
        ai = None
        ffrefs = []
        for ref in refs:
            if 'ai' in ref.name.lower():
                assert ai is None
                ai = ref
            else:
                ffrefs.append(ref)
        if ai is None:
            raise IOError("No Ab Initio reference found. Be sure to add the string 'ai' to its name.")
        ndofs = 3*system.natom
        href = ai.phess0.reshape([ndofs, ndofs]).copy()
        for ffref in ffrefs:
            href -= ffref.phess0.reshape([ndofs, ndofs])
        for index in xrange(valence.vlist.nv):
            #skip all masters (and its slaves) given in term_indices
            skip = False
            for imaster in term_indices:
                master = valence.data[imaster]
                if index==imaster or index in master.slaves:
                    skip = True
            if skip: continue
            #compute covalent contribution to reference
            num_hessian = False
            val = ForcePartValence(system)
            pars = valence.get_params(index)
            ics = valence.data[index].ics
            kind = valence.vlist.vtab[index]['kind']
            if kind==4:#Cosine
                chebychev = dihed_to_chebychev(pars, ics[0])
                if chebychev is not None:
                    val.add_term(chebychev)
                else:
                    num_hessian = True
                    args = tuple(pars) + tuple(ics)
                    val.add_term(Cosine(*args))
            else:
                args = tuple(pars) + tuple(ics)
                val.add_term(kind_to_term[kind](*args))
            ff = ForceField(system, [val])
            if num_hessian:
                hcov = estimate_cart_hessian(ff)
            else:
                hcov = np.zeros([ndofs, ndofs], float)
                ff.compute(hess=hcov)
            href -= hcov
        
        #construct the cost matrices A and B
        hcovs = []
        for i, iterm in enumerate(term_indices):
            if not valence.data[iterm].is_master():
                raise ValueError('term_indices should only contain masters')
            num_hessian = False
            kind = valence.vlist.vtab[iterm]['kind']
            fc = valence.get_params(iterm, only='fc')
            if not np.isnan(fc):
                self.init[i] = fc
            val = ForcePartValence(system)
            if kind==4:#Cosine
                self.upper[i] = 200.0*kjmol
                m, fc, rv = valence.get_params(iterm)
                for islave in valence.data[iterm].slaves:
                    ics = valence.data[islave].ics
                    chebychev = dihed_to_chebychev([m, 1.0, rv], ics[0])
                    if chebychev is not None:
                        val.add_term(chebychev)
                    else:
                        num_hessian = True
                        args = (m, 1.0, rv) + tuple(ics)
                        val.add_term(Cosine(*args))
            elif kind==3:#cross
                fc, rv0, rv1 = valence.get_params(iterm)
                for islave in valence.data[iterm].slaves:
                    ics = valence.data[islave].ics
                    args = (1.0,rv0, rv1) + tuple(ics)
                    val.add_term(Cross(*args))
            else:
                fc, rv = valence.get_params(iterm)
                for islave in valence.data[iterm].slaves:
                    ics = valence.data[islave].ics
                    args = (1.0,rv) + tuple(ics)
                    val.add_term(kind_to_term[kind](*args))
            ff = ForceField(system, [val])
            if num_hessian:
                hcov = estimate_cart_hessian(ff)
            else:
                hcov = np.zeros([ndofs, ndofs], float)
                ff.compute(hess=hcov)
            hcovs.append(hcov)
        
        for index1, hcov1 in enumerate(hcovs):
            self.B[index1] = np.sum(href*hcov1)
            self.A[index1,index1] = np.sum(hcov1*hcov1)
            for index2, hcov2 in enumerate(hcovs[:index1]):
                tmp = np.sum(hcov1*hcov2)
                self.A[index1,index2] = tmp
                self.A[index2,index1] = tmp

    def estimate(self, init=None, lower=None, upper=None):
        '''
            Estimate the force constants by minimizing the cost function

        '''
        if init is None:
            assert self.init is not None, 'No initial fcs defined'
            init = self.init.copy()
        if lower is None:
            assert self.lower is not None, 'No lower limit fcs defined'
            lower = self.lower.copy()
        if upper is None:
            assert self.upper is not None, 'No upper limit fcs defined'
            upper = self.upper.copy()
        x = boxqp(self.A, self.B, lower, upper, init)
        return x
