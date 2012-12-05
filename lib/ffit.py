#!/usr/bin/env python

import os, numpy as np
from cfit2 import *
from mfit2 import *

from fftable import FFTable

__all__ = ['FFitProgram', 'generate_terms']

def generate_terms(system):
    print ' FFIT  TERMS: auto generate harmonic force field terms'
    terms = []    
    for icname, ics in system.ics.iteritems():
        kind = icname.split('/')[0]
        types = icname.split('/')[1].split('.')
        if kind in ['bond', 'dist']:
            terms.append(HarmonicTerm(BondFilter(types[0], types[1]), dist))
        elif kind in ['bend', 'angle']:
            terms.append(HarmonicTerm(AngleFilter(types[0], types[1], types[2]), angle))
        elif kind in ['dihedral', 'dihed', 'torsion']:
            for ic in ics:
                n1 = len(system.neighbor_list[ic.indexes[1]])
                n2 = len(system.neighbor_list[ic.indexes[2]])
                if   4 in [n1, n2]: m = 3
                elif 3 in [n1, n2]: m = 2
                elif 2 in [n1, n2]: m = 1
                else: raise ValueError('Dihedral %i,%i,%i,%i has no atoms bonded to central pair' %(ic.indexes[0], ic.indexes[1], ic.indexes[2], ic.indexes[3]))
            terms.append(CosineTerm(DihedralFilter(types[0], types[1], types[2], types[3]), dihed, m, 0.0))
        else:
            return ValueError('Invalied IC type, recieved %s' %kind)
    if system.eirule>-1:
        terms.append( GaussianChargeTerm(PairFilter('*', '*', exclude_bonded=system.eirule), dist)  )
    return terms


class FFitProgram(DefaultProgram):
    def __init__(self, system):
        rules = []
        dn_out = '%s/out' %(os.getcwd())
        log._file = open('out-ffit2.log', 'w')
        fn_pars_init = 'int-pars.txt'
        lsqs = [HessianLSQ('Hess', 'int-system.chk', 1.0)]
        fix_patterns = ['q0']
        terms = generate_terms(system)
        DefaultProgram.__init__(self, dn_out, terms, rules, lsqs, fn_pars_init, fix_patterns=fix_patterns, free_patterns=[],
                                step_dl_min=1e-5, step_dl_max=1e2, max_iter=50, weight_mode=None, weight_ridge=0.1)

    def get_constraints(self, model):
        constraints = []
        def add_pos(par):
            constraints.append(LowerBound(par.longname, 0.0))
        for term in model.terms:
            if isinstance(term, HarmonicTerm):
                if term.pars[0].free:
                    add_pos(term.pars[0])
                if term.pars[1].free:
                    add_pos(term.pars[1]) 
            #if isinstance(term, CosineTerm):
            #    if term.pars[0].free:
            #        add_pos(term.pars[0])
        return constraints
    
    def run(self):
        print ' FFIT  RUN  : refine force field parameters using FFit2'
        model = DefaultProgram.run(self)
        dump_yaff_parameters(model.terms, model.training_set, 'pars.txt')
        fftab = FFTable.from_ffit2(model)
        return fftab
