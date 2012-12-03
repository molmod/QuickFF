#!/usr/bin/env python

import os, numpy as np
from cfit2 import *
from mfit2 import *

__all__ = ['FFitProgram', 'generate_terms']

def generate_terms(system):
    terms = []    
    for icname, ics in system.ics.iteritems():
        kind = icname.split('/')[0]
        types = icname.split('/')[1].split('.')
        if kind in ['bond', 'dist']:
            for ic in ics:
                terms.append(HarmonicTerm(BondFilter(ic.indexes[0], ic.indexes[1]), dist))
        elif kind in ['bend', 'angle']:
            for ic in ics:
                terms.append(HarmonicTerm(AngleFilter(ic.indexes[0], ic.indexes[1], ic.indexes[2]), angle))
        elif kind in ['dihedral', 'dihed', 'torsion']:
            for ic in ics:
                terms.append(CosineTerm(DihedralFilter(ic.indexes[0], ic.indexes[1], ic.indexes[2], ic.indexes[3]), dihed, 2, 0.0))
        else:
            return ValueError('Invalied IC type, recieved %s' %kind)
    #terms.append( GaussianChargeTerm(PairFilter('*', '*', exclude_bonded=0), dist)  )
    return terms


class FFitProgram(DefaultProgram):
    def __init__(self, system):
        rules = []
        dn_out = os.path.getcwd()
        fn_pars_init = 'pars_init.txt'
        lsqs = [HessianLSQ('Hess', system.fn_chk, 1.0)]
        terms = generate_terms(system)
        DefaultProgram.__init__(self, dn_out, terms, rules, lsqs, fn_pars_init, fix_patterns=[], free_patterns=[],
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
            if isinstance(term, CosineTerm):
                if term.pars[0].free:
                    add_pos(term.pars[0])
        return constraints
    
    def run(self):
        model = DefaultProgram.run(self)
        model.dump_pars('pars_fine.txt')
        return model
