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
            terms.append(CosineTerm(DihedralFilter(types[0], types[1], types[2], types[3]), dihed, 2, 0.0))
        else:
            return ValueError('Invalied IC type, recieved %s' %kind)
    #terms.append( GaussianChargeTerm(PairFilter('*', '*', exclude_bonded=0), dist)  )
    return terms


class FFitProgram(DefaultProgram):
    def __init__(self, system):
        rules = []
        dn_out = '%s/tmp' %(os.getcwd())
        log._file = open('ffit2.log', 'w')
        fn_pars_init = 'pars_init.txt'
        lsqs = [HessianLSQ('Hess', 'system.chk', 1.0)]
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
            if isinstance(term, CosineTerm):
                if term.pars[0].free:
                    add_pos(term.pars[0])
        return constraints
    
    def run(self):
        print ' FFIT  RUN  : refine force field parameters using FFit2'
        model = DefaultProgram.run(self)
        model.dump_pars('pars_fine.txt')
        os.system('rm -r %s' %self.dn_out)
        print ' FFIT  RUN  : constructing FFTable'
        icnames = []
        units = {}
        kdata = {}
        qdata = {}
        for rule in model.rules:
            for par in rule.pars:
                icname = '/'.join(par.prefix.split('/')[:2])
                kind = par.name[0].lower()
                if not icname in icnames: icnames.append(icname)
                unit = units.get(icname, {'k': None, 'q': None})
                unit[kind] = par.unit.replace('^', '**')
                units[icname] = unit
                if kind=='k':
                    kdata[icname] = np.array([par.value])
                elif kind=='q':
                    qdata[icname] = np.array([par.value])
                else:
                    raise ValueError('Invalid value for kind, recieved %s' %kind)
        fftab = FFTable(icnames, units)
        for icname in icnames:
            if icname.startswith('dihed/'):
                qdata[icname] = np.array([0.0])
                units[icname]['q'] = 'deg'
            fftab.add(icname, kdata[icname], qdata[icname])
        return fftab
