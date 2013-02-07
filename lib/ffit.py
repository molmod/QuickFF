#!/usr/bin/env python

import os, numpy as np
from molmod.units import *
from cfit2 import *

from mfit2.program import DefaultProgram as MFitDefaultProgram
from mfit2.terms import HarmonicTerm, CosineTerm, GaussianChargeTerm
from mfit2.yaff import dump_yaff_parameters as dump_myaff
from mfit2.lsq import HessianLSQ
from mfit2.ics import dist, angle, dihed

from zfit2.program import DefaultProgram as ZFitDefaultProgram
from zfit2.terms import SplitChargeTerm, FixedBackgroundChargeTerm
from zfit2.yaff import dump_yaff_parameters as dump_zyaff
from zfit2.lsq import ChargeLSQ

from fftable import FFTable

__all__ = ['ZFitProgram', 'MFitProgram']

def generate_mterms(system, icnames):
    print 'MFIT2  TERMS: auto generate harmonic force field terms'
    terms = []    
    for icname in icnames:
        ics = system.ics[icname]
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
                if   6 in [n1, n2]: m = 4
                elif 5 in [n1, n2]: m = 1
                elif 4 in [n1, n2]: m = 3
                elif 3 in [n1, n2]: m = 2
                elif 2 in [n1, n2]: m = 1
                else: raise ValueError('Dihedral %i,%i,%i,%i has no atoms bonded to central pair' %(ic.indexes[0], ic.indexes[1], ic.indexes[2], ic.indexes[3]))
            terms.append(CosineTerm(DihedralFilter(types[0], types[1], types[2], types[3]), dihed, m, 0.0))
        else:
            return ValueError('Invalied IC type, recieved %s' %kind)
    if system.eirule>-1:
        terms.append( GaussianChargeTerm(PairFilter('*', '*', exclude_bonded=system.eirule), dist)  )
    return terms


class ZFitProgram(ZFitDefaultProgram):
    def __init__(self, system):
        dn_out = '%s/out/zfit' %(os.getcwd())
        os.system('mkdir -p %s' %dn_out)
        log._file = open('%s/log.txt' %dn_out, 'w')
        terms = [TermGenerator(SplitChargeTerm, BondFilter('*', '*', allow_homonuclear=False))]
        terms.append(FixedBackgroundChargeTerm())
        rules = []
        lsqs  = [ChargeLSQ('AQ', 'int-system.chk', 1.0)]
        fn_pars_init = None
        self.scales = [1.0, 1.0, 1.0]
        if system.eirule > 0 or system.eirule==-1: self.scales[0] = 0.0
        if system.eirule > 1 or system.eirule==-1: self.scales[1] = 0.0
        if system.eirule > 2 or system.eirule==-1: self.scales[2] = 0.0
        ZFitDefaultProgram.__init__(self, dn_out, terms, rules, lsqs, fn_pars_init, fix_patterns=[], free_patterns=[],
                                step_dl_min=1e-5, step_dl_max=1e2, max_iter=50, weight_mode=None, weight_ridge=0.1)
    
    def run(self):
        print 'ZFIT2  RUN  : calculate split charges using FFit2'
        model = ZFitDefaultProgram.run(self)
        dump_zyaff(model.terms, model.training_set, 'pars_zyaff.txt', scales=self.scales)
        model.dump_pars('pars_zfit2.txt')
        return model.training_set[0].ac


class MFitProgram(MFitDefaultProgram):
    def __init__(self, system, icnames):
        dn_out = '%s/out/mfit' %(os.getcwd())
        os.system('mkdir -p %s' %dn_out)
        log._file = open('%s/log.txt' %dn_out, 'w')
        terms = generate_mterms(system, icnames)
        rules = []
        lsqs = [HessianLSQ('Hess', 'int-system.chk', 1.0)]
        fn_pars_init = 'int-pars.txt'
        fix_patterns = ['q0']
        MFitDefaultProgram.__init__(self, dn_out, terms, rules, lsqs, fn_pars_init, fix_patterns=fix_patterns, free_patterns=[],
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
        print 'MFIT2  RUN  : refine force field parameters using FFit2'
        model = MFitDefaultProgram.run(self)
        dump_myaff(model.terms, model.training_set, 'pars_myaff.txt')
        if os.path.isfile('pars_zyaff.txt'):
            os.system('cat pars_zyaff.txt pars_myaff.txt > pars_yaff.txt')
            os.system('rm pars_zyaff.txt pars_myaff.txt')
        else:
            os.system('mv pars_myaff.txt pars_yaff.txt')
        model.dump_pars('pars_mfit2.txt')
        fftab = FFTable.from_ffit2(model)
        return fftab
