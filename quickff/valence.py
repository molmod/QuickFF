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

from yaff.pes.ff import ForceField, ForcePartValence
from yaff.pes.parameters import *
from yaff.pes.vlist import ValenceList
from yaff.pes.vlist import Harmonic, Fues, Cosine, Cross
from yaff.pes.iclist import InternalCoordinateList
from yaff.pes.iclist import Bond, BendAngle, DihedAngle, OopDist
from yaff.pes.dlist import DeltaList
from yaff.sampling.harmonic import estimate_cart_hessian

from quickff.tools import dihed_to_chebychev, term_sort_atypes

import numpy as np

__all__ = ['ValenceFF']

class Term(object):
    '''
        A class to store easy-accessible information about a term included in 
        the valence force field
    '''
    def __init__(self, index, basename, kind, ics, tasks, units,master=None, slaves=None):
        self.index = index
        self.basename = basename
        self.kind = kind
        self.ics = ics
        self.tasks = tasks
        self.units = units
        self.master = master
        self.slaves = slaves
    
    def is_master(self):
        return self.master==self.index
    
    def to_string(self, valence):
        #check if current is master
        assert self.master is None or self.is_master(), \
            'Current term is not the master'
        #collect parameters
        npars = len(valence.get_params(self.index))
        pars = np.zeros([len(self.slaves)+1, npars], float)
        pars[0,:] = np.array(valence.get_params(self.index))
        for i, index in enumerate(self.slaves):
            pars[1+i,:] = np.array(valence.get_params(index))
        means = np.zeros(npars, float)
        stds = np.zeros(npars, float)
        for ipar in xrange(npars):
            means[ipar] = pars[:,ipar].mean()
            stds[ipar] = pars[:,ipar].std()
        #fmt is a list of tuples of the form: (name, index, string format)
        if self.kind==0:#harmonic
            fmt = [('fc = ', 0, '%6.1f'), ('rv  = ', 1, '%7.3f')]
        elif self.kind==3:#cross
            fmt = [('fc = ', 0, '%6.1f'), ('rv0 = ', 1, '%7.3f'), ('rv1 = ', 2, '%7.3f')]
        elif self.kind==4:#cosine
            fmt = [('fc = ', 1, '%6.1f'), ('rv  = ', 2, '%7.3f'), ('m = ', 0, '%i')]
        else:
            raise NotImplementedError
        line = self.basename+' '*(40-len(self.basename))
        for start, ipar, strfmt in fmt:
            par_line = start
            if not np.isnan(means[ipar]):
                par_line += strfmt %(means[ipar]/parse_unit(self.units[ipar]))
                par_line += '+-'#u' \u00B1 '
                par_line += strfmt %(stds[ipar]/parse_unit(self.units[ipar]))
                par_line += ' %s' %self.units[ipar]
            line += par_line + ' '*(35-len(par_line))
        return line


class ValenceFF(ForcePartValence):
    '''
        Class to collect all valence terms in the force field for which
        parameters need to be estimated.
    '''
    def __init__(self, system, specs=None):
        self.system = system
        self.terms = []
        ForcePartValence.__init__(self, system)
        self.init_terms(specs=specs)
        self.set_dihedral_multiplicities()
    
    def add_term(self, termpot, ics, atypes, tasks, units):
        '''
            Adds new term both to the Yaff vlist object and a new QuickFF
            list (self.terms) which holds all information about the term
            for easy access later in QuickFF.
        '''
        index = len(self.terms)
        #define the name
        if len(ics)==1:
            tmp = {
                (0,0): 'BondHarm/', (2,0): 'BendAHarm/', (4,4): 'Torsion/',
                (10,0): 'Oopdist/',
            }
            prefix = tmp[(ics[0].kind, termpot.kind)]
            suffix = ''
        else:
            assert len(ics)==2 and termpot.kind==3
            prefix = 'Cross/'
            if ics[0].kind==0 and ics[1].kind==0:
                suffix = '/bb' #first bond and second bond
            elif ics[0].kind==0 and ics[1].kind==2:
                if set(ics[0].index_pairs[0])==set(ics[1].index_pairs[0]):
                    suffix = '/b0a' #first bond and angle
                elif set(ics[0].index_pairs[0])==set(ics[1].index_pairs[1]):
                    suffix = '/b1a' #second bond and angle
                else:   
                    raise ValueError('Incompatible bond/angle given in cross term')
            else:
                raise ValueError('Incompatible ICs given in cross term')
        sorted_atypes = atypes[:]
        if len(ics)==1 and ics[0].kind==10:
            assert len(atypes)==4
            sorted_atypes = sorted(atypes[:3])+[atypes[3]]
        else:
            if atypes[0]>atypes[-1]:
                sorted_atypes = atypes[::-1]
        basename = prefix+'.'.join(sorted_atypes)+suffix
        #search for possible master and update slaves
        master = None
        slaves = None
        for i, term in enumerate(self.terms):
            if term.basename==basename:
                if term.is_master():
                    master = term.index
                    term.slaves.append(index)
                else:
                    assert master==term.master
        if master is None:
            master = index
            slaves = []
        #add term to self.terms and self.vlist.vtab
        self.terms.append(Term(
            index, basename, termpot.kind, ics, tasks,
            units, master=master, slaves=slaves
        ))
        _npars = {0: 2, 3: 3, 4: 3}#number of pars for each term kind
        args = [None,]*_npars[termpot.kind] + ics
        ForcePartValence.add_term(self, termpot(*args))
        assert len(self.terms)==self.vlist.nv
    
    def iter_masters(self, label=None):
        '''
            Iterate over all master terms (whos name possibly contain the given
            label) in the valence force field
        '''
        for term in self.terms:
            if label is None or label.lower() in term.basename.lower():
                if term.is_master():
                    yield term
    
    def init_terms(self, specs=None):
        '''
            specs
                A list of specifications to control the behaviour of certain
                terms for which parameters need to be estimated. The tasks
                of a term can be defined in this object.
        '''
        #TODO: include specs functionality
        if specs is not None:
            raise NotImplementedError('Specs functionality not yet implemented')
        ffatypes = [self.system.ffatypes[fid] for fid in self.system.ffatype_ids]
        #get the bond terms
        for bond in self.system.iter_bonds():
            bond, types = term_sort_atypes(ffatypes, bond, 'bond')
            units = ['kjmol/A**2', 'A']
            self.add_term(Harmonic, [Bond(*bond)], types, ['PT_ALL', 'HC_FC_DIAG'], units)
        #get the angle terms
        for angle in self.system.iter_angles():
            angle, types = term_sort_atypes(ffatypes, angle, 'angle')
            units = ['kjmol/rad**2', 'deg']
            self.add_term(Harmonic, [BendAngle(*angle)], types, ['PT_ALL', 'HC_FC_DIAG'], units)
        #get the dihedral terms
        for dihedral in self.system.iter_dihedrals():
            dihedral, types = term_sort_atypes(ffatypes, dihedral, 'dihedral')
            units = ['au', 'kjmol', 'deg']
            self.add_term(Cosine, [DihedAngle(*dihedral)], types, ['EQ_RV', 'HC_FC_DIAG'], units)
        #get the out-of-plane distance terms
        for oopdist in self.system.iter_oops():
            oopdist, types = term_sort_atypes(ffatypes, oopdist, 'dihedral')
            units = ['kjmol/A**2', 'A']
            self.add_term(Harmonic, [OopDist(*oopdist)], types, ['PT_ALL', 'HC_FC_DIAG'], units)
    
    def init_cross_terms(self, specs=None):
        ffatypes = [self.system.ffatypes[i] for i in self.system.ffatype_ids]
        for angle in self.system.iter_angles():
            angle, types = term_sort_atypes(ffatypes, angle, 'angle')
            bond0, btypes = term_sort_atypes(ffatypes, angle[:2], 'bond')
            bond1, btypes = term_sort_atypes(ffatypes, angle[1:], 'bond')
            #add stretch-stretch
            self.add_term(Cross,
                [Bond(*bond0), Bond(*bond1)],
                types, ['HC_FC_CROSS'], ['kjmol/A**2', 'A', 'A']
            )
            #add stretch0-bend
            self.add_term(Cross, 
                [Bond(*bond0), BendAngle(*angle)],
                types, ['HC_FC_CROSS'], ['kjmol/(A*rad)', 'A', 'deg']
            )
            #add stretch1-bend
            self.add_term(Cross, 
                [Bond(*bond1), BendAngle(*angle)],
                types, ['HC_FC_CROSS'], ['kjmol/(A*rad)', 'A', 'deg']
            )

    def set_dihedral_multiplicities(self):
        'Estimate the dihedral multiplicities from the local topology'
        for index in xrange(self.vlist.nv):
            term = self.vlist.vtab[index]
            ic = self.iclist.ictab[term['ic0']]
            if not ic['kind']==4: continue #only proceed for DihedAngle
            i = self.terms[index].ics[0].index_pairs[0][0]
            j = self.terms[index].ics[0].index_pairs[1][1]
            n1 = len(self.system.neighs1[i])
            n2 = len(self.system.neighs1[j])
            if set([n1,n2])==set([4,4]):   m = 3
            elif set([n1,n2])==set([3,4]): m = 6
            elif set([n1,n2])==set([2,4]): m = 3
            elif set([n1,n2])==set([3,3]): m = 2
            elif set([n1,n2])==set([2,3]): m = 2
            elif set([n1,n2])==set([2,2]): m = 1
            else:
                raise ValueError('No estimate for mult of %s(%i) found') \
                    %(self.terms[index].basename, index)
            self.set_params(index, m=m)
    
    def set_dihedral_rest_values(self):
        'Estimate the dihedral rest values from the local topology'
        #set the rest value for each master and slave seperately
        for index in xrange(self.vlist.nv):
            term = self.terms[index]
            vterm = self.vlist.vtab[index]
            ic = self.iclist.ictab[vterm['ic0']]
            if not ic['kind']==4: continue #only proceed for DihedAngle
            m = self.get_params(index, only='m')
            assert not np.isnan(m), 'Multiplicity not set for %s' %term.basename
            per = 360.0*deg/m
            psi0 = ic['value']%per
            if (psi0>=0 and psi0<=per/6.0) or (psi0>=5*per/6.0 and psi0<=per):
                rv = 0.0
            elif psi0>=2*per/6.0 and psi0<=4*per/6.0:
                rv = per/2.0
            else:
                raise ValueError('Could not determine rv of %s' %term.basename)
            self.set_params(index, rv0=rv)
        #check whether all slaves have identical rest value as the master
        for master in self.iter_masters(label='Torsion'):
            rv = self.get_params(master.index, only='rv')
            for slave in master.slaves:
                rvslave = self.get_params(slave, only='rv')
                assert rvslave==rv, "Slaves of %s do not have identical \
                    rest value as the master" %master.basename
    
    def calc_energy(self, pos):
        self.system.pos = pos.copy()
        self.dlist.forward()
        self.iclist.forward()
        self.vlist.forward()
        return self.compute()
    
    def get_hessian_contrib(self, index, fc=None):
        '''
            Get the contribution to the covalent hessian of term with given
            index (and its slaves). If fc is given, set the fc of the master
            and its slave to the given fc.
        '''
        val = ForcePartValence(self.system)
        kind = self.vlist.vtab[index]['kind']
        masterslaves = [index]+self.terms[index].slaves
        kind_to_term = {0: Harmonic, 2: Fues, 3: Cross, 4: Cosine}
        if kind==4:#Cosine
            m, k, rv = self.get_params(index)
            if fc is not None: k = fc
            for jterm in masterslaves:
                ics = self.terms[jterm].ics
                chebychev = None #dihed_to_chebychev([m, k, rv], ics[0])
                if chebychev is not None:
                    val.add_term(chebychev)
                else:
                    args = (m, k, rv) + tuple(ics)
                    val.add_term(Cosine(*args))
        elif kind==3:#cross
            k, rv0, rv1 = self.get_params(index)
            if fc is not None: k = fc
            for jterm in masterslaves:
                ics = self.terms[jterm].ics
                args = (k, rv0, rv1) + tuple(ics)
                val.add_term(Cross(*args))
        elif kind==0:#Harmonic:
            k, rv = self.get_params(index)
            if fc is not None: k = fc
            for jterm in masterslaves:
                ics = self.terms[jterm].ics
                args = (k, rv) + tuple(ics)
                val.add_term(kind_to_term[kind](*args))
        else:
            raise ValueError('Term kind %i not supported' %kind)
        ff = ForceField(self.system, [val])
        hcov = estimate_cart_hessian(ff)
        return hcov
    
    def set_params(self, term_index, fc=None, rv0=None, rv1=None, m=None):
        term = self.vlist.vtab[term_index]
        if term['kind'] in [0,2]:#['Harmonic', 'Fues']
            if fc is not None:  term['par0'] = fc
            if rv0 is not None: term['par1'] = rv0
        elif term['kind'] in [4]:#['Cosine']
            if m is not None:   term['par0'] = m
            if fc is not None:  term['par1'] = fc
            if rv0 is not None: term['par2'] = rv0
        elif term['kind'] in [3]:#['Cross']
            if fc is not None:  term['par0'] = fc
            if rv0 is not None: term['par1'] = rv0
            if rv1 is not None: term['par2'] = rv1
        else:
            raise NotImplementedError, \
                'set_params not implemented for Yaff %s term' %term['kind']
    
    def get_params(self, term_index, only='all'):
        term = self.vlist.vtab[term_index]
        if term['kind'] in [0,2]:#['Harmonic', 'Fues']
            if only.lower()=='all': return term['par0'], term['par1']
            elif only.lower()=='fc': return term['par0']
            elif only.lower()=='rv': return term['par1']
            else: raise ValueError('Invalid par kind definition %s' %only)
        elif term['kind'] in [4]:#['Cosine']
            if only.lower()=='all': return term['par0'], term['par1'], term['par2']
            elif only.lower()=='m': return term['par0']
            elif only.lower()=='fc': return term['par1']
            elif only.lower()=='rv': return term['par2']
            else: raise ValueError('Invalid par kind definition %s' %only)
        elif term['kind'] in [3]:#['Cross']
            if only.lower()=='all': return term['par0'], term['par1'], term['par2']
            elif only.lower()=='fc': return term['par0']
            elif only.lower()=='rv0': return term['par1']
            elif only.lower()=='rv1': return term['par2']
            else: raise ValueError('Invalid par kind definition %s' %only)
        else:
            raise NotImplementedError, \
                'set_params not implemented for Yaff %s term' %term['kind']

    def check_params(self, term, labels):
        '''
            Check whether the given term has all given pars defined in
            labels.
            
            **Arguments**
            
            term
                An instance of the Term class defining the term that has to be
                checked
            
            labels
                A list of strings defining which parameters should be checked.
                only arguments of the `only` option of Valence.get_params
                are allowed.
        '''
        for label in labels:
            value = self.get_params(term.index, only=label)
            assert not np.isnan(value), '%s of %s is not set' %(label, term.basename)
    
    def dump_screen(self):
        'Dump the parameters to stdout'
        sequence = ['bondharm', 'bendaharm', 'torsion', 'oopdist', 'cross']
        print ''
        for label in sequence:
            lines = []
            for term in self.iter_masters(label=label):
                lines.append(term.to_string(self))
            for line in sorted(lines):
                print '  '+line
        print ''

    def _bonds_to_yaff(self):
        'construct a bonds section of a yaff parameter file'
        prefix = 'BONDHARM'
        units = ParameterDefinition('UNIT', lines=['K kjmol/A**2', 'R0 A'])
        pars = ParameterDefinition('PARS')
        for i, master in enumerate(self.iter_masters(label=prefix)):
            ffatypes = master.basename.split('/')[1].split('.')
            K, q0 = self.get_params(master.index)
            if K<1e-6: continue
            pars.lines.append('%8s  %8s  %.10e  %.10e' %(
                ffatypes[0], ffatypes[1],
                K/parse_unit(master.units[0]), q0/parse_unit(master.units[1])
            ))
        return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})

    def _bends_to_yaff(self):
        'construct a bends section of a yaff parameter file'
        prefix = 'BENDAHARM'
        units = ParameterDefinition('UNIT', lines=['K kjmol/rad**2', 'THETA0 deg'])
        pars = ParameterDefinition('PARS')
        for i, master in enumerate(self.iter_masters(label=prefix)):
            ffatypes = master.basename.split('/')[1].split('.')
            K, q0 = self.get_params(master.index)
            if K<1e-6: continue
            pars.lines.append('%8s  %8s  %8s  %.10e  %.10e' %(
                ffatypes[0], ffatypes[1], ffatypes[2],
                K/parse_unit(master.units[0]), q0/parse_unit(master.units[1])
            ))
        return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})

    def _diheds_to_yaff(self):
        'construct a dihedral section of a yaff parameter file'
        prefix = 'TORSION'
        units = ParameterDefinition('UNIT', lines=['A kjmol', 'PHI0 deg'])
        pars = ParameterDefinition('PARS')
        for i, master in enumerate(self.iter_masters(label=prefix)):
            ffatypes = master.basename.split('/')[1].split('.')
            m, K, q0 = self.get_params(master.index)
            if K<1e-6: continue
            pars.lines.append('%8s  %8s  %8s  %8s  %1i %.10e  %.10e' %(
                ffatypes[0], ffatypes[1],  ffatypes[2], ffatypes[3], m,
                K/parse_unit(master.units[1]), q0/parse_unit(master.units[2])
            ))
        return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})

    def _opdists_to_yaff(self):
        'construct a opdist section of a yaff parameter file'
        prefix = 'OOPDIST'
        units = ParameterDefinition('UNIT', lines=['K kjmol/A**2', 'D0 A'])
        pars = ParameterDefinition('PARS')
        for i, master in enumerate(self.iter_masters(label=prefix)):
            ffatypes = master.basename.split('/')[1].split('.')
            K, q0 = self.get_params(master.index)
            if K<1e-6: continue
            pars.lines.append('%8s  %8s  %8s  %8s  %.10e  %.10e' %(
                ffatypes[0], ffatypes[1], ffatypes[2], ffatypes[3],
                K/parse_unit(master.units[0]), q0/parse_unit(master.units[1])
            ))
        return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})

    def _cross_to_yaff(self):
        'construct a cross section of a yaff parameter file'
        prefix = 'CROSS'
        units = ParameterDefinition(
            'UNIT', 
            lines=[
                'KSS kjmol/angstrom**2', 'KBS0 kjmol/(angstrom*rad)',
                'KBS1 kjmol/(angstrom*rad)', 'R0 angstrom', 'R1 angstrom',
                'THETA0 deg'
            ]
        )
        done = []
        pars = ParameterDefinition('PARS')
        for i, master in enumerate(self.iter_masters(label=prefix)):
            prefix, ffatypes, suffix = master.basename.split('/')
            label = prefix+'/'+ffatypes
            if label in done: continue
            for j, other in enumerate(self.iter_masters(label=label)):
                if 'bb' in other.basename:
                    bb = self.get_params(other.index)
                    kssunit = other.units[0]
                    r0unit = other.units[1]
                    r1unit = other.units[2]
                elif 'b0a' in other.basename:
                    b0a = self.get_params(other.index)
                    kbs0unit = other.units[0]
                    theta0unit = other.units[2]
                elif 'b1a' in other.basename:
                    b1a = self.get_params(other.index)
                    kbs1unit = other.units[0]
                else:
                    raise ValueError('Invalid Cross term %s' %other.basename)
            assert j==2, 'Exactly 3 %s terms should exist' %label
            assert bb[1]==b0a[1], 'Incompatible parameters in %s' %label
            assert bb[2]==b1a[1], 'Incompatible parameters in %s' %label
            assert b0a[2]==b1a[2], 'Incompatible parameters in %s' %label
            Kss, r0, r1 = bb
            Kbs0, r0, theta0 = b0a
            Kbs1, r1, theta0 = b1a
            ffatypes = ffatypes.split('.')
            pars.lines.append(
                '%8s  %8s  %8s  %.10e  %.10e  %.10e  %.10e  %.10e  %.10e' %(
                    ffatypes[0], ffatypes[1], ffatypes[2],
                    Kss/parse_unit(kssunit), Kbs0/parse_unit(kbs0unit),
                    Kbs1/parse_unit(kbs1unit), r0/parse_unit(r0unit),
                    r1/parse_unit(r1unit), theta0/parse_unit(theta0unit),
            ))
            done.append(label)
        return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})

    def dump_yaff(self, fn):
        sections = [
            self._bonds_to_yaff(), self._bends_to_yaff(), 
            self._diheds_to_yaff(), self._opdists_to_yaff(),
            self._cross_to_yaff(),
        ]
        f = open(fn, 'w')
        for section in sections:
            print >> f, '# %s' %section.prefix
            print >> f, '#-%s' %('-'*len(section.prefix))
            for line in section['UNIT'].lines:
                print >> f, '%s:UNIT  %s' %(section.prefix, line)
            print >> f, ''
            for line in section['PARS'].lines:
                print >> f, '%s:PARS  %s' %(section.prefix, line)
            print >> f, ''
            print >> f, ''
        f.close()
