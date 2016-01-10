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

from yaff.pes.ff import ForcePartValence
from yaff.pes.parameters import *
from yaff.pes.vlist import ValenceList
from yaff.pes.vlist import Harmonic, Fues, Cosine, Cross
from yaff.pes.iclist import InternalCoordinateList
from yaff.pes.iclist import Bond, BendAngle, DihedAngle, OopDist
from yaff.pes.dlist import DeltaList

import numpy as np

__all__ = ['ValenceFF']

class Term(object):
    '''
        A class to store some extra easy-accessible information compared to
        what is stored in the Yaff classes.
    '''
    def __init__(self, index, kind, basename, tasks, ics, units, master=None, slaves=None):
        self.index = index
        self.kind = kind
        self.basename = basename
        self.tasks = tasks
        self.ics = ics
        self.units = units
        self.master = master
        self.slaves = slaves
    
    def is_master(self):
        return self.index==self.master
    
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
        line = '\t'+self.basename+' '*(40-len(self.basename))
        for start, ipar, strfmt in fmt:
            par_line = start
            if not np.isnan(means[ipar]):
                par_line += strfmt %(means[ipar]/parse_unit(self.units[ipar]))
                par_line += u' \u00B1 '
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
        self.data = []
        ForcePartValence.__init__(self, system)
        self.init_terms(specs=specs)
        self.set_dihedral_multiplicities()
    
    def add_data(self, kind, basename, tasks, ics, units):
        index = len(self.data)
        master = None
        slaves = None
        for i, term in enumerate(self.data):
            if term.basename==basename:
                if term.is_master():
                    master = term.index
                    term.slaves.append(index)
                else:
                    assert master==term.master
        if master is None:
            master = index
            slaves = []
        term = Term(
            index, kind, basename, tasks, ics, 
            units, master=master, slaves=slaves
        )
        self.data.append(term)
    
    def iter_masters(self, label=None):
        for term in self.data:
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
            types = [ffatypes[i] for i in bond]
            if types[-1]<types[0]:
                types = types[::-1]
                bond = bond[::-1]
            ic = Bond(bond[0], bond[1])
            basename = 'BondHarm/%s.%s' %(types[0],types[1])
            units = ['kjmol/A**2', 'A']
            self.add_data(0, basename, ['PT_ALL', 'HC_FC_DIAG'], [ic], units)
            self.add_term(Harmonic(None, None, ic))
        #get the angle terms
        for angle in self.system.iter_angles():
            types = [ffatypes[i] for i in angle]
            if types[-1]<types[0]:
                types = types[::-1]
                angle = angle[::-1]
            ic = BendAngle(angle[0], angle[1], angle[2])
            basename = 'BendAHarm/%s.%s.%s' %(types[0],types[1],types[2])
            units = ['kjmol/rad**2', 'deg']
            self.add_data(0,basename, ['PT_ALL', 'HC_FC_DIAG'], [ic], units)
            self.add_term(Harmonic(None, None, ic))
        #get the dihedral terms
        for dihedral in self.system.iter_dihedrals():
            types = [ffatypes[i] for i in dihedral]
            if types[-1]<types[0]:
                types = types[::-1]
                dihedral = dihedral[::-1]
            ic = DihedAngle(dihedral[0], dihedral[1], dihedral[2], dihedral[3])
            basename = 'Torsion/%s.%s.%s.%s' %(types[0],types[1],types[2],types[3])
            units = ['au', 'kjmol', 'deg']
            self.add_data(4, basename, ['EQ_RV', 'HC_FC_DIAG'], [ic], units)
            self.add_term(Cosine(None, None, None, ic))
        #get the out-of-plane distance terms
        for oopdist in self.system.iter_oops():
            types = sorted([ffatypes[i] for i in oopdist[:3]])
            types.append(ffatypes[oopdist[3]])
            ic = OopDist(oopdist[0], oopdist[1], oopdist[2], oopdist[3])
            basename = 'Oopdist/%s.%s.%s.%s' %(types[0],types[1],types[2],types[3])
            units = ['kjmol/A**2', 'A']
            self.add_data(0, basename, ['PT_ALL', 'HC_FC_DIAG'], [ic], units)
            self.add_term(Harmonic(None, None, ic))
    
    def _get_bond(self, atom0, atom1):
        '''
            Method to extract the Bond IC between the given atoms from the
            existing ICList. Return error if it does not yet exist.
        '''
        bond = Bond(atom0, atom1)
        rows_signs = bond.get_rows_signs(self.dlist)
        key = (bond.kind,) + sum(rows_signs, ())
        if key not in self.iclist.lookup.keys():
            bond = Bond(atom1, atom0)
            rows_signs = bond.get_rows_signs(self.dlist)
            key = (bond.kind,) + sum(rows_signs, ())
            assert key in self.iclist.lookup.keys(), \
                'No bond between atoms %i and %i detected' %(atom0, atom1)
        return bond
    
    def init_cross_terms(self, specs=None):
        ffatypes = [self.system.ffatypes[i] for i in self.system.ffatype_ids]
        for angle in self.system.iter_angles():
            types = [ffatypes[i] for i in angle]
            if types[-1]<types[0]:
                types = types[::-1]
                angle = angle[::-1]
            bond0 = self._get_bond(angle[0], angle[1])
            bond1 = self._get_bond(angle[1], angle[2])
            bend = BendAngle(angle[0], angle[1], angle[2])
            #add stretch-stretch
            basename = 'CrossBondBond/%s.%s.%s' %(types[0],types[1],types[2])
            units = ['kjmol/A**2', 'A', 'A']
            self.add_data(3, basename, ['HC_FC_CROSS'], [bond0,bond1], units)
            self.add_term(Cross(None, None, None, bond0, bond1))
            #add stretch0-bend
            basename = 'CrossBond0Bend/%s.%s.%s' %(types[0],types[1],types[2])
            units = ['kjmol/(A*rad)', 'A', 'deg']
            self.add_data(3, basename, ['HC_FC_CROSS'], [bond0, bend], units)
            self.add_term(Cross(None, None, None, bond0, bend))
            #add stretch1-bend
            basename = 'CrossBond1Bend/%s.%s.%s' %(types[0],types[1],types[2])
            units = ['kjmol/(A*rad)', 'A', 'deg']
            self.add_data(3, basename, ['HC_FC_CROSS'], [bond1, bend], units)
            self.add_term(Cross(None, None, None, bond1, bend))

    def set_dihedral_multiplicities(self):  
        for index in xrange(self.vlist.nv):
            term = self.vlist.vtab[index]
            ic = self.iclist.ictab[term['ic0']]
            if not ic['kind']==4: continue #only proceed for DihedAngle
            i = self.data[index].ics[0].index_pairs[0][0]
            j = self.data[index].ics[0].index_pairs[1][1]
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
                    %(self.data[index].basename, index)
            self.set_params(index, m=m)
    
    def calc_energy(self, pos):
        self.system.pos = pos.copy()
        self.dlist.forward()
        self.iclist.forward()
        self.vlist.forward()
        return self.compute()
    
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

    def dump_screen(self):
        sequence = ['bondharm', 'bendaharm', 'torsion', 'oopdist', 'cross']
        print ''
        for label in sequence:
            lines = []
            for term in self.iter_masters(label=label):
                lines.append(term.to_string(self))
            for line in sorted(lines):
                print line
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
        pars = ParameterDefinition('PARS')
        for i, master in enumerate(self.iter_masters(label=prefix)):
            ffatypes = master.basename.split('/')[1].split('.')
            K, rv0, rv1 = self.get_params(master.index)
            if K<1e-6: continue
            if 'bondbond' in master.basename.lower():
                pars.lines.append(
                    '%8s  %8s  %8s  %.10e  %.10e  %.10e  %.10e  %.10e  %.10e' %(
                        ffatypes[0], ffatypes[1], ffatypes[2],
                        K/parse_unit(master.units[0]), 0.0, 0.0,
                        rv0/parse_unit(master.units[1]),
                        rv1/parse_unit(master.units[2]), 0.0,
                ))
            elif 'bond0bend' in master.basename.lower():
                pars.lines.append(
                    '%8s  %8s  %8s  %.10e  %.10e  %.10e  %.10e  %.10e  %.10e' %(
                        ffatypes[0], ffatypes[1], ffatypes[2],
                        0.0, K/parse_unit(master.units[0]), 0.0,
                        rv0/parse_unit(master.units[1]),
                        0.0, rv1/parse_unit(master.units[2]),
                ))
            elif 'bond1bend' in master.basename.lower():
                pars.lines.append(
                    '%8s  %8s  %8s  %.10e  %.10e  %.10e  %.10e  %.10e  %.10e' %(
                        ffatypes[0], ffatypes[1], ffatypes[2],
                        0.0, 0.0, K/parse_unit(master.units[0]),
                        0.0, rv0/parse_unit(master.units[1]),
                        rv1/parse_unit(master.units[2]),
                ))
            else:
                raise ValueError
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
