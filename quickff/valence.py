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
from yaff.pes.vlist import Harmonic, PolyFour, Fues, Cosine, Cross
from yaff.pes.iclist import InternalCoordinateList
from yaff.pes.iclist import Bond, BendAngle, DihedCos, DihedAngle, OopDist, SqOopDist
from yaff.pes.dlist import DeltaList
from yaff.sampling.harmonic import estimate_cart_hessian

from quickff.tools import term_sort_atypes, get_multiplicity, get_restvalue

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
        #convert term pars to string
        line = self.basename+' '*(40-len(self.basename))
        if self.kind==0:#harmonic
            fc, rv = pars.mean(axis=0)
            dfc, drv = pars.std(axis=0)
            par_line = 'fc = %7.1f +- %5.1f %s' %(
                fc/parse_unit(self.units[0]), dfc/parse_unit(self.units[0]),
                self.units[0]
            )
            line += par_line + ' '*(35-len(par_line))
            par_line = 'rv  = %7.3f +- %7.3f %s' %(
                rv/parse_unit(self.units[1]), drv/parse_unit(self.units[1]),
                self.units[1]
            )
            line += par_line
        elif self.kind==1 and self.ics[0].kind==3:#PolyFour for torsc2harm
            fcs = 0.5*pars[:,3]
            rvs = np.arccos(np.sqrt(-pars[:,1]/(2*pars[:,3])))
            fc, dfc, rv, drv = fcs.mean(), fcs.std(), rvs.mean(), rvs.std()
            par_line = 'fc = %7.1f +- %5.1f %s' %(
                fc/parse_unit(self.units[3]), dfc/parse_unit(self.units[3]),
                self.units[3]
            )
            line += par_line + ' '*(35-len(par_line))
            par_line = 'rv  = %7.3f +- %7.3f %s' %(
                rv/deg, drv/deg, 'deg'
            )
            line += par_line
        elif self.kind==3:#cross
            fc, rv0, rv1 = pars.mean(axis=0)
            dfc, drv0, drv1 = pars.std(axis=0)
            par_line = 'fc = %7.1f +- %5.1f %s' %(
                fc/parse_unit(self.units[0]),
                dfc/parse_unit(self.units[0]),
                self.units[0]
            )
            line += par_line + ' '*(35-len(par_line))
            par_line = 'rv0 = %7.3f +- %7.3f %s' %(
                rv0/parse_unit(self.units[1]),
                drv0/parse_unit(self.units[1]),
                self.units[1]
            )
            line += par_line + ' '*(35-len(par_line))
            par_line = 'rv1 = %7.3f +- %7.3f %s' %(
                rv1/parse_unit(self.units[2]),
                drv1/parse_unit(self.units[2]),
                self.units[2]
            )
            line += par_line
        elif self.kind==4:#cosine
            m, fc, rv = pars.mean(axis=0)
            dm, dfc, drv = pars.std(axis=0)
            par_line = 'fc = %7.1f +- %5.1f %s' %(
                fc/parse_unit(self.units[1]),
                dfc/parse_unit(self.units[1]),
                self.units[1]
            )
            line += par_line + ' '*(35-len(par_line))
            par_line = 'rv  = %7.3f +- %7.3f %s' %(
                rv/parse_unit(self.units[2]),
                drv/parse_unit(self.units[2]),
                self.units[2]
            )
            line += par_line + ' '*(35-len(par_line))
            par_line = 'm = %1i +- %.1f' %(
                m/parse_unit(self.units[0]),
                dm/parse_unit(self.units[0]),
            )
            line += par_line
        else:
            raise NotImplementedError
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
        self.init_bond_terms()
        self.init_bend_terms()
        self.init_dihedral_terms()
        self.init_oop_terms()
    
    def add_term(self, pot, ics, atypes, tasks, units):
        '''
            Adds new term both to the Yaff vlist object and a new QuickFF
            list (self.terms) which holds all information about the term
            for easy access later in QuickFF.
        '''
        index = len(self.terms)
        #define the name
        if len(ics)==1:
            tmp = {
                (0,0) : 'BondHarm/', (2,0) : 'BendAHarm/' ,
                (4,4) : 'Torsion/' , (3,1) : 'TorsC2Harm/',
                (10,0): 'Oopdist/' , (11,0): 'SqOopdist/' ,
            }
            prefix = tmp[(ics[0].kind, pot.kind)]
            suffix = ''
        else:
            assert len(ics)==2 and pot.kind==3
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
        basename = prefix+'.'.join(atypes)+suffix
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
        term = Term(
            index, basename, pot.kind, ics, tasks,
            units, master=master, slaves=slaves
        )
        self.terms.append(term)
        if pot.kind==1:#all 4 parameters of PolyFour are given as 1 tuple
            args = [(None,)*len(units)] + ics
        else:
            args = [None,]*len(units) + ics
        ForcePartValence.add_term(self, pot(*args))
        return term
    
    def iter_masters(self, label=None):
        '''
            Iterate over all master terms (whos name possibly contain the given
            label) in the valence force field
        '''
        for term in self.terms:
            if label is None or label.lower() in term.basename.lower():
                if term.is_master():
                    yield term
    
    def init_bond_terms(self):
        ffatypes = [self.system.ffatypes[fid] for fid in self.system.ffatype_ids]
        #get the bond terms
        for bond in self.system.iter_bonds():
            bond, types = term_sort_atypes(ffatypes, bond, 'bond')
            units = ['kjmol/A**2', 'A']
            self.add_term(Harmonic, [Bond(*bond)], types, ['PT_ALL', 'HC_FC_DIAG'], units)

    def init_bend_terms(self):
        ffatypes = [self.system.ffatypes[fid] for fid in self.system.ffatype_ids]
        #get the angle terms
        for angle in self.system.iter_angles():
            angle, types = term_sort_atypes(ffatypes, angle, 'angle')
            units = ['kjmol/rad**2', 'deg']
            self.add_term(Harmonic, [BendAngle(*angle)], types, ['PT_ALL', 'HC_FC_DIAG'], units)

    def init_dihedral_terms(self, verbose=False):
        '''
            Estimate the dihedral potentials from the local topology. The
            dihedral potential will be one of the two following possibilities:
            
                If a multiplicity m can be found such that the equilibrium value
                of all instances of the torsion are within `thresshold` of 0 deg
                or per/2 with per = 180deg/m, the following potential will be
                chosen:
                
                    0.5*K*(1-cos(m*psi-m*psi0)) with psi0 = 0 or 360/(2*m) 
                
                If no such multiplicity can be found, but one can found a rest
                value psi0 such that the equilibrium values of all instances of
                the torsion are within `thresshold` of psi0, psi0-180deg, -psi0
                and 180deg-psi0, the following potential will be chosen:
                
                    0.5*K*(cos(2*psi)-cos(2*psi0))**2
                    
                    which is equal to a Yaff PolyFour term
                    
                    a0*cos(psi) + a1*cos(psi)^2 + a2*cos(psi)^3 + a3*cos(psi)^4
                    
                    with a0=0, a1=K*-4*cos(psi0)**2, a2=0, a3=K*2
        '''
        #get all dihedrals
        from molmod.ic import dihed_angle
        ffatypes = [self.system.ffatypes[fid] for fid in self.system.ffatype_ids]
        dihedrals = {}
        for dihedral in self.system.iter_dihedrals():
            dihedral, types = term_sort_atypes(ffatypes, dihedral, 'dihedral')
            if types in dihedrals.keys():
                dihedrals[types].append(dihedral)
            else:
                dihedrals[types] = [dihedral]
        #loop over all distinct dihedral types
        for types, diheds in dihedrals.iteritems():
            psi0s = np.zeros(len(diheds), float)
            ms = np.zeros(len(diheds), float)
            for i, dihed in enumerate(diheds):
                rs = np.array([self.system.pos[j] for j in dihed])
                psi0s[i] = dihed_angle(rs)[0]
                n1 = len(self.system.neighs1[dihed[1]])
                n2 = len(self.system.neighs1[dihed[2]])
                ms[i] = get_multiplicity(n1, n2)
            if None in ms or ms.std()>1e-3:
                print 'WARNING: no dihedral potential for %s (no multiplicity found)' %('.'.join(types))
                continue
            m = int(np.round(ms.mean()))
            rv = get_restvalue(psi0s, m, thresshold=5*deg)
            if rv is not None:
                #a regular Cosine term is used for the dihedral potential
                for dihed in diheds:
                    term = self.add_term(Cosine, [DihedAngle(*dihed)], types, ['HC_FC_DIAG'], ['au', 'kjmol', 'deg'])
                    self.set_params(term.index, rv0=rv, m=m)
            elif m==2:
                #try to determine a term harmonic in the cos(2*psi)
                #fold all dihedral angles to the range [0,90]
                folded = []
                for psi0 in psi0s:
                    if abs(psi0)<90*deg:
                        folded.append(abs(psi0))
                    else:
                        folded.append(180*deg-abs(psi0))
                if np.array(folded).std()<5*deg:
                    #a PolyFour(Cos(2*psi)) is added
                    for dihed in diheds:
                        term = self.add_term(PolyFour, [DihedCos(*dihed)], types, ['HC_FC_DIAG'], ['kjmol', 'kjmol', 'kjmol', 'kjmol'])
                        a1 = -4.0*np.cos(np.array(folded).mean())**2*1.0*kjmol
                        self.set_params(term.index, a0=0.0, a1=a1, a2=0.0, a3=2.0*kjmol)
                else:
                    #no dihedral potential could be determine, hence it is ignored
                    print 'WARNING: no dihedral potential for %s (attempted TORSC2HARM but folded angles to far appart, %.3e deg)' %('.'.join(types), np.array(folded).std()/deg)
                    continue
            else:
                #no dihedral potential could be determine, hence it is ignored
                print 'WARNING: no dihedral potential for %s (no TORSION or TORSC2HARM attempted)' %('.'.join(types))
                continue

    def init_oop_terms(self):
        #get all dihedrals
        from molmod.ic import opbend_dist
        ffatypes = [self.system.ffatypes[fid] for fid in self.system.ffatype_ids]
        opdists = {}
        for opdist in self.system.iter_oops():
            opdist, types = term_sort_atypes(ffatypes, opdist, 'opdist')
            if types in opdists.keys():
                opdists[types].append(opdist)
            else:
                opdists[types] = [opdist]
        #loop over all distinct opdist types
        for types, oops in opdists.iteritems():
            d0s = np.zeros(len(oops), float)
            for i, oop in enumerate(oops):
                rs = np.array([#mind the order, is different in molmod and yaff
                    self.system.pos[oop[3]],
                    self.system.pos[oop[0]],
                    self.system.pos[oop[1]],
                    self.system.pos[oop[2]],
                ])
                d0s[i] = abs(opbend_dist(rs)[0])
            if d0s.mean()<5e-2:
                #add regular term harmonic in oopdist
                for oop in oops:
                    term = self.add_term(Harmonic, [OopDist(*oop)], types, ['HC_FC_DIAG'], ['kjmol/A**2', 'A'])
                    self.set_params(term.index, rv0=0.0)
            else:
                #add term harmonic in square of oopdist
                for oop in oops:
                    self.add_term(Harmonic, [SqOopDist(*oop)], types, ['PT_ALL', 'HC_FC_DIAG'], ['kjmol/A**4', 'A**2'])
    
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
    
    def calc_energy(self, pos):
        old =  self.system.pos.copy()
        self.system.pos = pos.copy()
        self.dlist.forward()
        self.iclist.forward()
        self.vlist.forward()
        energy = self.compute()
        self.system.pos = old
        return energy
    
    def get_hessian_contrib(self, index, fc=None):
        '''
            Get the contribution to the covalent hessian of term with given
            index (and its slaves). If fc is given, set the fc of the master
            and its slave to the given fc.
        '''
        val = ForcePartValence(self.system)
        kind = self.vlist.vtab[index]['kind']
        masterslaves = [index]+self.terms[index].slaves
        kind_to_term = {0: Harmonic, 1: PolyFour, 2: Fues, 3: Cross, 4: Cosine}
        if kind==4:#Cosine
            m, k, rv = self.get_params(index)
            if fc is not None: k = fc
            for jterm in masterslaves:
                ics = self.terms[jterm].ics
                args = (m, k, rv) + tuple(ics)
                val.add_term(Cosine(*args))
        elif kind==3:#cross
            k, rv0, rv1 = self.get_params(index)
            if fc is not None: k = fc
            for jterm in masterslaves:
                ics = self.terms[jterm].ics
                args = (k, rv0, rv1) + tuple(ics)
                val.add_term(Cross(*args))
        elif kind==1:#Polyfour
            a0, a1, a2, a3 = list(self.get_params(index))
            s = 1.0
            if fc is not None:
                s = 2.0*fc/a3
            for jterm in masterslaves:
                ics = self.terms[jterm].ics
                args = ([s*a0,s*a1,s*a2,s*a3],)+tuple(ics)
                val.add_term(PolyFour(*args))
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
    
    def set_params(self, term_index, fc=None, rv0=None, rv1=None, m=None,
            a0=None, a1=None, a2=None, a3=None):
        term = self.vlist.vtab[term_index]
        if term['kind'] in [0,2]:#['Harmonic', 'Fues']
            if fc is not None:  term['par0'] = fc
            if rv0 is not None: term['par1'] = rv0
        elif term['kind'] in [1]:#['PolyFour']
            if a0 is not None: term['par0'] = a0
            if a1 is not None: term['par1'] = a1
            if a2 is not None: term['par2'] = a2
            if a3 is not None: term['par3'] = a3
            if fc is not None:
                pars = self.get_params(term_index)
                assert pars[3] is not None
                term['par1'] = 2.0*fc/pars[3]*pars[1]
                term['par3'] = 2.0*fc
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
        elif term['kind'] in [1]:#['PolyFour']
            if only.lower()=='all': return term['par0'], term['par1'], term['par2'], term['par3']
            elif only.lower()=='a0': return term['par0']
            elif only.lower()=='a1': return term['par1']
            elif only.lower()=='a2': return term['par2']
            elif only.lower()=='a3': return term['par3']
            elif only.lower()=='fc': return 0.5*term['par3']
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
        sequence = ['bondharm', 'bendaharm', 'torsion', 'torsc2harm', 'oopdist', 'cross']
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
            if K<1e-6*kjmol/angstrom**2: continue
            pars.lines.append('%8s  %8s  %.10e  %.10e' %(
                ffatypes[0], ffatypes[1], K/(kjmol/angstrom**2), q0/angstrom
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
            if K<1e-6*kjmol: continue
            pars.lines.append('%8s  %8s  %8s  %.10e  %.10e' %(
                ffatypes[0], ffatypes[1], ffatypes[2], K/kjmol, q0/deg
            ))
        return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})

    def _torsions_to_yaff(self):
        'construct a dihedral section of a yaff parameter file'
        prefix = 'TORSION'
        units = ParameterDefinition('UNIT', lines=['A kjmol', 'PHI0 deg'])
        pars = ParameterDefinition('PARS')
        for i, master in enumerate(self.iter_masters(label=prefix)):
            ffatypes = master.basename.split('/')[1].split('.')
            m, K, q0 = self.get_params(master.index)
            if K<1e-6*kjmol: continue
            pars.lines.append('%8s  %8s  %8s  %8s  %1i %.10e  %.10e' %(
                ffatypes[0], ffatypes[1],  ffatypes[2], ffatypes[3], m,
                K/kjmol, q0/deg
            ))
        return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})
    
    def _torsc2harm_to_yaff(self):
        prefix = 'TORSC2HARM'
        units = ParameterDefinition('UNIT', lines=['A kjmol', 'COS0 au'])
        pars = ParameterDefinition('PARS')
        for i, master in enumerate(self.iter_masters(label=prefix)):
            ffatypes = master.basename.split('/')[1].split('.')
            a = self.get_params(master.index)
            K = 0.5*a[3]
            cos0 = np.arccos(np.sqrt(-0.5*a[1]/a[3]))
            if K<1e-6*kjmol: continue
            pars.lines.append('%8s  %8s  %8s  %8s  %.10e  %.10e' %(
                ffatypes[0], ffatypes[1],  ffatypes[2], ffatypes[3],
                K/kjmol, cos0
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
                K/(kjmol/angstrom**2), q0/angstrom
            ))
        return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})
    
    def _sqopdists_to_yaff(self):
        'construct a opdist section of a yaff parameter file'
        prefix = 'SQOOPDIST'
        units = ParameterDefinition('UNIT', lines=['K kjmol/A**4', 'D0 A**2'])
        pars = ParameterDefinition('PARS')
        for i, master in enumerate(self.iter_masters(label=prefix)):
            ffatypes = master.basename.split('/')[1].split('.')
            K, q0 = self.get_params(master.index)
            if K<1e-6: continue
            pars.lines.append('%8s  %8s  %8s  %8s  %.10e  %.10e' %(
                ffatypes[0], ffatypes[1], ffatypes[2], ffatypes[3],
                K/(kjmol/angstrom**4), q0/angstrom**2
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
                elif 'b0a' in other.basename:
                    b0a = self.get_params(other.index)
                elif 'b1a' in other.basename:
                    b1a = self.get_params(other.index)
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
                '%8s  %8s  %8s  % .10e  % .10e  % .10e  %.10e  %.10e  %.10e' %(
                    ffatypes[0], ffatypes[1], ffatypes[2],
                    Kss/(kjmol/angstrom**2), Kbs0/(kjmol/angstrom),
                    Kbs1/(kjmol/angstrom), r0/angstrom, r1/angstrom, theta0/deg,
            ))
            done.append(label)
        return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})

    def dump_yaff(self, fn):
        sections = [
            self._bonds_to_yaff(), self._bends_to_yaff(), 
            self._torsions_to_yaff(), self._torsc2harm_to_yaff(),
            self._opdists_to_yaff(), self._cross_to_yaff(),
        ]
        f = open(fn, 'w')
        for section in sections:
            if len(section['PARS'].lines)==0: continue
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
