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

from quickff.tools import term_sort_atypes, get_multiplicity, get_restvalue, \
    digits
from quickff.log import log

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

    def get_atoms(self):
        'Get the ordered list of indexes of the atoms involved'
        atoms = None
        if self.kind==3:#cross
            if self.ics[0].kind==0 and self.ics[1].kind in [1,2]:
                ic = self.ics[1]
            elif self.ics[0].kin in [1,2] and self.ics[1].kind==0:
                ic = self.ics[0]
            elif self.ics[0].kind==0 and self.ics[1].kind==0:
                a0 = self.ics[0].index_pairs[0]
                a1 = self.ics[1].index_pairs[0]
                if   a0[1]==a1[0]: return [a0[0], a0[1], a1[1]]
                elif a0[0]==a1[1]: return [a1[0], a1[1], a0[1]]
        else:
            ic = self.ics[0]
        if ic.kind==0:#Bond
            atoms = ic.index_pairs[0]
        elif ic.kind in [1,2]: #bend
            a0 = ic.index_pairs[0]
            a1 = ic.index_pairs[1]
            atoms = [a0[1], a0[0], a1[1]]
        elif ic.kind in [3,4]: #dihedral
            a0 = ic.index_pairs[0]
            a1 = ic.index_pairs[1]
            a2 = ic.index_pairs[2]
            atoms = [a0[1], a0[0], a1[1], a2[1]]
        elif ic.kind in [10,11]: #oopdist
            a0 = ic.index_pairs[0]
            a1 = ic.index_pairs[1]
            a2 = ic.index_pairs[2]
            atoms = [a0[0], a0[1], a1[1], a2[1]]
        if atoms is None:
            raise ValueError('get_atoms not supported for term %s' %self.basename)
        else:
            return atoms

    def to_string(self, valence, max_name=38, max_line=72):
        #check if current is master
        assert self.master is None or self.is_master(), \
            'Current term is not the master'
        #collect parameters
        npars = len(valence.get_params(self.index))
        pars = np.zeros([len(self.slaves)+1, npars], float)
        pars[0,:] = np.array(valence.get_params(self.index))
        for i, index in enumerate(self.slaves):
            pars[1+i,:] = np.array(valence.get_params(index))
        #set default config (applicable for harmonic terms)
        means = pars.mean(axis=0)
        stds = pars.std(axis=0)
        formats = [
            'fc = %%5s %s %%4s' %(u"\u00B1"),
            'rv = %%5s %s %%4s' %(u"\u00B1"),
        ]
        ndigits = [(5,4), (5,4)]
        units = self.units
        #set special config
        if self.kind==1 and self.ics[0].kind==3:#PolyFour for torsc2harm
            fcs = 0.5*pars[:,3].mean()
            rvs = pars[:,0].mean()
            means = fcs.means(), rvs.mean()
            stds = fcs.std(), rvs.std()
            units = [self.units[3], 'deg']
        elif self.kind==3:#cross
            formats = [
                'fc = %%4s %s %%2s' %(u"\u00B1"),
                'rv0 = %%4s %s %%3s' %(u"\u00B1"),
                'rv0 = %%4s %s %%3s' %(u"\u00B1")
            ]
            ndigits = [(4,2), (4,3), (4,3)]
        elif self.kind==4:#cosine
            m, fc, rv = pars.mean(axis=0)
            dm, dfc, drv = pars.std(axis=0)
            means = fc, rv, m
            stds = dfc, drv, np.nan
            formats = [
                'fc = %%4s %s %%3s' %(u"\u00B1"),
                'rv = %%4s %s %%3s' %(u"\u00B1"),
                'm = %1s%0s'
            ]
            units = [self.units[1], self.units[2], 'au']
            ndigits = [(4,3), (4,3), (1,0)]
        #convert term pars to string
        line = '%s (%s)' %(
            self.basename[:max_line],
            '  '.join([unit.replace('**','^') for unit in self.units])
        )
        line += ' '*(max_line-len(line))
        for fmt, mean, std, ndigit, unit in zip(formats, means, stds, ndigits, units):
            smean = digits(mean/parse_unit(unit), ndigit[0])
            sstd = digits(std/parse_unit(unit), ndigit[1])
            line += '    ' + fmt %(smean, sstd)
        return line


class ValenceFF(ForcePartValence):
    '''
        Class to collect all valence terms in the force field for which
        parameters need to be estimated.
    '''
    def __init__(self, system, specs=None):
        '''
            **Arguments**

            system
                an instance of the Yaff System class containing all system
                properties

            **Keyword Arguments**

            specs
                Not yet implemented
        '''
        with log.section('VAL', 2, timer='Initializing'):
            log.dump('Initializing valence force field')
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

            **Arguments**

            pot
                an instance of ValenceTerm from `yaff.pes.vlist.py` representing
                the potential chosen for the term

            ics
                list of InternalCoordinate instances from `yaff.pes.iclist.py`

            atypes
                ordered string of atom types (required to assign name to the
                term)

            tasks
                List of strings defining all tasks that have to be performed for
                this term. Possible tasks are: PT_ALL, HC_FC_DIAG, HC_FC_CROSS,
                EQ_RV.

            units
                Units for all parameters in this term (ordered in the same
                way as parameters are stored in `yaff.vlist.vtab`)
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

    def modify_term(self, term_index, pot, ics, basename, tasks, units):
        '''
            Modify the term with given index to a new valence term.
        '''
        #TODO: only allow masters and automatically also modify the slaves
        with log.section('VAL', 2):
            #modify in valence.terms
            old_term = self.terms[term_index]
            new_term = Term(
                term_index, basename, pot.kind, ics, tasks,
                units, master=old_term.master, slaves=old_term.slaves
            )
            self.terms[term_index] = new_term
            #modify in valence.vlist.vtab
            vterm = self.vlist.vtab[term_index]
            if pot.kind==1:#all 4 parameters of PolyFour are given as 1 tuple
                args = [(None,)*len(units)] + ics
            else:
                args = [None,]*len(units) + ics
            new = pot(*args)
            vterm['kind'] = new.kind
            for i in xrange(len(new.pars)):
                vterm['par%i'%i] = new.pars[i]
            ic_indexes = new.get_ic_indexes(self.iclist)
            for i in xrange(len(ic_indexes)):
                vterm['ic%i'%i] = ic_indexes[i]

    def iter_terms(self, label=None):
        '''
            Iterate over all terms in the valence force field. If label is
            given, only iterate over terms with the label its name.
        '''
        for term in self.terms:
            if label is None:
                yield term
            elif label.startswith('/'):
                trimmed = label[1:]
                if trimmed.lower() in term.basename.lower()[:len(trimmed)]:
                    yield term
            elif label.lower() in term.basename.lower():
                yield term

    def iter_masters(self, label=None):
        '''
            Iterate over all master terms in the valence force field. If label
            is given, only iterate of the terms with the label in its name.
        '''
        for term in self.iter_terms(label=label):
            if term.is_master():
                yield term

    def init_bond_terms(self):
        '''
            Initialize all bond terms in the system based on the bonds attribute
            of the system instance. All bond terms are given harmonic
            potentials.
        '''
        with log.section('VAL', 3, 'Initializing'):
            ffatypes = [self.system.ffatypes[fid] for fid in self.system.ffatype_ids]
            #get the bond terms
            nbonds = 0
            for bond in self.system.iter_bonds():
                bond, types = term_sort_atypes(ffatypes, bond, 'bond')
                units = ['kjmol/A**2', 'A']
                self.add_term(Harmonic, [Bond(*bond)], types, ['PT_ALL', 'HC_FC_DIAG'], units)
                nbonds += 1
        log.dump('Added %i Harmonic bond terms' %nbonds)

    def init_bend_terms(self):
        '''
            Initialize all bend terms in the system based on the bends attribute
            of the system instance. All bend terms are given harmonic
            potentials.
        '''
        with log.section('VAL', 3, 'Initializing'):
            ffatypes = [self.system.ffatypes[fid] for fid in self.system.ffatype_ids]
            #get the angle terms
            nbends = 0
            for angle in self.system.iter_angles():
                angle, types = term_sort_atypes(ffatypes, angle, 'angle')
                units = ['kjmol/rad**2', 'deg']
                self.add_term(Harmonic, [BendAngle(*angle)], types, ['PT_ALL', 'HC_FC_DIAG'], units)
                nbends += 1
        log.dump('Added %i Harmonic bend terms' %nbends)

    def init_dihedral_terms(self, thresshold=20*deg):
        '''
            Initialize the dihedral potentials from the local topology. The
            dihedral potential will be one of the two following possibilities:

                The multiplicity m is determined from the local topology, i.e.
                the number of neighbors of the central two atoms in the dihedral

                If the equilibrium value of all instances of the torsion are
                within `thresshold` of 0 deg or per/2 with per = 180deg/m,
                the following potential will be chosen:

                    0.5*K*(1-cos(m*psi-m*psi0)) with psi0 = 0 or 360/(2*m)
        '''
        with log.section('VAL', 3, 'Initializing'):
            #get all dihedrals
            from molmod.ic import dihed_angle, _dihed_angle_low
            ffatypes = [self.system.ffatypes[fid] for fid in self.system.ffatype_ids]
            dihedrals = {}
            for dihedral in self.system.iter_dihedrals():
                dihedral, types = term_sort_atypes(ffatypes, dihedral, 'dihedral')
                if types in dihedrals.keys():
                    dihedrals[types].append(dihedral)
                else:
                    dihedrals[types] = [dihedral]
            #loop over all distinct dihedral types
            ncos = 0
            for types, diheds in dihedrals.iteritems():
                psi0s = np.zeros(len(diheds), float)
                ms = np.zeros(len(diheds), float)
                for i, dihed in enumerate(diheds):
                    if self.system.cell.nvec>0:
                        d10 = self.system.pos[dihed[0]] - self.system.pos[dihed[1]]
                        d12 = self.system.pos[dihed[2]] - self.system.pos[dihed[1]]
                        d23 = self.system.pos[dihed[3]] - self.system.pos[dihed[2]]
                        self.system.cell.mic(d10)
                        self.system.cell.mic(d12)
                        self.system.cell.mic(d23)
                        psi0s[i] = _dihed_angle_low(d10, d12, d23, 0)[0]
                    else:
                        rs = np.array([self.system.pos[j] for j in dihed])
                        psi0s[i] = dihed_angle(rs)[0]
                    n1 = len(self.system.neighs1[dihed[1]])
                    n2 = len(self.system.neighs1[dihed[2]])
                    ms[i] = get_multiplicity(n1, n2)
                if np.isnan(ms).any() or ms.std()>1e-3:
                    ms_string = str(ms)
                    if np.isnan(ms).any():
                        ms_string = 'nan'
                    log.warning('missing dihedral for %s (m is %s)' %('.'.join(types), ms_string))
                    continue
                m = int(np.round(ms.mean()))
                rv = get_restvalue(psi0s, m, thresshold=thresshold, mode=1)
                if rv is not None:
                    #a regular Cosine term is used for the dihedral potential
                    for dihed in diheds:
                        term = self.add_term(Cosine, [DihedAngle(*dihed)], types, ['HC_FC_DIAG'], ['au', 'kjmol', 'deg'])
                        self.set_params(term.index, rv0=rv, m=m)
                        ncos += 1
                else:
                    #no dihedral potential could be determine, hence it is ignored
                    log.warning('missing dihedral for %s (could not determine rest value from %s)' %('.'.join(types), str(psi0s/deg)))
                    continue
        log.dump('Added %i Cosine dihedral terms' %ncos)

    def init_oop_terms(self, thresshold_zero=5e-2*angstrom):
        '''
            Initialize all out-of-plane terms in the system based on the oops
            attribute of the system instance. All oops are given harmonic
            potentials.
        '''
        with log.section('VAL', 3, 'Initializing'):
            #get all dihedrals
            from molmod.ic import opbend_dist, _opdist_low
            ffatypes = [self.system.ffatypes[fid] for fid in self.system.ffatype_ids]
            opdists = {}
            for opdist in self.system.iter_oops():
                opdist, types = term_sort_atypes(ffatypes, opdist, 'opdist')
                if types in opdists.keys():
                    opdists[types].append(opdist)
                else:
                    opdists[types] = [opdist]
            #loop over all distinct opdist types
            nharm = 0
            nsq = 0
            for types, oops in opdists.iteritems():
                d0s = np.zeros(len(oops), float)
                for i, oop in enumerate(oops):
                    if self.system.cell.nvec>0:
                        d01 = self.system.pos[oop[1]]-self.system.pos[oop[0]]
                        d02 = self.system.pos[oop[2]]-self.system.pos[oop[0]]
                        d03 = self.system.pos[oop[3]]-self.system.pos[oop[0]]
                        self.system.cell.mic(d01)
                        self.system.cell.mic(d02)
                        self.system.cell.mic(d03)
                        d0s[i] = abs(_opdist_low(d01, d02, d03, 0)[0])
                    else:
                        rs = np.array([#mind the order, is(or was) wrongly documented in molmod
                            self.system.pos[oop[0]],
                            self.system.pos[oop[1]],
                            self.system.pos[oop[2]],
                            self.system.pos[oop[3]],
                        ])
                        d0s[i] = abs(opbend_dist(rs)[0])
                if d0s.mean()<thresshold_zero: #TODO: check this thresshold
                    #add regular term harmonic in oopdist
                    for oop in oops:
                        term = self.add_term(Harmonic, [OopDist(*oop)], types, ['HC_FC_DIAG'], ['kjmol/A**2', 'A'])
                        self.set_params(term.index, rv0=0.0)
                        nharm += 1
                else:
                    #add term harmonic in square of oopdist
                    log.dump('Mean absolute value of OopDist %s is %.3e A, used SQOOPDIST' %('.'.join(types), d0s.mean()/angstrom))
                    for oop in oops:
                        self.add_term(Harmonic, [SqOopDist(*oop)], types, ['PT_ALL', 'HC_FC_DIAG'], ['kjmol/A**4', 'A**2'])
                        nsq += 1
        log.dump('Added %i Harmonic and %i SquareHarmonic out-of-plane distance terms' %(nharm, nsq))

    def init_cross_terms(self, specs=None):
        '''
            Initialize cross terms between bonds and bends.
        '''
        with log.section('VAL', 3, 'Initializing'):
            ffatypes = [self.system.ffatypes[i] for i in self.system.ffatype_ids]
            nbend = 0
            for angle in self.system.iter_angles():
                angle, types = term_sort_atypes(ffatypes, angle, 'angle')
                skip_angle = False
                for term in self.iter_masters('.'.join(types)):
                    if len(term.ics)>1: continue
                    ickind = term.ics[0].kind
                    potkind = term.kind
                    if ickind==1 or (ickind==2 and potkind!=0):
                        skip_angle = True
                        log.dump('Skipped stretch-angle cross term for %s due to incompatible diagonal bend term with potkind=%i and ickind=%i' %('.'.join(types), potkind, ickind))
                        break
                bond0, btypes = term_sort_atypes(ffatypes, angle[:2], 'bond')
                bond1, btypes = term_sort_atypes(ffatypes, angle[1:], 'bond')
                #add stretch-stretch
                self.add_term(Cross,
                    [Bond(*bond0), Bond(*bond1)],
                    types, ['HC_FC_CROSS'], ['kjmol/A**2', 'A', 'A']
                )
                if not skip_angle:
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
                nbend += 1
        log.dump('Added %i cross terms for angle patterns' %nbend)

    def apply_constraints(self, constraints):
        '''
            Routine to apply equality constraints in the force field fitting

            **Arguments**

            contraints
                A dictionairy containing (master: slaves) definitions in which
                master is a string defining the master basename and slaves is a
                list of strings defining the slave basenames.
        '''
        for mastername, slavenames in constraints.iteritems():
            masters = [term for term in self.iter_masters(mastername)]
            assert len(masters)==1, 'Master %s is not uniquely defined' %mastername
            master = masters[0]
            for slavename in slavenames:
                for slave in self.iter_terms(slavename):
                    slave.basename = master.basename
                    slave.master = master.index
                    slave.slaves = []
                    master.slaves.append(slave.index)

    def calc_energy(self, pos):
        old =  self.system.pos.copy()
        self.system.pos = pos.copy()
        self.dlist.forward()
        self.iclist.forward()
        energy = self.vlist.forward()
        #energy = self.compute()
        self.system.pos = old
        self.dlist.forward()
        self.iclist.forward()
        self.vlist.forward()
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
            if fc is not None:
                a3 = 2.0*fc
                a1 = -4.0*fc*np.cos(a0)**2
            for jterm in masterslaves:
                ics = self.terms[jterm].ics
                args = ([0.0,a1,0.0,a3],)+tuple(ics)
                val.add_term(PolyFour(*args))
        elif kind==0:#Harmonic:
            k, rv = self.get_params(index)
            if fc is not None: k = fc
            for jterm in masterslaves:
                ics = self.terms[jterm].ics
                args = (k, rv) + tuple(ics)
                val.add_term(Harmonic(*args))
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
            if fc is not None or rv0 is not None:
                if fc is None:  fc = self.get_params(term_index, only='fc')
                if rv0 is None: rv0 = self.get_params(term_index, only='rv')
                term['par0'] = rv0
                term['par1'] = -4.0*fc*np.cos(rv0)**2
                term['par2'] = 0.0
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
            elif only.lower()=='rv': return term['par0']
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
            raise NotImplementedError(
                'get_params not implemented for Yaff %s term' % term['kind'])

    def is_negligible(self, term_index):
        """Return True if the given term can be neglected (e.g. in parameter files)."""
        # Note that the units may not be strictly correct: below per angstrom, radian
        # squared, etc should be used. Because that would not affect the order of
        # magnitude of the threshold, it is not worth adding such complications to the
        # code.
        term = self.vlist.vtab[term_index]
        if term['kind'] in [0, 2]:  # ['Harmonic', 'Fues']
            return abs(term['par0']) < 1e-6*kjmol
        elif term['kind'] in [1]:   # ['PolyFour']
            # Not sure how to handle this one...
            # For now, never neglect.
            return False
        elif term['kind'] in [4]:   # ['Cosine']
            return abs(term['par1']) < 1e-6*kjmol
        elif term['kind'] in [3]:   # ['Cross']
            return abs(term['par0']) < 1e-6*kjmol
        else:
            raise NotImplementedError(
                'is_negligible not implemented for Yaff %s term' % term['kind'])

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
            if label=='all':
                for i, par in enumerate(value):
                    assert not np.isnan(par), 'Par%i of %s is not set' %(i, term.basename)
            else:
                assert not np.isnan(value), '%s of %s is not set' %(label, term.basename)

    def dump_logger(self, print_level=1):
        if log.log_level<print_level: return
        with log.section('', print_level):
            sequence = [
                'bondharm',
                'bendaharm', 'bendcharm', 'bendcos',
                'torsion', 'torsc2harm', 'dihedharm',
                'oopdist', 'cross'
            ]
            log.dump('')
            for label in sequence:
                lines = []
                for term in self.iter_masters(label=label):
                    lines.append(term.to_string(self))
                for line in sorted(lines):
                    log.dump(line)
                    log.dump('')

    def _bonds_to_yaff(self):
        'construct a BONDHARM section of a yaff parameter file'
        prefix = 'BONDHARM'
        units = ParameterDefinition('UNIT', lines=['K kjmol/A**2', 'R0 A'])
        pars = ParameterDefinition('PARS')
        for i, master in enumerate(self.iter_masters(label=prefix)):
            if self.is_negligible(i): continue
            ffatypes = master.basename.split('/')[1].split('.')
            K, q0 = self.get_params(master.index)
            pars.lines.append('%8s  %8s  %.10e  %.10e' %(
                ffatypes[0], ffatypes[1], K/(kjmol/angstrom**2), q0/angstrom
            ))
        return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})

    def _bendaharm_to_yaff(self):
        'construct a BENDAHARM section of a yaff parameter file'
        prefix = 'BENDAHARM'
        units = ParameterDefinition('UNIT', lines=['K kjmol/rad**2', 'THETA0 deg'])
        pars = ParameterDefinition('PARS')
        for i, master in enumerate(self.iter_masters(label=prefix)):
            if self.is_negligible(i): continue
            ffatypes = master.basename.split('/')[1].split('.')
            K, q0 = self.get_params(master.index)
            pars.lines.append('%8s  %8s  %8s  %.10e  %.10e' %(
                ffatypes[0], ffatypes[1], ffatypes[2], K/kjmol, q0/deg
            ))
        return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})

    def _bendcharm_to_yaff(self):
        'construct a BENDCHARM section of a yaff parameter file'
        prefix = 'BENDCHARM'
        units = ParameterDefinition('UNIT', lines=['K kjmol', 'COS0 au'])
        pars = ParameterDefinition('PARS')
        for i, master in enumerate(self.iter_masters(label=prefix)):
            if self.is_negligible(i): continue
            ffatypes = master.basename.split('/')[1].split('.')
            K, q0 = self.get_params(master.index)
            pars.lines.append('%8s  %8s  %8s  %.10e  %.10e' %(
                ffatypes[0], ffatypes[1], ffatypes[2], K/kjmol, q0
            ))
        return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})

    def _bendcos_to_yaff(self):
        'construct a BENDCOS section of a yaff parameter file'
        prefix = 'BENDCOS'
        units = ParameterDefinition('UNIT', lines=['A kjmol', 'PHI0 deg'])
        pars = ParameterDefinition('PARS')
        for i, master in enumerate(self.iter_masters(label=prefix)):
            if self.is_negligible(i): continue
            ffatypes = master.basename.split('/')[1].split('.')
            m, K, q0 = self.get_params(master.index)
            pars.lines.append('%8s  %8s  %8s  %1i %.10e  %.10e' %(
                ffatypes[0], ffatypes[1], ffatypes[2], m, K/kjmol, q0/deg
            ))
        return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})

    def _torsions_to_yaff(self):
        'construct a TORSION section of a yaff parameter file'
        prefix = 'TORSION'
        units = ParameterDefinition('UNIT', lines=['A kjmol', 'PHI0 deg'])
        pars = ParameterDefinition('PARS')
        for i, master in enumerate(self.iter_masters(label=prefix)):
            if self.is_negligible(i): continue
            ffatypes = master.basename.split('/')[1].split('.')
            m, K, q0 = self.get_params(master.index)
            pars.lines.append('%8s  %8s  %8s  %8s  %1i %.10e  %.10e' %(
                ffatypes[0], ffatypes[1],  ffatypes[2], ffatypes[3], m,
                K/kjmol, q0/deg
            ))
        return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})

    def _torsc2harm_to_yaff(self):
        'construct a TORSC2HARM section of a yaff parameter file'
        prefix = 'TORSC2HARM'
        units = ParameterDefinition('UNIT', lines=['A kjmol', 'COS0 au'])
        pars = ParameterDefinition('PARS')
        for i, master in enumerate(self.iter_masters(label=prefix)):
            if self.is_negligible(i): continue
            ffatypes = master.basename.split('/')[1].split('.')
            a = self.get_params(master.index)
            K = 0.5*a[3]
            cos0 = np.arccos(np.sqrt(-0.5*a[1]/a[3]))
            pars.lines.append('%8s  %8s  %8s  %8s  %.10e  %.10e' %(
                ffatypes[0], ffatypes[1],  ffatypes[2], ffatypes[3],
                K/kjmol, cos0
            ))
        return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})

    def _dihedharm_to_yaff(self):
        'construct a DIHEDHARM section of a yaff parameter file'
        prefix = 'DIHEDHARM'
        units = ParameterDefinition('UNIT', lines=['K kjmol/rad**2', 'PHI0 deg'])
        pars = ParameterDefinition('PARS')
        for i, master in enumerate(self.iter_masters(label=prefix)):
            if self.is_negligible(i): continue
            ffatypes = master.basename.split('/')[1].split('.')
            K, phi0 = self.get_params(master.index)
            pars.lines.append('%8s  %8s  %8s  %8s  %.10e  %.10e' %(
                ffatypes[0], ffatypes[1],  ffatypes[2], ffatypes[3],
                K/kjmol, phi0/deg
            ))
        return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})

    def _opdists_to_yaff(self):
        'construct a OOPDIST section of a yaff parameter file'
        prefix = 'OOPDIST'
        units = ParameterDefinition('UNIT', lines=['K kjmol/A**2', 'D0 A'])
        pars = ParameterDefinition('PARS')
        for i, master in enumerate(self.iter_masters(label=prefix)):
            if self.is_negligible(i): continue
            if 'sqoopdist' in master.basename.lower(): continue
            ffatypes = master.basename.split('/')[1].split('.')
            K, q0 = self.get_params(master.index)
            pars.lines.append('%8s  %8s  %8s  %8s  %.10e  %.10e' %(
                ffatypes[0], ffatypes[1], ffatypes[2], ffatypes[3],
                K/(kjmol/angstrom**2), q0/angstrom
            ))
        return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})

    def _sqopdists_to_yaff(self):
        'construct a SQOOPDIST section of a yaff parameter file'
        prefix = 'SQOOPDIST'
        units = ParameterDefinition('UNIT', lines=['K kjmol/A**4', 'D0 A**2'])
        pars = ParameterDefinition('PARS')
        for i, master in enumerate(self.iter_masters(label=prefix)):
            if self.is_negligible(i): continue
            ffatypes = master.basename.split('/')[1].split('.')
            K, q0 = self.get_params(master.index)
            pars.lines.append('%8s  %8s  %8s  %8s  %.10e  %.10e' %(
                ffatypes[0], ffatypes[1], ffatypes[2], ffatypes[3],
                K/(kjmol/angstrom**4), q0/angstrom**2
            ))
        return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})

    def _cross_to_yaff(self):
        'construct a CROSS section of a yaff parameter file'
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
            if self.is_negligible(i): continue
            prefix, ffatypes, suffix = master.basename.split('/')
            label = prefix+'/'+ffatypes+'/'
            if label in done: continue
            bb, b0a, b1a = None, None, None
            for j, other in enumerate(self.iter_masters(label=label)):
                if 'bb' in other.basename:
                    bb = self.get_params(other.index)
                elif 'b0a' in other.basename:
                    b0a = self.get_params(other.index)
                elif 'b1a' in other.basename:
                    b1a = self.get_params(other.index)
                else:
                    raise ValueError('Invalid Cross term %s' %other.basename)
            Kss, Kbs0, Kbs1 = 0.0, 0.0, 0.0
            r0, r1, theta0 = 0.0, 0.0, 0.0
            if bb is not None: Kss, r0, r1 = bb
            if b0a is not None: Kbs0, r0, theta0 = b0a
            if b1a is not None: Kbs1, r1, theta0 = b1a
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
            self._bonds_to_yaff(),
            self._bendaharm_to_yaff(), self._bendcharm_to_yaff(),
            self._bendcos_to_yaff(),
            self._torsions_to_yaff(), self._torsc2harm_to_yaff(),
            self._dihedharm_to_yaff(),
            self._opdists_to_yaff(),self._sqopdists_to_yaff(),
            self._cross_to_yaff(),
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
