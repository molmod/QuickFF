# -*- coding: utf-8 -*-
# QuickFF is a code to quickly derive accurate force fields from ab initio input.
# Copyright (C) 2012 - 2019 Louis Vanduyfhuys <Louis.Vanduyfhuys@UGent.be>
# Steven Vandenbrande <Steven.Vandenbrande@UGent.be>,
# Jelle Wieme <Jelle.Wieme@UGent.be>,
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

from __future__ import unicode_literals, absolute_import

from molmod.units import *
from molmod.ic import bend_angle, _bend_angle_low, dihed_angle, _dihed_angle_low

from yaff.pes.ff import ForceField, ForcePartValence
from yaff.pes.vlist import *
from yaff.pes.iclist import InternalCoordinateList
from yaff.pes.iclist import Bond, BendAngle, BendCos, DihedCos, DihedAngle, \
    OopDist, SqOopDist
from yaff.pes.dlist import DeltaList
from yaff.sampling.harmonic import estimate_cart_hessian

from quickff.tools import term_sort_atypes, get_multiplicity, get_restvalue, \
    digits
from quickff.log import log

import numpy as np, re

__all__ = ['ValenceFF']

class Term(object):
    '''
        A class to store easy-accessible information about a term included in
        the valence force field
    '''
    def __init__(self, index, basename, kind, ics, tasks, units,master=None, slaves=None, diag_term_indexes=[]):
        self.index = index
        self.basename = basename
        self.kind = kind
        self.ics = ics
        self.tasks = tasks
        self.units = units
        self.master = master
        self.slaves = slaves
        self.diag_term_indexes=diag_term_indexes

    def is_master(self):
        return self.master==self.index

    def get_atoms(self):
        'Get the ordered list of indexes of the atoms involved'
        atoms = None
        ic = None
        if self.kind==3:#cross
            #check if one of ics is dihedral
            for testic in self.ics:
                if testic.kind in [3,4,12,13,14,15]:
                    ic = testic
                    break
            #check if one of ics is angle
            if ic is None:
                for testic in self.ics:
                    if testic.kind in [1,2]:
                        ic = testic
                        break
            #only bonds in ics, either r01,r12 or r01,r23
            if ic is None:
                assert self.ics[0].kind==0 and self.ics[1].kind==0, 'IC other than bond, bendangle, bendcos, dihedangle or dihedcos in cross term. Not implemented!'
                a0 = self.ics[0].index_pairs[0]
                a1 = self.ics[1].index_pairs[0]
                if   a0[1]==a1[0]:
                    return [a0[0], a0[1], a1[1]]
                elif a0[0]==a1[1]:
                    return [a1[0], a1[1], a0[1]]
                else:
                    return [a0[0], a0[1], a1[0], a1[1]]
        else:
            ic = self.ics[0]
        assert ic is not None
        if ic.kind==0:#Bond
            atoms = ic.index_pairs[0]
        elif ic.kind in [1,2]: #bend
            a0 = ic.index_pairs[0]
            a1 = ic.index_pairs[1]
            atoms = [a0[1], a0[0], a1[1]]
        elif ic.kind in [3,4,12,13,14,15]: #dihedral
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
            'fc = %%5s %s %%4s' %("+-"),
            'rv = %%5s %s %%4s' %("+-"),
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
                'fc = %%4s %s %%2s' %("+-"),
                'rv0 = %%4s %s %%3s' %("+-"),
                'rv1 = %%4s %s %%3s' %("+-")
            ]
            ndigits = [(4,2), (4,3), (4,3)]
        elif self.kind==4:#cosine
            m, fc, rv = pars.mean(axis=0)
            dm, dfc, drv = pars.std(axis=0)
            means = fc, rv, m
            stds = dfc, drv, np.nan
            formats = [
                'fc = %%4s %s %%3s' %("+-"),
                'rv = %%4s %s %%3s' %("+-"),
                'm = %1s%0s'
            ]
            units = [self.units[1], self.units[2], 'au']
            ndigits = [(4,3), (4,3), (1,0)]
        elif self.kind in [5, 6, 7, 8, 9]:#chebychev
            sign = pars[:,1].mean()
            fcs = pars[:,0]
            means = fcs.mean(), sign
            stds = fcs.std(), np.nan
            formats = [
                'fc  = %%4s %s %%3s' %("+-"),
                'sgn = %3s%0s',
            ]
            units = [self.units[0], 'au']
            ndigits = [(4,3), (3,0)]
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
    def __init__(self, system, settings):
        '''
            **Arguments**

            system
                an instance of the Yaff System class containing all system
                properties

            settings
                an instance of `Settings` containing all QuickFF settings
        '''
        with log.section('VAL', 2, timer='Initializing'):
            log.dump('Initializing valence force field')
            self.system = system
            self.settings = settings
            self.terms = []
            ForcePartValence.__init__(self, system)
            if self.settings.do_bonds:
                self.init_bond_terms()
            if self.settings.do_bends:
                self.init_bend_terms()
            if self.settings.do_dihedrals:
                self.init_dihedral_terms()
            if self.settings.do_oops:
                self.init_oop_terms()

    def add_term(self, pot, ics, basename, tasks, units, diag_term_indexes=[]):
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

            basename
                base name for the term, identical for master and each slave

            tasks
                List of strings defining all tasks that have to be performed for
                this term. Possible tasks are: PT_ALL, HC_FC_DIAG, HC_FC_CROSS,
                EQ_RV.

            units
                Units for all parameters in this term (ordered in the same
                way as parameters are stored in `yaff.vlist.vtab`)

            **Optional arguments**

            diag_term_indexes
                Indexes of the diagonal terms that correspond to the ics
                composing an off-diagonal term. Empty if the term is diagonal.
        '''
        index = len(self.terms)
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
            units, master=master, slaves=slaves, diag_term_indexes=diag_term_indexes
        )
        self.terms.append(term)
        if pot.kind==1:#all 4 parameters of PolyFour are given as 1 tuple
            args = [(None,)*len(units)] + ics
        elif pot.kind in [5,6,7,8,9]: #Chebychev
            args = [None] + ics #sign will first be defaulted as -1, but is set to the correct value in init_dihedral_terms
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
            elif pot.kind in [5,6,7,8,9]: #chebychev potentials
                args = [None,]+ics
            else:
                args = [None,]*len(units) + ics
            new = pot(*args)
            vterm['kind'] = new.kind
            for i in range(len(new.pars)):
                vterm['par%i'%i] = new.pars[i]
            ic_indexes = new.get_ic_indexes(self.iclist)
            for i in range(len(ic_indexes)):
                vterm['ic%i'%i] = ic_indexes[i]

    def iter_terms(self, label=None, use_re=False):
        '''
            Iterate over all terms in the valence force field. If label is
            given, only iterate over terms that matches label.
        '''
        for term in self.terms:
            if label is None:
                yield term
            elif use_re:
                pattern = re.compile(label, re.IGNORECASE)
                if pattern.match(term.basename):
                    yield term
                else:
                    continue
            else:
                if label.lower() in term.basename.lower():
                    yield term
                else:
                    continue

    def iter_masters(self, label=None, use_re=False):
        '''
            Iterate over all master terms in the valence force field. If label
            is given, only iterate of the terms with the label in its name.
        '''
        for term in self.iter_terms(label=label, use_re=use_re):
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

            #list of bonds which should be excluded
            if self.settings.excl_bonds is not None:
                excl_bonds = self.settings.excl_bonds.split(',')

            for bond in self.system.iter_bonds():
                skip = False
                bond, types = term_sort_atypes(ffatypes, bond, 'bond')

                if self.settings.excl_bonds is not None:
                    bond_opt1 = '.'.join(types)
                    bond_opt2 = '.'.join(types[::-1])
                    for excl in excl_bonds:
                        pattern = re.compile(excl, re.IGNORECASE)
                        if pattern.match(bond_opt1) or pattern.match(bond_opt2):
                            skip = True

                if not skip:
                    units = ['kjmol/A**2', 'A']
                    basename = self.settings.bond_term+'/'+'.'.join(types)
                    if self.settings.bond_term.lower()   == 'bondharm':
                        pot = Harmonic
                    elif self.settings.bond_term.lower() == 'bondmm3':
                        pot = MM3Quartic
                    elif self.settings.bond_term.lower() == 'bondfues':
                        pot = Fues
                    else:
                        raise ValueError('Bond kind %s not supported' %self.settings.bond_term)
                    term = self.add_term(pot, [Bond(*bond)], basename, ['PT_ALL', 'HC_FC_DIAG'], units)
                    nbonds += 1
                else:
                    log.dump('Excluded %s bond'%'.'.join(types))
        log.dump('Added %i bond terms' %nbonds)

    def init_bend_terms(self, thresshold=10*deg):
        '''
            Initialize all bend terms in the system based on the bends attribute
            of the system instance. All bend terms are given harmonic
            potentials either in the angle or the cosine of the angle.
        '''
        with log.section('VAL', 3, 'Initializing'):
            #list of bends which should be excluded
            if self.settings.excl_bends is not None:
                excl_bends = self.settings.excl_bends.split(',')

            #get the angle terms
            ffatypes = [self.system.ffatypes[fid] for fid in self.system.ffatype_ids]
            angles = {}
            for angle in self.system.iter_angles():
                angle, types = term_sort_atypes(ffatypes, angle, 'dihedral')
                if types in list(angles.keys()):
                    angles[types].append(angle)
                else:
                    angles[types] = [angle]
            #loop over all distinct angle types
            nabends = 0
            ncbends = 0
            nsqbends = 0
            for types, bends in angles.items():
                potkind = None
                rvs = []
                for i, bend in enumerate(bends):
                    if self.system.cell.nvec>0:
                        d10 = self.system.pos[bend[0]] - self.system.pos[bend[1]]
                        d12 = self.system.pos[bend[2]] - self.system.pos[bend[1]]
                        self.system.cell.mic(d10)
                        self.system.cell.mic(d12)
                        rvs.append(_bend_angle_low(d10, d12, 0)[0])
                    else:
                        rs = np.array([self.system.pos[j] for j in bend])
                        rvs.append(bend_angle(rs)[0])
                #sort rvs in rvs in [90-thresshold, 90+thresshold], rvs in
                #[180-thresshold,180] and others
                rvs90 = [rv for rv in rvs if 90*deg-thresshold<rv<90*deg+thresshold]
                rvs180 = [rv for rv in rvs if 180*deg-thresshold<rv<180*deg]
                rvsother = [rv for rv in rvs if rv not in rvs90 and rv not in rvs180]
                #detect whether rvs are centered around 90 and 180, then
                #use 1-cos(4*theta) term
                if len(rvs90)>0 and len(rvs180)>0 and len(rvsother)==0:
                    log.dump('%s has equilibrium values around 90/180 deg, using 1-cos(4*theta) potential' %('.'.join(types)))
                    potkind='squarebend'
                elif len(rvs180)>0 and len(rvs90)==0 and len(rvsother)==0:
                    log.dump('%s has equilibrium values around 180 deg, using 1+cos(theta) potential' %('.'.join(types)))
                    potkind='lincos'
                else:
                    potkind='angleharm'
                #add term
                for bend in bends:
                    skip = False

                    if self.settings.excl_bends is not None:
                        bend_opt1 = '.'.join(types)
                        bend_opt2 = '.'.join(types[::-1])
                        for excl in excl_bends:
                            pattern = re.compile(excl, re.IGNORECASE)
                            if pattern.match(bend_opt1) or pattern.match(bend_opt2):
                                skip = True

                    if not skip:
                        if potkind=='lincos':
                            basename = 'BendCheby1/'+'.'.join(types)
                            term = self.add_term(Chebychev1, [BendCos(*bend)], basename, ['HC_FC_DIAG'], ['kjmol', 'au'])
                            self.set_params(term.index, sign=1)
                            ncbends += 1
                        elif potkind=='squarebend':
                            basename = 'BendCheby4/'+'.'.join(types)
                            term = self.add_term(Chebychev4, [BendCos(*bend)], basename, ['HC_FC_DIAG'], ['kjmol', 'au'])
                            self.set_params(term.index, sign=-1)
                            nsqbends += 1
                        elif potkind=='angleharm':
                            basename = self.settings.bend_term+'/'+'.'.join(types)
                            if self.settings.bend_term.lower()   == 'bendaharm':
                                pot = Harmonic
                            elif self.settings.bend_term.lower() == 'bendmm3':
                                pot = MM3Bend
                            else:
                                raise ValueError('Bond kind %s not supported' %self.settings.bend_term)
                            self.add_term(pot, [BendAngle(*bend)], basename, ['PT_ALL', 'HC_FC_DIAG'], ['kjmol/rad**2', 'deg'])
                            nabends += 1
                        else:
                            raise ValueError('')
                    else:
                        log.dump('Excluded %s bend'%'.'.join(types))
        log.dump('Added %i bend terms (an)harmonic in the angle, %i bend terms with 1+cos(angle) potential and %i bend terms with 1-cos(4*theta) potential.' %(nabends, ncbends, nsqbends))

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
            #list of dihedrals which should be excluded
            if self.settings.excl_dihs is not None:
                excl_dihs = self.settings.excl_dihs.split(',')

            #get all dihedrals
            ffatypes = [self.system.ffatypes[fid] for fid in self.system.ffatype_ids]
            dihedrals = {}
            for dihedral in self.system.iter_dihedrals():
                dihedral, types = term_sort_atypes(ffatypes, dihedral, 'dihedral')
                if types in list(dihedrals.keys()):
                    dihedrals[types].append(dihedral)
                else:
                    dihedrals[types] = [dihedral]
            #loop over all distinct dihedral types
            ncheb = 0
            ncos = 0
            for types, diheds in dihedrals.items():
                psi0s = np.zeros(len(diheds), float)
                ms = np.zeros(len(diheds), float)
                bendskip = False
                for i, dihed in enumerate(diheds):
                    if self.system.cell.nvec>0:
                        d10 = self.system.pos[dihed[0]] - self.system.pos[dihed[1]]
                        d12 = self.system.pos[dihed[2]] - self.system.pos[dihed[1]]
                        d23 = self.system.pos[dihed[3]] - self.system.pos[dihed[2]]
                        self.system.cell.mic(d10)
                        self.system.cell.mic(d12)
                        self.system.cell.mic(d23)
                        #check if bending angle is not 180 deg
                        bend012 = _bend_angle_low(d10, d12, 0)[0]
                        bend123 = _bend_angle_low(d12, d23, 0)[0]
                        if bend012>175*deg or bend123>175*deg:
                            bendskip = True
                            continue
                        psi0s[i] = _dihed_angle_low(d10, d12, d23, 0)[0]
                    else:
                        rs = np.array([self.system.pos[j] for j in dihed])
                        bend012 = bend_angle(rs[:3])[0]
                        bend123 = bend_angle(rs[1:4])[0]
                        if bend012>175*deg or bend123>175*deg:
                            bendskip = True
                            continue
                        psi0s[i] = dihed_angle(rs)[0]
                    n1 = len(self.system.neighs1[dihed[1]])
                    n2 = len(self.system.neighs1[dihed[2]])
                    ms[i] = get_multiplicity(n1, n2)
                if bendskip:
                    log.warning('Dihedral for %s contains bend angle close to 180 deg, skipping' %('.'.join(types)))
                    continue
                if np.isnan(ms).any() or ms.std()>1e-3:
                    ms_string = str(ms)
                    if np.isnan(ms).any():
                        ms_string = 'nan'
                    log.warning('missing dihedral for %s (m is %s)' %('.'.join(types), ms_string))
                    continue
                m = int(np.round(ms.mean()))
                rv = get_restvalue(psi0s, m, thresshold=thresshold, mode=1)
                if rv is not None:
                    do_chebychev = True
                    chebypot = None
                    if m==1:
                        chebypot = Chebychev1
                        if   abs(rv)          <1e-6: sign = -1
                        elif abs(rv-180.0*deg)<1e-6: sign = 1
                        else: do_chebychev = False
                    elif m==2:
                        chebypot = Chebychev2
                        if   abs(rv)         <1e-6: sign = -1
                        elif abs(rv-90.0*deg)<1e-6: sign = 1
                        else: do_chebychev = False
                    elif m==3:
                        chebypot = Chebychev3
                        if   abs(rv)         <1e-6: sign = -1
                        elif abs(rv-60.0*deg)<1e-6: sign = 1
                        else: do_chebychev = False
                    elif m==4:
                        chebypot = Chebychev4
                        if   abs(rv)         <1e-6: sign = -1
                        elif abs(rv-45.0*deg)<1e-6: sign = 1
                        else: do_chebychev = False
                    elif m==6:
                        chebypot = Chebychev6
                        if   abs(rv)         <1e-6: sign = -1
                        elif abs(rv-30.0*deg)<1e-6: sign = 1
                        else: do_chebychev = False
                    else:
                        do_chebychev = False
                    for dihed in diheds:
                        skip = False

                        if self.settings.excl_dihs is not None:
                            dih_opt1 = '.'.join(types)
                            dih_opt2 = '.'.join(types[::-1])
                            for excl in excl_dihs:
                                pattern = re.compile(excl, re.IGNORECASE)
                                if pattern.match(dih_opt1) or pattern.match(dih_opt2):
                                    skip = True

                        if not skip:
                            if do_chebychev:
                                assert chebypot is not None
                                basename = 'TorsCheby%i/' %m+'.'.join(types)
                                term = self.add_term(chebypot, [DihedCos(*dihed)], basename, ['HC_FC_DIAG'], ['kjmol', 'au'])
                                self.set_params(term.index, sign=sign)
                                ncheb += 1
                            else:
                                basename = 'Torsion/'+'.'.join(types)
                                term = self.add_term(Cosine, [DihedAngle(*dihed)], basename, ['HC_FC_DIAG'], ['au', 'kjmol', 'deg'])
                                self.set_params(term.index, rv0=rv, m=m)
                                ncos += 1
                        else:
                            log.dump('Excluded %s dihedral'%'.'.join(types))
                else:
                    #no dihedral potential could be determine, hence it is ignored
                    log.warning('missing dihedral for %s (could not determine rest value from %s)' %('.'.join(types), str(psi0s/deg)))
                    continue
            log.dump('Added %i Cosine dihedral terms (of which %i are described using Chebychev terms)' %(ncos+ncheb, ncheb))

    def init_oop_terms(self, thresshold_zero=5e-2*angstrom):
        #TODO: make settings option ofo thresshold_zero
        '''
            Initialize all out-of-plane terms in the system based on the oops
            attribute of the system instance. All oops are given harmonic
            potentials.
        '''
        with log.section('VAL', 3, 'Initializing'):
            #list of oopds which should be excluded
            if self.settings.excl_oopds is not None:
                excl_oopds = self.settings.excl_oopds.split(',')

            #get all out-of-plane distances
            from molmod.ic import opbend_dist, _opdist_low
            ffatypes = [self.system.ffatypes[fid] for fid in self.system.ffatype_ids]
            opdists = {}
            for opdist in self.system.iter_oops():
                opdist, types = term_sort_atypes(ffatypes, opdist, 'opdist')
                if types in list(opdists.keys()):
                    opdists[types].append(opdist)
                else:
                    opdists[types] = [opdist]
            #loop over all distinct opdist types
            nharm = 0
            nsq = 0
            for types, oops in opdists.items():
                skip = False

                if self.settings.excl_oopds is not None:
                    oopd_opt1 = '.'.join(types)
                    oopd_opt2 = '.'.join([types[0],types[2],types[1],types[3]])
                    oopd_opt3 = '.'.join([types[1],types[0],types[2],types[3]])
                    oopd_opt4 = '.'.join([types[1],types[2],types[0],types[3]])
                    oopd_opt5 = '.'.join([types[2],types[0],types[1],types[3]])
                    oopd_opt6 = '.'.join([types[2],types[1],types[0],types[3]])
                    for excl in excl_oopds:
                        pattern = re.compile(excl, re.IGNORECASE)
                        if pattern.match(oopd_opt1) or pattern.match(oopd_opt2) or pattern.match(oopd_opt3) or pattern.match(oopd_opt4) or pattern.match(oopd_opt5) or pattern.match(oopd_opt6):
                            skip = True
                if not skip:
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
                            basename = 'Oopdist/'+'.'.join(types)
                            term = self.add_term(Harmonic, [OopDist(*oop)], basename, ['HC_FC_DIAG'], ['kjmol/A**2', 'A'])
                            self.set_params(term.index, rv0=0.0)
                            nharm += 1
                    else:
                        #add term harmonic in square of oopdist
                        log.dump('Mean absolute value of OopDist %s is %.3e A, used SQOOPDIST' %('.'.join(types), d0s.mean()/angstrom))
                        for oop in oops:
                            basename = 'SqOopdist/'+'.'.join(types)
                            self.add_term(Harmonic, [SqOopDist(*oop)], basename, ['PT_ALL', 'HC_FC_DIAG'], ['kjmol/A**4', 'A**2'])
                            nsq += 1
                else:
                    log.dump('Excluded %s oopd'%'.'.join(types))
        log.dump('Added %i Harmonic and %i SquareHarmonic out-of-plane distance terms' %(nharm, nsq))

    def get_term_index(self, label):
        # Find index of master term(s) matching given label
        pattern = re.compile(label, re.IGNORECASE)
        candidates = []
        for iterm, term in enumerate(self.terms):
            if pattern.match(term.basename) and term.is_master():
                candidates.append(iterm)
        assert len(candidates)<2, 'Multiple masters found for %s: %s' %(label, ','.join([self.terms[iterm].basename for iterm in candidates]))
        if len(candidates)==0:
            if label.startswith('^Bond') or label.startswith('^Bend') or label.startswith('^Tors'):
                sublabels = label.split('|')
                prefix0, types0 = sublabels[0].lstrip('^').rstrip('$').split('/')
                label  = '^'+prefix0+'/'+'\.'.join(types0.split('.')[::-1])+'$'
                if len(sublabels)>1:
                    prefix1, types1 = sublabels[1].lstrip('^').rstrip('$').split('/')
                    label += '|'
                    label += '^'+prefix1+'/'+'\.'.join(types1.split('.')[::-1])+'$'
                pattern = re.compile(label, re.IGNORECASE)
                for iterm, term in enumerate(self.terms):
                    if pattern.match(term.basename) and term.is_master():
                        candidates.append(iterm)
            assert len(candidates)<2, 'Multiple masters found for %s: %s' %(label, ','.join([self.terms[iterm].basename for iterm in candidates]))
        assert (len(candidates))==1, "Could not find term for %s" % (label)
        return candidates[0]


    def init_cross_angle_terms(self):
        '''
            Initialize cross terms between bonds and bends.
        '''
        with log.section('VAL', 3, 'Initializing'):
            #list of bonds which should be excluded
            if self.settings.excl_bonds is not None:
                excl_bonds = self.settings.excl_bonds.split(',')

            #list of bends which should be excluded
            if self.settings.excl_bends is not None:
                excl_bends = self.settings.excl_bends.split(',')

            ffatypes = [self.system.ffatypes[i] for i in self.system.ffatype_ids]
            #add cross terms for angle patterns
            nss = 0
            nsa = 0
            for angle in self.system.iter_angles():
                skip = False
                angle, types = term_sort_atypes(ffatypes, angle, 'angle')
                anglekind = None

                if self.settings.excl_bends is not None:
                    bend_opt1 = '.'.join(types)
                    bend_opt2 = '.'.join(types[::-1])
                    for excl in excl_bends:
                        pattern = re.compile(excl, re.IGNORECASE)
                        if pattern.match(bend_opt1) or pattern.match(bend_opt2):
                            skip = True

                if not skip:
                    for term in self.iter_masters('^.*/'+'\.'.join(types)+'$', use_re=True):
                        if len(term.get_atoms())!=3: continue
                        if len(term.ics)>1: continue
                        assert anglekind is None, '2 masters detected for angle %s' %('.'.join(types))
                        anglekind = term.ics[0].kind
                    assert anglekind is not None, 'No master found for angle %s' %('.'.join(types))
                    bond0, btypes0 = term_sort_atypes(ffatypes, angle[:2], 'bond')
                    bond1, btypes1 = term_sort_atypes(ffatypes, angle[1:], 'bond')

                    if self.settings.excl_bonds is not None:
                        bond_opt1 = '.'.join(btypes0)
                        bond_opt2 = '.'.join(btypes0[::-1])
                        bond_opt3 = '.'.join(btypes1)
                        bond_opt4 = '.'.join(btypes1[::-1])
                        for excl in excl_bonds:
                            pattern = re.compile(excl, re.IGNORECASE)
                            if pattern.match(bond_opt1) or pattern.match(bond_opt2) or pattern.match(bond_opt3) or pattern.match(bond_opt4):
                                skip = True

                    if not skip:
                        # Find indexes of diagonal terms corresponding to ics
                        # appearing here
                        diag_term_indexes = []
                        for btypes in [btypes0,btypes1]:
                            label = '^Bond.*/%s$' %('.'.join(btypes))
                            diag_term_indexes.append(self.get_term_index(label))
                        label = '^Bend.*/%s$' %('.'.join(types))
                        diag_term_indexes.append(self.get_term_index(label))
                        #add stretch-stretch
                        if self.settings.do_cross_ASS:
                            self.add_term(
                                Cross, [Bond(*bond0), Bond(*bond1)],
                                'Cross/'+'.'.join(types)+'/bb', ['HC_FC_CROSS_ASS'], ['kjmol/A**2', 'A', 'A'],
                                diag_term_indexes=diag_term_indexes[:2],
                            )
                            nss += 1
                        #add stretch-bends
                        if self.settings.do_cross_ASA:
                            if anglekind == 2:
                                basename = 'Cross/'+'.'.join(types)
                                ic = BendAngle(*angle)
                                unit = 'deg'
                            #elif anglekind == 1:
                            #TODO: this is switched off, does not make much difference and
                            #gives issues in the case of diagonal BendCheby4 terms
                            #    basename = 'CrossCBend/'+'.'.join(types)
                            #    ic = BendCos(*angle)
                            #    unit = 'au'
                            else:
                                log.dump('Skipped stretch-angle cross term for %s due to incompatible diagonal bend term with ickind=%i' %('.'.join(types), anglekind))
                                continue
                            self.add_term(
                                Cross, [Bond(*bond0), ic],
                                basename+'/b0a', ['HC_FC_CROSS_ASA'], ['kjmol/A', 'A', unit],
                                diag_term_indexes=[diag_term_indexes[0],diag_term_indexes[2]],
                            )
                            self.add_term(
                                Cross, [Bond(*bond1), ic],
                                basename+'/b1a', ['HC_FC_CROSS_ASA'], ['kjmol/A', 'A', unit],
                                diag_term_indexes=[diag_term_indexes[1],diag_term_indexes[2]],
                            )
                            nsa += 2
            log.dump('Added %i stretch-stretch and %i stretch-angle cross terms from angle patterns' %(nss, nsa))

    def init_cross_dihed_terms(self):
        '''
            Initialize cross terms between diheds and bonds,bends.
        '''
        from yaff.pes.iclist import DihedCos2, DihedCos3, DihedCos4, DihedCos6

        with log.section('VAL', 3, 'Initializing'):
            ffatypes = [self.system.ffatypes[i] for i in self.system.ffatype_ids]
            #add cross terms for dihedral patterns
            nss = 0
            nsd = 0
            naa = 0
            nad = 0
            for dihed in self.system.iter_dihedrals():
                bond01, btypes01 = term_sort_atypes(ffatypes, dihed[0:2], 'bond')
                bond12, btypes12 = term_sort_atypes(ffatypes, dihed[1:3], 'bond')
                bond23, btypes23 = term_sort_atypes(ffatypes, dihed[2:4], 'bond')
                angle012, atypes012 = term_sort_atypes(ffatypes, dihed[0:3], 'angle')
                angle123, atypes123 = term_sort_atypes(ffatypes, dihed[1:4], 'angle')
                dihed, types = term_sort_atypes(ffatypes, dihed, 'dihedral')
                # Find indexes of diagonal terms corresponding to ics
                # appearing here
                diag_term_indexes = []
                for btypes in [btypes01,btypes12,btypes23]:
                    label = '^Bond.*/%s$' %('.'.join(btypes))
                    diag_term_indexes.append(self.get_term_index(label))
                for atypes in [atypes012,atypes123]:
                    label = '^Bend.*/%s$' %('.'.join(atypes))
                    diag_term_indexes.append(self.get_term_index(label))
                label = '^Tors.*/%s$' %('.'.join(types))
                diag_term_indexes.append(self.get_term_index(label))
                #get multiplicity of dihedral term to determine which DihedCos
                m, DihedIC = None, None
                for term in self.iter_masters('^.*/'+'\.'.join(types)+'$', use_re=True):
                    if term.kind == 4:
                        assert DihedIC is None and m is None
                        m=self.get_params(term.index, only='m')
                        if   m==1: DihedIC = DihedCos
                        elif m==2: DihedIC = DihedCos2
                        elif m==3: DihedIC = DihedCos3
                        elif m==4: DihedIC = DihedCos4
                        elif m==6: DihedIC = DihedCos6
                    elif term.kind == 5:
                        assert DihedIC is None and m is None
                        m, DihedIC = 1, DihedCos
                    elif term.kind == 6:
                        assert DihedIC is None and m is None
                        m, DihedIC = 2, DihedCos2
                    elif term.kind == 7:
                        assert DihedIC is None and m is None
                        m, DihedIC = 3, DihedCos3
                    elif term.kind == 8:
                        assert DihedIC is None and m is None
                        m, DihedIC = 4, DihedCos4
                    elif term.kind == 9:
                        assert DihedIC is None and m is None
                        m, DihedIC = 6, DihedCos6
                if m is None or DihedIC is None:
                    log.dump('No multiplicity found for %s, skipping' %'.'.join(types))
                    continue
                #get type of angle012 and angle123
                angle012_type = None
                angle123_type = None
                for term in self.iter_masters('^.*/'+'\.'.join(atypes012)+'$', use_re=True):
                    if len(term.get_atoms())!=3: continue
                    if len(term.ics)>1: continue
                    assert angle012_type is None, 'Two masters found for angle %s' %(str(types[:4]))
                    angle012_type = term.ics[0].kind
                for term in self.iter_masters('^.*/'+'\.'.join(atypes123)+'$', use_re=True):
                    if len(term.get_atoms())!=3: continue
                    if len(term.ics)>1: continue
                    assert angle123_type is None, 'Two masters found for angle %s' %(str(types[1:]))
                    angle123_type = term.ics[0].kind

                #add stretch-stretch term:
                if self.settings.do_cross_DSS:

                    if self.settings.excl_bonds is not None:
                        raise NotImplementedError
                    if self.settings.excl_bends is not None:
                        raise NotImplementedError
                    if self.settings.excl_dihs is not None:
                        raise NotImplementedError

                    basename = 'CrossBondDih%i/'%m+'.'.join(types)
                    self.add_term(
                        Cross, [Bond(*bond01), Bond(*bond23)],
                        basename+'/bb', ['HC_FC_CROSS_DSS'], ['kjmol/A**2', 'A', 'A'],
                        diag_term_indexes=[diag_term_indexes[0],diag_term_indexes[2]]
                    )
                    nss += 1
                #add stretch-dihedral terms:
                if self.settings.do_cross_DSD:
                    basename = 'CrossBondDih%i/'%m+'.'.join(types)
                    self.add_term(
                        Cross, [Bond(*bond01), DihedIC(*dihed)],
                        basename+'/b0d', ['HC_FC_CROSS_DSD'], ['kjmol/A', 'A', 'au'],
                        diag_term_indexes=[diag_term_indexes[0],diag_term_indexes[5]]
                    )
                    self.add_term(
                        Cross, [Bond(*bond12), DihedIC(*dihed)],
                        basename+'/b1d', [], ['kjmol/A', 'A', 'au'],
                        diag_term_indexes=[diag_term_indexes[1],diag_term_indexes[5]]
                    )
                    self.add_term(
                        Cross, [Bond(*bond23), DihedIC(*dihed)],
                        basename+'/b2d', ['HC_FC_CROSS_DSD'], ['kjmol/A', 'A', 'au'],
                        diag_term_indexes=[diag_term_indexes[2],diag_term_indexes[5]]
                    )
                    nsd += 3
                #add angle-angle term:
                if self.settings.do_cross_DAA:

                    if self.settings.excl_bonds is not None:
                        raise NotImplementedError
                    if self.settings.excl_bends is not None:
                        raise NotImplementedError
                    if self.settings.excl_dihs is not None:
                        raise NotImplementedError

                    assert angle012_type is not None, 'No master found for angle012 in %s' %('.'.join(types))
                    assert angle123_type is not None, 'No master found for angle123 in %s' %('.'.join(types))
                    #add angle-angle term:
                    if angle012_type == 2 and angle123_type == 2:
                        self.add_term(
                            Cross, [BendAngle(*angle012), BendAngle(*angle123)],
                            'CrossBendDih%i/'%m+'.'.join(types)+'/aa',
                            ['HC_FC_CROSS_DAA'], ['kjmol', 'deg', 'deg'],
                            diag_term_indexes=[diag_term_indexes[3],diag_term_indexes[4]]
                        )
                        nad += 1
                    else:
                        log.dump('Skipped angle-angle cross term for %s due to incompatible bend kinds' %('.'.join(types)))
                        continue

                if self.settings.do_cross_DAD:
                    #add angle-dihedral terms:
                    if angle012_type == 2:
                        self.add_term(
                            Cross, [BendAngle(*angle012), DihedIC(*dihed)],
                            'CrossBendDih%i/'%m+'.'.join(types)+'/a0d',
                            ['HC_FC_CROSS_DAD'], ['kjmol', 'deg', 'au'],
                            diag_term_indexes=[diag_term_indexes[3],diag_term_indexes[5]]
                        )
                        nad += 1
                    elif angle012_type == 1:
                        self.add_term(
                            Cross, [BendCos(*angle012), DihedIC(*dihed)],
                            'CrossCBendDih%i/'%m+'.'.join(types)+'/a0d',
                            ['HC_FC_CROSS_DAD'], ['kjmol', 'au', 'au'],
                            diag_term_indexes=[diag_term_indexes[4],diag_term_indexes[5]]
                        )
                        nad += 1
                    else:
                        log.dump('Skipped angle-dihedral cross term for %s due to incompatible diagonal bend term with ickind=%i' %('.'.join(types), angle012_type))
                        continue
                    if angle123_type == 2:
                        self.add_term(
                            Cross, [BendAngle(*angle123), DihedIC(*dihed)],
                            'CrossBendDih%i/'%m+'.'.join(types)+'/a1d',
                            ['HC_FC_CROSS_DAD'], ['kjmol', 'deg', 'au'],
                            diag_term_indexes=[diag_term_indexes[3],diag_term_indexes[5]]
                        )
                        nad += 1
                    elif angle123_type == 1:
                        self.add_term(
                            Cross, [BendCos(*angle123), DihedIC(*dihed)],
                            'CrossCBendDih%i/'%m+'.'.join(types)+'/a1d',
                            ['HC_FC_CROSS_DAD'], ['kjmol', 'au', 'au'],
                            diag_term_indexes=[diag_term_indexes[4],diag_term_indexes[5]]
                        )
                        nad += 1
                    else:
                        log.dump('Skipped angle-dihedral cross term for %s due to incompatible diagonal bend term with ickind=%i' %('.'.join(types), angle132_type))
                        continue
        log.dump('Added %i stretch-stretch, %i stretch-dihedral, %i angle-angle and %i angle-dihedral cross terms from dihedral patterns' %(nss, nsd, naa, nad))

    def apply_constraints(self, constraints):
        '''
            Routine to apply equality constraints in the force field fitting

            **Arguments**

            contraints
                A dictionairy containing (master: slaves) definitions in which
                master is a string defining the master basename and slaves is a
                list of strings defining the slave basenames.
        '''
        for mastername, slavenames in constraints.items():
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
        if kind in [5,6,7,8,9]:#Chebychev
            potentials={5: Chebychev1, 6: Chebychev2, 7: Chebychev3, 8: Chebychev4, 9: Chebychev6}
            k, sign = self.get_params(index, only='all')
            if fc is not None: k = fc
            for jterm in masterslaves:
                ics = self.terms[jterm].ics
                pot = potentials[kind]
                args = (k,) + tuple(ics)
                val.add_term(pot(*args,sign=sign))
        elif kind==4:#Cosine
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
        elif kind in [0,2,11,12]:#[Harmonic,Fues,MM3Quartic,MM3Bend]
            potentials={0:Harmonic,2:Fues,11:MM3Quartic,12:MM3Bend}
            k, rv = self.get_params(index)
            if fc is not None: k = fc
            for jterm in masterslaves:
                ics = self.terms[jterm].ics
                args = (k, rv) + tuple(ics)
                val.add_term(potentials[kind](*args))
        else:
            raise ValueError('Term kind %i not supported' %kind)
        ff = ForceField(self.system, [val])
        hcov = estimate_cart_hessian(ff)
        return hcov

    def set_params(self, term_index, fc=None, rv0=None, rv1=None, m=None,
            a0=None, a1=None, a2=None, a3=None, sign=None, ediss=None, exp=None):
        term = self.vlist.vtab[term_index]
        if term['kind'] in [0,2,11,12]:#['Harmonic', 'Fues', 'MM3Quartic', 'MM3Bend']
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
        elif term['kind'] in [5,6,7,8,9]: #Chebychevs
            if fc is not None: term['par0'] = fc
            if sign is not None: term['par1'] = sign
        elif term['kind'] in [3]:#['Cross']
            if fc is not None:  term['par0'] = fc
            if rv0 is not None: term['par1'] = rv0
            if rv1 is not None: term['par2'] = rv1
        elif term['kind'] in [14]: #Morse
            if ediss is not None: term['par0'] = ediss
            if fc is not None:
                if exp is not None: raise IOError('When fc is set in Morse, exp cannot be set, but will be adapted to Ediss and fc')
                term['par1'] = np.sqrt(fc/(2.0*term['par0']))
            if exp is not None:
                if fc is not None: raise IOError('When exp is set in Morse, fc cannot be set, but will be adapted to Ediss and exp')
                term['par1'] = exp
            if rv0 is not None:   term['par2'] = rv0
        else:
            raise NotImplementedError('set_params not implemented for Yaff %s term' %term['kind'])

    def get_params(self, term_index, only='all'):
        term = self.vlist.vtab[term_index]
        if term['kind'] in [0,2,11,12]:#['Harmonic', 'Fues', 'MM3Quartic', 'MM3Bend']
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
        elif term['kind'] in [5,6,7,8,9]: #Chebychevs
            if only.lower()=='all': return term['par0'], term['par1'],
            elif only.lower()=='fc': return term['par0']
            elif only.lower()=='sign': return term['par1']
            else: raise ValueError('Invalid par kind definition %s' %only)
        elif term['kind'] in [3]:#['Cross']
            if only.lower()=='all': return term['par0'], term['par1'], term['par2']
            elif only.lower()=='fc': return term['par0']
            elif only.lower()=='rv0': return term['par1']
            elif only.lower()=='rv1': return term['par2']
            else: raise ValueError('Invalid par kind definition %s' %only)
        elif term['kind'] in [14]: #Morse
            if only.lower()=='all': return term['par0'], term['par1'], term['par2']
            elif only.lower()=='ediss': return term['par0']
            elif only.lower()=='exp': return term['par1']
            elif only.lower()=='rv0': return term['par2']
            elif only.lower()=='fc': return 2.0*term['par0']*term['par1']**2
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
        if term['kind'] in [0, 2, 11, 12]:  # ['Harmonic', 'Fues', 'MM3Quartic', 'MM3Bend']
            return abs(term['par0']) < 1e-6*kjmol
        elif term['kind'] in [14]: #['Morse']
            return abs(2.0*term['par0']*term['par1']**2) < 1e-6*kjmol
        elif term['kind'] in [1,13]:   # ['PolyFour', 'BondDoubleWell']
            # Not sure how to handle this one...
            # For now, never neglect.
            return False
        elif term['kind'] in [4]:   # ['Cosine']
            return abs(term['par1']) < 1e-6*kjmol
        elif term['kind'] in [3]:   # ['Cross']
            return abs(term['par0']) < 1e-6*kjmol
        elif term['kind'] in [5,6,7,8,9]: #chebychev
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
                'bondharm', 'bondmm3', 'bondfues',
                'bendaharm', 'bendmm3',
                'bendcheby1', 'bendcheby4', 'bendcharm',
                'torscheby', 'torsion', 'torsc2harm', 'dihedharm',
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
