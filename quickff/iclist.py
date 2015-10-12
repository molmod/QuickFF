# -*- coding: utf-8 -*-
# QuickFF is a code to quickly derive accurate force fields from ab initio input.
# Copyright (C) 2012 - 2015 Louis Vanduyfhuys <Louis.Vanduyfhuys@UGent.be>
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

'''List of the internal coordinates that appear in the covalent force-field expression.
'''


import numpy as np

from yaff import InternalCoordinateList, Bond, BendAngle, DihedAngle, OopDist

__all__ = ['ICList']

class ICList(InternalCoordinateList):
    '''
    A list of internal coordinates to which the QuickFF procedure will be
    applied. This is a child class from the Yaff Internal CoordinateList, which
    will handle all computations. The extra layer implemented here just adds
    a bit more information relevant to QuickFF.
    '''
    def __init__(self, dlist):
        self.icnames = []
        self.kunits = []
        self.qunits = []
        self.icname_ids = []
        self.ics = []
        super(ICList, self).__init__(dlist)

    def sort_ffatypes(self, ffatypes, kind=''):
        '''
        Sort the atom types of an internal coordinate. This is used to obtain
        a string that allows to distinguish types of internal coordinates.
        '''
        if kind != 'opdist':
            if ffatypes[0] < ffatypes[-1]:
                return ffatypes
            elif ffatypes[0]==ffatypes[-1]:
                if len(ffatypes)==4:
                    if ffatypes[1]<ffatypes[2]:
                        return ffatypes
                    else:
                        return ffatypes[::-1]
                else:
                    return ffatypes
            else:
                return ffatypes[::-1]
        else:
            result = sorted(ffatypes[:3])
            result.append(ffatypes[3])
            return result

    def finalize_icnames(self):
        # Turn some lists into arrays
        self.icnames = np.asarray(self.icnames)
        self.icname_ids = np.asarray(self.icname_ids)
        assert self.icname_ids.shape[0]==self.nic

    def add_ic(self, atoms, kind, skip_ics):
        '''
        Add an internal coordinate to the list. Next to adding it to the
        Yaff InternalCoordinateList, some more information about the ic is
        stored here.
        '''
        name = kind + '/'+'.'.join(self.sort_ffatypes([self.dlist.system.ffatypes[self.dlist.system.ffatype_ids[at]] for at in atoms], kind=kind))
        if kind=='bond':
            ic = Bond(atoms[0],atoms[1])
            qunit, kunit ='A', 'kjmol/A**2'
        elif kind=='angle':
            ic = BendAngle(atoms[0], atoms[1], atoms[2])
            qunit, kunit ='deg', 'kjmol/rad**2'
        elif kind=='dihed':
            ic = DihedAngle(atoms[0], atoms[1], atoms[2], atoms[3])
            qunit, kunit ='deg', 'kjmol'
        elif kind=='opdist':
            ic = OopDist(atoms[0], atoms[1], atoms[2], atoms[3])
            qunit, kunit ='A', 'kjmol/A**2'
        else: raise NotImplementedError
        if name in skip_ics: return
        if not name in self.icnames:
            self.icnames.append(name)
            self.kunits.append(kunit)
            self.qunits.append(qunit)
        self.icname_ids.append(self.icnames.index(name))
        super(ICList, self).add_ic(ic)
        self.ics.append(ic)
