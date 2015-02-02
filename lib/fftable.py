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

from molmod.units import deg, angstrom, kjmol, parse_unit, rad
import numpy as np

from quickff.tools import statistics


__all__ = ['DataArray', 'FFTable']


class DataArray(object):
    '''
        A class to store data representing measurements (or calculations) of
        the same observable. Some basic statistics will be done automatically
        and can be accessed at all time.
    '''
    def __init__(self, data=None, unit='au'):
        if data is None:
            self.data = None
        else:
            self.data = np.array(data)
        self.unit = unit
        self.mean, self.std, self.num = statistics(self.data)

    def append(self, value):
        'Add value to array and update statistics'
        if self.data is None:
            self.data = np.array([value])
        else:
            assert isinstance(self.data, np.ndarray)
            self.data = np.append(self.data, value)
        self.mean, self.std, self.num = statistics(self.data)

    def __len__(self):
        'Get the number of elements in the array'
        return len(self.data)

    def string(self):
        'Convert the statistics of the array into a human readable string'
        if self.num == 0:
            return ''
        else:
            s = u'%9.3f \u00B1 %9.3f ' % (
                self.mean/parse_unit(self.unit),
                self.std/parse_unit(self.unit)) + \
                self.unit + ' '*(15-len(self.unit)
            )
            return s.encode('utf-8')


class FFTable(object):
    '''
        A class to read, store and dump Force Field parameters. The parameters
        of the force field term for an ic of kind icname can be accessed as
        follows:

            k, q0 = fftable[icname]

        The statistics of a parameter (e.q. the force constant of the term
        related to icname) can be accessed as follows:

            k_mean = fftable.pars[icname]['k'].mean
            k_std = fftable.pars[icname]['k'].std
            k_num = fftable.pars[icname]['k'].num
    '''
    def __init__(self):
        self.pars = {}

    def add(self, icname, ks, q0s, **kwargs):
        'Add force field parameters for ic with name icname'
        assert isinstance(ks, DataArray)
        assert isinstance(q0s, DataArray)
        self.pars[icname] = {'k': ks, 'q0': q0s}
        for kw, data in kwargs.iteritems():
            assert isinstance(data, DataArray)
            self.pars[icname][kw] = data

    def __getitem__(self, icname):
        'Get the mean force constant and restvalue of ic with name icname'
        k  = self.pars[icname]['k'].mean
        q0 = self.pars[icname]['q0'].mean
        return k, q0

    def print_screen(self):
        'Print force field parameters of all ics to a human readable table'
        maxlength = max([len(icname) for icname in self.pars.keys()]) + 2
        for icname, pars in sorted(self.pars.iteritems()):
            descr = icname + ' '*(maxlength-len(icname))
            print '    %s   K = %s    q0 = %s' % (
                descr, pars['k'].string(), pars['q0'].string()
            )

    def dump_ffit2(self, fn, mode='w'):
        'Dump force field parameters to a file in FFit2 format'
        f = open(fn, mode)
        print >> f, '# Parameters'
        print >> f, '# -----------------------------------------------------------------------------#------'
        print >> f, '# longname                                                   unit        value # fx/fr'
        print >> f, '# -----------------------------------------------------------------------------#------'
        for icname in sorted(self.pars.keys()):
            kind = icname.split('/')[0]
            atypes = icname.split('/')[1]
            k, q0 = self[icname]
            if kind == 'bond':
                name = 'bond/%s/harm/dist/K' % atypes
                print >> f, '  %50s %12s % 12.6f #  free' % (
                     name + ' '*(50-len(name)),
                    'kjmol/A^2', k/(kjmol/angstrom**2)
                )
                name = 'bond/%s/harm/dist/q0' % atypes
                print >> f, '  %50s %12s % 12.6f #  free' % (
                     name + ' '*(50-len(name)),
                    'A', q0/angstrom
                )
            elif kind == 'angle':
                name = 'angle/%s/harm/angle/K' % atypes
                print >> f, '  %50s %12s % 12.6f #  free' % (
                    name + ' '*(50-len(name)),
                    'kjmol/rad^2', k/(kjmol/rad**2)
                )
                name = 'angle/%s/harm/angle/q0' % atypes
                print >> f, '  %50s %12s % 12.6f #  free' % (
                    name + ' '*(50-len(name)),
                    'deg', q0/deg
                )
            elif kind == 'dihed':
                if 'm' in self.pars[icname].keys():
                    m = self.pars[icname]['m'].mean
                    name = 'dihed/%s/cos-m%i-0/dihed/K' % (atypes, m)
                    print >> f, '  %50s %12s % 12.6f #  free' % (
                        name + ' '*(50-len(name)),
                        'kjmol', k/kjmol
                    )
                else:
                    name = 'dihed/%s/harm/cdihed' % (atypes)
                    print >> f, '  %50s %12s % 12.6f #  free' % (
                        name+'/K' + ' '*(48-len(name)),
                        'kjmol', k/kjmol
                    )
                    print >> f, '  %50s %12s % 12.6f #  free' % (
                        name+'/q0' + ' '*(47-len(name)),
                        'au', q0
                    )
            elif kind == 'opdist':
                name = 'opdist/%s/harm/opdist/K' % atypes
                print >> f, '  %50s %12s % 12.6f #  free' % (
                     name + ' '*(50-len(name)),
                    'kjmol/A^2', k/(kjmol/angstrom**2)
                )
                name = 'opdist/%s/harm/opdist/q0' % atypes
                print >> f, '  %50s %12s % 12.6f #  free' % (
                     name + ' '*(50-len(name)),
                    'A', q0/angstrom
                )
        print >> f, '# -----------------------------------------------------------------------------#------'
        f.close()

    def dump_yaff(self, fn, mode='w'):
        'Dump force field parameters to a file in Yaff format'
        f = open(fn, mode)
        print >> f, '# Bond stretch'
        print >> f, '# ============'
        print >> f, ''
        print >> f, '# Mathematical form depends on the kind selected below. Few kinds are supported:'
        print >> f, '# - BONDHARM: 0.5*K*(r-R0)**2'
        print >> f, '# - BONDFUES: 0.5*K*R0**2*(1+(R0/r)*((R0/r)-2.0))'
        print >> f, ''
        print >> f, '# The actual parameters and their units may depend on the kind.'
        print >> f, 'BONDHARM:UNIT K kjmol/angstrom**2'
        print >> f, 'BONDHARM:UNIT R0 angstrom'
        print >> f, ''
        print >> f, '# -----------------------------------------------------------------'
        print >> f, '# KEY         ffatype0 ffatype1  K                 R0'
        print >> f, '# -----------------------------------------------------------------'
        for icname in sorted(self.pars.keys()):
            if not icname.startswith('bond'): continue
            atypes = icname.split('/')[1].split('.')
            k, q0 = self[icname]
            print >> f, 'BONDHARM:PARS %8s %8s % .10e % .10e' % (
                atypes[0], atypes[1], k/(kjmol/angstrom**2), q0/angstrom
            )
        print >> f, ''
        print >> f, '# Angle bending'
        print >> f, '# ============='
        print >> f, ''
        print >> f, '# Mathematical form depends on the kind selected below. Few kinds are supported:'
        print >> f, '# - BENDAHARM: 0.5*K*(theta-THETA0)**2'
        print >> f, '# - BENDCHARM: 0.5*K*(cos(theta)-cos(THETA0))**2'
        print >> f, '# - UBHARM: 0.5*K*(r-R0)**2'
        print >> f, '# where theta is the bending angle and r is the distance between the non-bonded'
        print >> f, '# pair of atoms.'
        print >> f, ''
        print >> f, '# The actual parameters and their units may depend on the kind.'
        print >> f, 'BENDAHARM:UNIT K kjmol/rad**2'
        print >> f, 'BENDAHARM:UNIT THETA0 deg'
        print >> f, ''
        print >> f, '# ---------------------------------------------------------------------------'
        print >> f, '# KEY          ffatype0 ffatype1 ffatype2  K                 THETA0/COS0/R0'
        print >> f, '# ---------------------------------------------------------------------------'
        for icname in sorted(self.pars.keys()):
            if not icname.startswith('angle'): continue
            atypes = icname.split('/')[1].split('.')
            k, q0 = self[icname]
            if q0 > 180.0*deg:
                q0 = 2*np.pi - q0
            print >> f, 'BENDAHARM:PARS %8s %8s %8s % .10e % .10e' % (
                atypes[0], atypes[1], atypes[2], k/(kjmol/rad**2), q0/deg
            )
        print >> f, ''
        print >> f, '# Torsional terms'
        print >> f, '# ==============='
        print >> f, ''
        print >> f, '# The following mathemetical for is supported:'
        print >> f, '#  - TORSION:   0.5*A*(1-COS(M*(PHI-PHI0)))'
        print >> f, '#  - TORSCHARM: 0.5*A*(COS(PHI)-COS0)**2'
        print >> f, ''
        print >> f, '# The actual parameters and their units may depend on the kind.'
        print >> f, 'TORSION:UNIT A kjmol'
        print >> f, 'TORSION:UNIT PHI0 deg'
        print >> f, 'TORSCHARM:UNIT A kjmol'
        print >> f, 'TORSCHARM:UNIT COS0 au'
        print >> f, ''
        print >> f, '# -------------------------------------------------------------------------------------'
        print >> f, '# KEY          ffatype0 ffatype1 ffatype2 ffatype4  M  A                 PHI0/COS0'
        print >> f, '# -------------------------------------------------------------------------------------'
        for icname in sorted(self.pars.keys()):
            if not icname.startswith('dihed'): continue
            atypes = icname.split('/')[1].split('.')
            k, q0 = self[icname]
            if 'm' in self.pars[icname].keys():
                m = self.pars[icname]['m'].mean
                print >> f, 'TORSION:PARS   %8s %8s %8s %8s %2i % .10e % .10e' % (
                    atypes[0], atypes[1], atypes[2], atypes[3],
                    m, k/kjmol, q0/deg
                )
            else:
                print >> f, 'TORSCHARM:PARS %8s %8s %8s %8s    % .10e % .10e' % (
                    atypes[0], atypes[1], atypes[2], atypes[3],
                    k/kjmol, q0
                )
        print >> f, ''
        print >> f, '# Out-of-plane terms'
        print >> f, '# ==============='
        print >> f, ''
        print >> f, '# The following mathemetical for is supported:'
        print >> f, '#  - OPDIST: 0.5*K*(d - d0)^2'
        print >> f, ''
        print >> f, '# The actual parameters and their units may depend on the kind.'
        print >> f, 'OOPDIST:UNIT K kjmol/angstrom**2'
        print >> f, 'OOPDIST:UNIT D0 angstrom'
        print >> f, ''
        print >> f, '# -------------------------------------------------------------------------------------'
        print >> f, '# KEY        ffatype0 ffatype1 ffatype2 ffatype4  K                 D0'
        print >> f, '# -------------------------------------------------------------------------------------'
        for icname in sorted(self.pars.keys()):
            if not icname.startswith('opdist'): continue
            atypes = icname.split('/')[1].split('.')
            k, q0 = self[icname]
            print >> f, 'OOPDIST:PARS %8s %8s %8s %8s % .10e % .10e' % (
                atypes[0], atypes[1], atypes[2], atypes[3],
                k/(kjmol/angstrom**2), q0/angstrom
            )
        print >> f, ''
        f.close()
