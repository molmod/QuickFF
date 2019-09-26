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

from __future__ import print_function, absolute_import
from io import IOBase
from quickff.log import log
from molmod.units import parse_unit
import os

try:
    from importlib.resources import path
except ImportError:
    from importlib_resources import path


__all__  = ['Settings']

def is_not_none(key, value):
    if value is None:
        raise IOError('Setting for key %s should be specified, is now None.' %(key))

def has_value(values):
    if values is None: return
    values = [v.lower() for v in values]
    def check(key, value):
        if value is None: return
        if value.lower() not in values:
            raise IOError('Setting for key %s should be one of %s. Got %s' %(key, str(values), value))
    return check


def is_float(key, value):
    if value is None: return
    try:
        value = float(value)
    except ValueError:
        raise IOError('Setting for key %s should be of type float. Got %s.' %(key, str(value)))


def is_bool(key, value):
    if not isinstance(value, bool):
        raise IOError('Setting for key %s should be of type bool. Got %s.' %(key, str(value)))


def is_string(key, value):
    if value is None: return
    if not isinstance(value, str):
        raise IOError('Setting for key %s should be of type string. Got %s.' %(key, str(value)))


def is_list_strings(key, value):
    if value is None: return
    if ',' in value:
        value = value.split(',')
        for i,v in enumerate(value):
            if not isinstance(v,str):
                raise IOError('Setting for key %s should be a string or a list of strings. Element %i is %s.' %(key, i, v))
    else:
        if not isinstance(value, str):
            raise IOError('Setting for key %s should be a string or a list of strings. Got %s.' %(key, str(value)))


def is_nonexisting_file_name(key, value):
    if value is None: return
    if os.path.isfile(value):
        raise IOError('Setting for key %s should be non-existing file name, got %s which already exists.' %(key, value))


def is_existing_file_name(key, value):
    if value is None: return
    if not os.path.isfile(value):
        raise IOError('Setting for key %s should be existing file name, got %s which does not exist.' %(key, value))


key_checks = {
    'fn_yaff'               : [is_string, is_nonexisting_file_name],
    'fn_charmm22_prm'       : [is_string, is_nonexisting_file_name],
    'fn_charmm22_psf'       : [is_string, is_nonexisting_file_name],
    'fn_sys'                : [is_string, is_nonexisting_file_name],
    'plot_traj'             : [is_string, has_value(['None', 'Final', 'All'])],
    'xyz_traj'              : [is_bool],
    'fn_traj'               : [is_string],
    'log_level'             : [is_not_none, is_string, has_value(['silent','low','medium','high','highest'])],
    'log_file'              : [is_string, is_nonexisting_file_name],
    'program_mode'          : [is_not_none, has_value(['DeriveFF','MakeTrajectories','PlotTrajectories'])],
    'only_traj'             : [is_not_none, is_string],
    'ffatypes'              : [is_list_strings],
    'ei'                    : [is_string, is_existing_file_name],
    'ei_rcut'               : [is_float],
    'vdw'                   : [is_string, is_existing_file_name],
    'vdw_rcut'              : [is_float],
    'covres'                : [is_string, is_existing_file_name],
    'excl_bonds'            : [is_list_strings],
    'excl_bends'            : [is_list_strings],
    'excl_dihs'             : [is_list_strings],
    'excl_oopds'            : [is_list_strings],
    'do_hess_mass_weighting': [is_bool],
    'do_hess_negfreq_proj'  : [is_bool],
    'do_cross_svd'          : [is_bool],
    'cross_svd_rcond'       : [is_float],
    'pert_traj_tol'         : [is_float],
    'pert_traj_energy_noise': [is_float],
    'do_bonds'              : [is_bool],
    'do_bends'              : [is_bool],
    'do_dihedrals'          : [is_bool],
    'do_oops'               : [is_bool],
    'do_cross_ASS'          : [is_bool],
    'do_cross_ASA'          : [is_bool],
    'do_cross_DSS'          : [is_bool],
    'do_cross_DSD'          : [is_bool],
    'do_cross_DAA'          : [is_bool],
    'do_cross_DAD'          : [is_bool],
    'consistent_cross_rvs'  : [is_bool],
    'remove_dysfunctional_cross' : [is_bool],
    'bond_term'             : [is_not_none, is_string, has_value(['bondharm','bondfues','bondmm3'])],
    'bend_term'             : [is_not_none, is_string, has_value(['bendaharm','bendmm3'])],
    'do_squarebend'         : [is_bool],
    'do_bendclin'           : [is_bool],
    'do_sqoopdist_to_oopdist': [is_bool],
}


class Settings(object):
    'Class to control the behaviour of a Quickff run'
    def __init__(self, fn=None, **kwargs):
        '''
            **Keyword Arguments**

            fn
                file name of a config file from which settings can be read.
                Each line contains a single setting (except the lines starting
                with a #) and should have the following syntax

                    key: value

            kwargs
                each setting can also be parsed directly through the use of
                keyword arguments in the __init__ constructor with the syntax
                key=value. The settings parsed through the use of kwargs
                overwrite those in the config file.

        '''
        #first read general RC settings from .quickffrc file
        with path('quickff.data', 'quickffrc') as fn_default:
            self.read_config_file(fn_default)
        #if a config file is provided, read settings from this file and
        #overwrite the default RC settings
        if fn is not None:
            self.read_config_file(fn)
        #if settings are defined through keyword arguments to this init
        #constructor, read these settings and overwrite general RC as
        #wel as config file settings
        for key, value in kwargs.items():
            #don't impose keyword argument that hasn't been set
            if value is None: continue
            key = key.lstrip().rstrip()
            if key=='suffix':
                continue
            self.set(key, value)
        if 'suffix' in list(kwargs.keys()) and kwargs['suffix'] is not None:
            self._set_suffix(kwargs['suffix'])
        self.check()
        self._set_log()

        with log.section('SETT', 2, 'Initializing'):
            self.dump_log()

    def read_config_file(self, fn):
        with open(fn, 'r') as f:
            for iline, line in enumerate(f.readlines()):
                line = line.split('#')[0].lstrip().rstrip().rstrip('\n')
                if len(line)==0:
                    continue
                elif not ':' in line:
                    raise IOError('Line %i in %s does not contain a colon' %(iline, fn))
                else:
                    key, value = line.split(':')
                    key = key.lstrip().rstrip()
                    value = value.lstrip().rstrip()
                    value = value.lstrip().rstrip()
                    if value.lower()=='none':
                        value = None
                    elif value.lower()=='true':
                        value = True
                    elif value.lower()=='false':
                        value = False
                    else:
                        dtypes = [int, float]
                        found_dtype = False
                        for dtype in dtypes:
                            try:
                                value = dtype(value)
                                found_dtype = True
                                break
                            except ValueError:
                                try:
                                    value = parse_unit(value)
                                    found_dtype = True
                                    break
                                except ValueError:
                                    pass
                    self.set(key, value)

    def _set_suffix(self, suffix):
        for key, fn in self.__dict__.items():
            if fn is None or not key.startswith('fn_') or key=='fn_traj': continue
            prefix, extension = fn.split('.')
            self.__dict__[key] = '%s%s.%s' %(prefix, suffix, extension)

    def _set_log(self):
        log.set_level(self.log_level)
        f = self.log_file
        if f is not None and (isinstance(f, str) or isinstance(f, IOBase)):
            log.write_to_file(f)

    def set(self, key, value):
        if key not in list(key_checks.keys()):
            IOError('Key %s is not allowed in settings routine' %key)
        if isinstance(value, str) and value.lower()=='default':
            return
        self.__dict__[key] = value

    def check(self):
        for key, value in self.__dict__.items():
            for check_function in key_checks[key]:
                check_function(key,value)

    def dump_log(self):
        sorted_keys = sorted(self.__dict__.keys())
        with log.section('SETT', 3):
            for key in sorted_keys:
                value = str(self.__dict__[key])
                log.dump('%s  %s' %(key+' '*(30-len(key)), value))

    def dump_file(self, fn):
        sorted_keys = sorted(self.__dict__.keys())
        with open(fn, 'w') as f:
            for key in sorted_keys:
                value = str(self.__dict__[key])
                print('%s:   %s' %(key+' '*(30-len(key)), value), file=f)
