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

from quickff.log import log
from molmod.units import parse_unit

__all__  = ['Settings']

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
        with log.section('SETT', 2, 'Initializing'):
            #first read general RC settings from .quickffrc file
            from context import context
            self.read_config_file(context.get_fn('quickffrc'))
            #if a config file is provided, read settings from this file and
            #overwrite the general RC settings
            if fn is not None:
                self.read_config_file(fn)
            #if settings are defined through keyword arguments to this init
            #constructor, read these settings and overwrite general RC as 
            #wel as config file settings
            for key, value in kwargs.iteritems():
                key = key.lstrip().rstrip()
                if key=='suffix':
                    continue
                self.set(key, value)
            if 'suffix' in kwargs.keys() and kwargs['suffix'] is not None:
                self._set_suffix(kwargs['suffix'])
            self._set_log()
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
        for key, fn in self.__dict__.iteritems():
            if not key.startswith('fn_') or key=='fn_traj': continue
            prefix, extension = fn.split('.')
            self.__dict__[key] = '%s%s.%s' %(prefix, suffix, extension)

    def _set_log(self):
        log.set_level(self.log_level)
        f = self.log_file
        if f is not None and (isinstance(f, str) or isinstance(f, file)):
            log.write_to_file(f)

    def set(self, key, value):
        if isinstance(value, str) and value.lower()=='default':
            return
        if value is None and key in self.__dict__.keys():
            return
        self.__dict__[key] = value

    def dump_log(self):
        sorted_keys = sorted(self.__dict__.keys())
        with log.section('SETT', 3):
            for key in sorted_keys:
                value = str(self.__dict__[key])
                log.dump('%s  %s' %(key+' '*(30-len(key)), value))
