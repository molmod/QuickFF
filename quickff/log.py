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

from __future__ import print_function, absolute_import, unicode_literals
from io import IOBase
import os, sys, datetime, getpass, atexit
import numpy, scipy, matplotlib

__all__ = ['version', 'log']

version = '2.2.4'

header = """
________________/\\\\\\_________/\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\__/\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\_________________
______________/\\\\\\\\/\\\\\\\\_____\\/\\\\\\///////////__\\/\\\\\///////////_________________
_____________/\\\\\\//\\////\\\\\\___\\/\\\\\\_____________\\/\\\\\\___________________________
_____________/\\\\\\______\\//\\\\\\__\\/\\\\\\\\\\\\\\\\\\\\\\_____\\/\\\\\\\\\\\\\\\\\\\\\\__________________
_____________\\//\\\\\\______/\\\\\\___\\/\\\\\\///////______\\/\\\\\\///////__________________
_______________\\///\\\\\\\\/\\\\\\\\/____\\/\\\\\\_____________\\/\\\\\\________________________
__________________\\////\\\\\\//______\\/\\\\\\_____________\\/\\\\\\_______________________
______________________\\///\\\\\\\\\\\\___\\/\\\\\\_____________\\/\\\\\\______________________
_________________________\\//////____\\///______________\\///______________________

                            Welcom to QuickFF
   a Python package to quickly derive force fields from ab initio input data

                                 Written by
     Louis Vanduyfhuys(1)*, Steven Vandenbrande(1) and Toon Verstraelen(1)

         (1) Center for Molecular Modeling, Ghent University Belgium.
                   * mailto: Louis.Vanduyfhuys@UGent.be
"""

footer = """
__/\\\\\\__________________________________________________________________/\\\\\\____
  \\ \\\\\\                                                                 \\ \\\\\\
   \\ \\\\\\     End of file. Thanks for using QuickFF! Come back soon!!     \\ \\\\\\
____\\///__________________________________________________________________\\///__
"""

def splitstring(string, length, separators=[' ','/','_']):
    result = []
    remainder = string
    while len(remainder)>length:
        i=length-1
        while remainder[i] not in separators and i>0:
            i -= 1
        result.append(remainder[:i+1])
        remainder = remainder[i+1:]
    if len(remainder)>0:
        result.append(remainder)
    return result


class Section(object):
    def __init__(self, logger, label, level, timer_description):
        self.logger = logger
        self.old_label = logger.label
        self.new_label = label
        self.old_level = self.logger.section_level
        self.new_level = level
        self.timer_description = timer_description

    def __enter__(self):
        if self.new_label!=self.old_label and self.logger.log_level>0:
            self.logger.add_blank_line = True
        self.logger.label = self.new_label
        self.logger.section_level = self.new_level
        if self.timer_description is not None:
            self.begin = datetime.datetime.now()
            self.end = None

    def __exit__(self, type, value, traceback):
        if self.new_label!=self.old_label and self.logger.log_level>0:
            self.logger.add_blank_line = True
        self.logger.label = self.old_label
        self.logger.section_level = self.old_level
        if self.timer_description is not None:
            self.end = datetime.datetime.now()
            for i, (description, time) in enumerate(self.logger.timetable):
                if description==self.timer_description:
                    self.logger.timetable[i][1] = time + self.end-self.begin
                    return
            self.logger.timetable.append([self.timer_description, self.end-self.begin])


class Logger(object):
    def __init__(self, level, _f=sys.stdout, max_label_length=8, line_length=80):
        self.set_level(level)
        self._f = _f
        self.mll = max_label_length
        self.ll = line_length
        self._active = False
        self.label = 'QFF'
        self.add_blank_line = False
        self.timetable = []

    def set_level(self, level):
        if isinstance(level, int):
            if level>=0 and level<=4:
                self.log_level = level
            else:
                raise ValueError('Integer level should be between 0 and 4 (boundaries included).')
        elif isinstance(level, str):
            allowed = ['silent', 'low', 'medium', 'high', 'highest']
            if level.lower() in allowed:
                self.log_level = allowed.index(level.lower())
            else:
                raise ValueError('String level should be silent, low, medium, high or highest.')
        self.section_level = None

    def write_to_file(self, f):
        if isinstance(f, str):
            self._f = open(f, 'w')
        elif isinstance(f, IOBase):
            self._f = f
        else:
            raise ValueError('File argument f should be a string representing a filename or a File instance')

    def section(self, label, level, timer=None):
        '''
            Construct a section instance for use in with statements to control
            section printing and timing.
        '''
        return Section(self, label, level, timer)

    def dump(self, message, new_line=True):
        if self.section_level<=self.log_level:
            if not self._active:
                self._active = True
                self.print_header()
            assert self.label is not None
            if new_line and self.add_blank_line:
                print('', file=self._f)
                self.add_blank_line = False
            line = ''
            for piece in splitstring(message, self.ll-self.mll):
                line += ' ' + self.label[:self.mll-2] + ' '
                line += ' '*(self.mll-2 - len(self.label[:self.mll-2]))
                line += piece
                line += '\n'
            line = line.rstrip('\n')
            print(line, file=self._f)

    def warning(self, message, new_line=True):
        '''
            Warnings are printed whenever log_level is higher than 0, i.e. not
            in silent mode.
        '''
        if self.log_level>0:
            if not self._active:
                self._active = True
                self.print_header()
            assert self.label is not None
            if new_line and self.add_blank_line:
                print('', file=self._f)
                self.add_blank_line = False
            line = ''
            for piece in splitstring('WARNING: '+message, self.ll-self.mll):
                line += ' ' + self.label[:self.mll-2] + ' '
                line += ' '*(self.mll-2 - len(self.label[:self.mll-2]))
                line += piece #.encode('utf-8')
                line += '\n'
            line = line.rstrip('\n')
            print(line, file=self._f)

    def print_header(self):
        if self.log_level>0:
            print(header, file=self._f)
            print('', file=self._f)
        mll = self.mll
        self.mll = 20
        with self.section('USER', 1): self.dump(getpass.getuser(), new_line=False)
        with self.section('MACHINE', 1): self.dump(' '.join(os.uname()), new_line=False)
        with self.section('TIME', 1): self.dump(datetime.datetime.now().isoformat().replace('T', ' '), new_line=False)
        with self.section('QUICKFF VERSION', 1): self.dump(version.replace('\n', ''), new_line=False)
        with self.section('PYTHON VERSION', 1): self.dump(sys.version.replace('\n', ''), new_line=False)
        with self.section('NUMPY VERSION', 1): self.dump(numpy.__version__, new_line=False)
        with self.section('SCIPY VERSION', 1): self.dump(scipy.__version__, new_line=False)
        with self.section('MATPLOTLIB VERSION', 1): self.dump(matplotlib.__version__, new_line=False)
        with self.section('CURRENT DIR', 1): self.dump(os.getcwd(), new_line=False)
        with self.section('COMMAND LINE', 1): self.dump(' '.join(sys.argv), new_line=False)
        self.mll = mll
        if self.log_level>0:
            print('', file=self._f)
            print('~'*80, file=self._f)

    def exit(self):
        if self._active:
            self.print_timetable()
            self.print_footer()
        self.close()

    def print_footer(self):
        if self.log_level>0:
            print(footer, file=self._f)

    def print_timetable(self):
        if self.log_level>0:
            print('~'*80, file=self._f)
            print('', file=self._f)
        with self.section('TIMING', 1):
            for label, time in self.timetable:
                line = '%30s  ' %(label+' '*(30-len(label)))
                line += str(time)
                self.dump(line)

    def close(self):
        if isinstance(self._f, IOBase):
            self._f.close()


log = Logger('medium')
atexit.register(log.exit)
