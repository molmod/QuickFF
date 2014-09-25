# -*- coding: utf-8 -*-
#QuickFF is a code to quickly derive accurate force fields from ab initio input.
#Copyright (C) 2012 - 2014 Louis Vanduyfhuys <Louis.Vanduyfhuys@UGent.be>
#Steven Vandenbrande <Steven.Vandenbrande@UGent.be>,
#Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center for Molecular Modeling
#(CMM), Ghent University, Ghent, Belgium; all rights reserved unless otherwise
#stated.
#
#This file is part of QuickFF.
#
#QuickFF is free software; you can redistribute it and/or
#modify it under the terms of the GNU General Public License
#as published by the Free Software Foundation; either version 3
#of the License, or (at your option) any later version.
#
#QuickFF is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, see <http://www.gnu.org/licenses/>
#
#--

'''Defines the context in which QuickFF is used.

   This module controls global parameters that are purely technical, e.g. the
   location of example files needed for testing. It is certainly not meant to
   keep track of input parameters for a computation.

   This module contains a context object, an instance of the :class:`Context`
   class. For now, its functionality is rather limited. It tries to figure
   out the location of the share directory. It is assumed that the share
   directory is called ``share``. If the share directory does not exist,
   an error is raised.
'''


import os
from glob import glob


__all__ = ['context', 'Context']


class Context(object):
    '''Finds out where the share directory is located etc.

       The share directory contains examples directories with Gaussian FCHK,
       Horton HDF5, and text files with quickff reference output.
    '''
    def __init__(self):
        self.share_dir = os.getenv('QFFSHARE')
        if self.share_dir is None:
            fn_share_dir = os.path.join(os.path.dirname(__file__), 'share_dir.txt')
            if os.path.isfile(fn_share_dir):
                with open(fn_share_dir) as f:
                    self.share_dir = os.path.join(f.read().strip(), 'share/quickff')
        if self.share_dir is None:
            self.share_dir = './share'
        self.share_dir = os.path.abspath(self.share_dir)
        if not os.path.isdir(self.share_dir):
            raise IOError('Can not find the shared files. The directory %s does not exist.' % self.share_dir)

    def get_fn(self, filename):
        '''Return the full path to the given filename in the share directory.'''
        return os.path.join(self.share_dir, filename)

    def glob(self, pattern):
        '''Return all files in the share directory that match the given pattern.'''
        return glob(self.get_fn(pattern))


context = Context()
