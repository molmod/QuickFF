#!/usr/bin/env python
# -*- coding: utf-8 -*-
# QuickFF is a python library to quickly derive force fields from ab initio training data.
# Copyright (C) 2013 Louis Vanduyfhuys <Louis.Vanduyfhuys@UGent.be>,
# Center for Molecular Modeling (CMM), Ghent University, Ghent, Belgium;
# all rights reserved unless otherwise stated.
#
# This file is part of QuickFF.
#
# QuickFFis free software; you can redistribute it and/or
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


import glob
from distutils.core import setup

setup(
    name='QuickFF',
    version='0.1',
    description='Python library to quickly derive force fields from ab initio training data.',
    author='Louis Vanduyfhuys',
    author_email='Louis.Vanduyfhuys@UGent.be',
    url='http://molmod.ugent.be/code/',
    package_dir = {'quickff': 'lib'},
    packages=['quickff'],
    scripts=['scripts/qff-est.py', 'scripts/qff-traj.py'],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Topic :: Science/Engineering :: Molecular Science'
    ],
)
