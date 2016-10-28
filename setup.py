#!/usr/bin/env python
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


from glob import glob
import os, sys
from distutils.core import setup
from distutils.command.install_data import install_data

class my_install_data(install_data):
    """Add a share_dir.txt file that points to the root for the shared files.
       It is otherwise impossible to figure out the location of these data
       files at runtime.
    """
    def run(self):
        # Do the normal install_data
        install_data.run(self)
        # Create the file share_dir.txt. It's exact content is only known
        # at installation time.
        dist = self.distribution
        libdir = dist.command_obj["install_lib"].install_dir
        for name in dist.packages:
            if '.' not in name:
                destination = os.path.join(libdir, name, "share_dir.txt")
                print "Creating %s" % destination
                if not self.dry_run:
                    f = file(destination, "w")
                    print >> f, self.install_dir
                    f.close()

def find_all_data_files(dn):
    result = []
    for root, dirs, files in os.walk(dn):
        if len(files) > 0:
            files = [os.path.join(root, fn) for fn in files]
            result.append(('share/quickff/' + root[6:], files))
    return result

setup(
    name='QuickFF',
    version='2.1.2',
    description='Python library to quickly derive force fields from ab initio training data.',
    author='Louis Vanduyfhuys',
    author_email='Louis.Vanduyfhuys@UGent.be',
    url='http://molmod.ugent.be/code/',
    package_dir = {'quickff': 'quickff'},
    packages=['quickff'],
    cmdclass = {'install_data': my_install_data},
    data_files=find_all_data_files('share'),
    scripts=['scripts/qff.py', 'scripts/qff-input-ei.py'],
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
