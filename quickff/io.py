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

'''Readers ab initio vibrational calculations.
'''

import numpy as np
import xml.etree.ElementTree as ET

from molmod.periodic import periodic
from molmod.units import angstrom, electronvolt, amu
from molmod.io.fchk import FCHKFile

from yaff.pes.ext import PairPotEI
from yaff.pes.ff import ForcePartPair

from quickff.log import log

__all__ = ['VASPRun', 'read_abinitio', 'make_yaff_ei']


class VASPRun(object):
    #TODO Figure out the logic begind vasprun.xml to parse it in a more structured manner
    #TODO Test on more files
    '''
        Load information from a vasprun.xml file
    '''
    def __init__(self, filename, field_labels=[]):
        '''
            **Arguments**

            filename
                Filename of vasprun.xml

            **Optional Arguments**

            field_labels
                List of things we want to read. If an empty list is provided,
                only numbers, masses, initial positions and initial cell
                vectors are read.
        '''
        # Link between field labels and tags in the xml file
        tag_dictionary = {
            'rvecs_init': (".//structure[@name='initialpos']/crystal/varray[@name='basis']", angstrom),
            'pos_init'  : (".//structure[@name='initialpos']/varray[@name='positions']", 1.0),
            'gradient'  : (".//varray[@name='forces']", -electronvolt/angstrom),
            'hessian'   : (".//dynmat/varray[@name='hessian']", -electronvolt/angstrom**2/amu),
                }
        self.fields = {}
        self.tree = ET.parse(filename)
        self.root = self.tree.getroot()
        assert self.root.tag=='modeling', "Root tag is not modeling, this is not a standard vasprun.xml file"
        if not 'rvecs_init' in field_labels: field_labels.append('rvecs_init')
        if not 'pos_init' in field_labels: field_labels.append('pos_init')
        # Read atomic numbers
        self.fields['numbers'] = np.asarray([periodic[atom.find('c').text.strip()].number for atom in self.root.findall(".//*[@name='atoms']/set")[0]])
        # Read atomtypes
        atomtypes = np.asarray([int(atom.findall('c')[1].text) for atom in self.root.findall(".//*[@name='atoms']/set")[0]]) - 1
        # Read atomic masses for atomtype
        masses = np.asarray([float(atom.findall('c')[2].text) for atom in self.root.findall(".//*[@name='atomtypes']/set")[0]])*amu
        self.fields['masses'] = np.zeros(self.fields['numbers'].shape)
        for iatype in xrange(masses.shape[0]):
            self.fields['masses'][atomtypes==iatype] = masses[iatype]
        # Read SCF energies
        self.fields['energies'] = np.array([float(step.find('energy/i[@name="e_fr_energy"]').text)*electronvolt\
                 for step in self.root.findall('.//calculation')])
        # Read all requested arrays
        for label in field_labels:
            if not label in tag_dictionary.keys():
                raise NotImplementedError, "Failed to read %s from xml file" % label
            self.fields[label] = self.read_array(tag_dictionary[label][0], unit=tag_dictionary[label][1])
        # Convert fractional to Cartesian coordinates
        self.fields['pos_init'] = np.dot(self.fields['pos_init'], self.fields['rvecs_init'])
        # Hessian is mass-weighted, we want the pure second derivatives
        if 'hessian' in self.fields.keys():
            m3 = np.sqrt(np.array(sum([[m,m,m] for m in self.fields['masses']],[])))
            self.fields['hessian'] = m3.reshape((-1,1))*self.fields['hessian']*m3

    def read_array(self, tag, unit=1.0):
        result = []
        for match in self.root.findall(tag):
            result.append([])
            for line in match.findall('v'):
                result[-1].append([float(w) for w in line.text.split()])
        if len(result)==1: result = result[0]
        return np.asarray(result)*unit


def read_abinitio(fn):
    '''
        Wrapper to read all information from an ab initio calculation that
        QuickFF needs. Currently Gaussian .fchk and VASP .xml files are
        supported.
    '''
    extension = fn.split('.')[-1]
    if extension=='fchk':
        fchk = FCHKFile(fn)
        numbers = fchk.fields.get('Atomic numbers')
        energy = fchk.fields.get('Total Energy')
        coords = fchk.fields.get('Current cartesian coordinates').reshape([len(numbers), 3])
        grad = fchk.fields.get('Cartesian Gradient').reshape([len(numbers), 3])
        hess = fchk.get_hessian().reshape([len(numbers), 3, len(numbers), 3])
        masses = None
        rvecs = None
        pbc = [0,0,0]
    elif extension=='xml':
        vasprun = VASPRun(fn,field_labels=['hessian','gradient'])
        numbers = vasprun.fields['numbers']
        coords = vasprun.fields['pos_init']
        energy = vasprun.fields['energies'][0]
        grad = vasprun.fields['gradient'][0]
        hess = vasprun.fields['hessian'].reshape((len(numbers),3,len(numbers),3 ))
        masses = vasprun.fields['masses']
        rvecs = vasprun.fields['rvecs_init']
        pbc = [1,1,1]
    else: raise NotImplementedError
    return numbers, coords, energy, grad, hess, masses, rvecs, pbc


def make_yaff_ei(fn, charges, bcis=None, radii=None):
    assert charges is not None or bcis is not None, 'Either charges or bcis should be parsed'
    f = open(fn, 'w')
    print >> f, '#Fixed charges'
    print >> f, '#---------------'
    print >> f, ''
    print >> f, 'FIXQ:UNIT Q0 e'
    print >> f, 'FIXQ:UNIT P e'
    print >> f, 'FIXQ:UNIT R angstrom'
    print >> f, 'FIXQ:SCALE 1 1.0'
    print >> f, 'FIXQ:SCALE 2 1.0'
    print >> f, 'FIXQ:SCALE 3 1.0'
    print >> f, 'FIXQ:DIELECTRIC 1.0'
    print >> f, ''
    if charges is not None or radii is not None:
        print >> f, '# Atomic parameters'
        print >> f, '# ----------------------------------------------------'
        print >> f, '# KEY        label  Q_0A              R_A'
        print >> f, '# ----------------------------------------------------'
        if charges is not None:
            ffatypes = charges.keys()
        else:
            ffatypes = radii.keys()
        for ffatype in ffatypes:
            charge, radius = 0.0, 0.0
            if charges is not None: charge = charges[ffatype]
            if radii is not None: radius = radii[ffatype]
            print >> f, 'FIXQ:ATOM %8s % 13.10f  %12.10f' %(ffatype, charge, radius/angstrom)
    if bcis is not None:
        print >> f, '# Bond parameters'
        print >> f, '# ----------------------------------------------------'
        print >> f, '# KEY         label0   label1           P_AB          '
        print >> f, '# ----------------------------------------------------'
        for key, bci in bcis.iteritems():
            ffatype1, ffatype2 = key.split('.')
            print >> f, 'FIXQ:BOND  %8s  %8s  % 12.10f' %(ffatype1, ffatype2, bci)
    f.close()


def read_bci_constraints(fn):
    '''
        Read constraints for a charge to bci fit. The constraints should be
        written to a file in the following format:
        
            master0: slave00,slave01,slave02: sign
            master1: slave10,slave11,slave12: sign
        
        There should be a new line for each master and format is insensitive
        towards spaces (: and , serve as seperators). Lines starting with #
        are ignored (i.e. # is the comment identifier). Sign indicates if a
        sign switch should be introduced when mapping the slave bci to the
        master bci.
    '''
    constraints = {}
    with open(fn, 'r') as f:
        for line in f.readlines():
            if line.startswith('#'): continue
            master, suffix, sign = line.split(':')
            slaves = [(slave.strip(), float(sign)) for slave in suffix.split(',')]
            label = master.strip()
            if label in constraints.keys():
                for slave in slaves: constraints[label].append(slave)
            else:
                constraints[label] = slaves
    return constraints
