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

'''Readers ab initio vibrational calculations.
'''

import numpy as np
import xml.etree.ElementTree as ET

from molmod.periodic import periodic
from molmod.units import angstrom, electronvolt, amu, kcalmol, kjmol, deg
from molmod.io.fchk import FCHKFile

from yaff.pes.ext import PairPotEI
from yaff.pes.ff import ForcePartPair
from yaff.pes.parameters import *

from quickff.log import log


__all__ = ['VASPRun', 'read_abinitio', 'make_yaff_ei', 'dump_charmm22_prm',
           'dump_charmm22_psf', 'dump_yaff']


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
        for iatype in range(masses.shape[0]):
            self.fields['masses'][atomtypes==iatype] = masses[iatype]
        # Read SCF energies
        self.fields['energies'] = np.array([float(step.find('energy/i[@name="e_fr_energy"]').text)*electronvolt\
                 for step in self.root.findall('.//calculation')])
        # Read all requested arrays
        for label in field_labels:
            if not label in list(tag_dictionary.keys()):
                raise NotImplementedError("Failed to read %s from xml file" % label)
            self.fields[label] = self.read_array(tag_dictionary[label][0], unit=tag_dictionary[label][1])
        # Convert fractional to Cartesian coordinates
        self.fields['pos_init'] = np.dot(self.fields['pos_init'], self.fields['rvecs_init'])
        # Hessian is mass-weighted, we want the pure second derivatives
        if 'hessian' in list(self.fields.keys()):
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


def read_abinitio(fn, do_hess=True):
    '''
        Wrapper to read all information from an ab initio calculation that
        QuickFF needs. Currently Gaussian .fchk and VASP .xml files are
        supported.

        **Optional Arguments**

        do_hess
            Extract the hessian from the ab initio output. For qff-input-ei.py,
            it is interesting to be able to switch this off.
    '''
    extension = fn.split('.')[-1]
    if extension=='fchk':
        fchk = FCHKFile(fn)
        numbers = fchk.fields.get('Atomic numbers')
        energy = fchk.fields.get('Total Energy')
        coords = fchk.fields.get('Current cartesian coordinates').reshape([len(numbers), 3])
        grad = fchk.fields.get('Cartesian Gradient').reshape([len(numbers), 3])
        if do_hess:
            hess = fchk.get_hessian().reshape([len(numbers), 3, len(numbers), 3])
        else:
            hess = None
        masses = None
        rvecs = None
        pbc = [0,0,0]
    elif extension=='xml':
        vasprun = VASPRun(fn,field_labels=['hessian','gradient'])
        numbers = vasprun.fields['numbers']
        coords = vasprun.fields['pos_init']
        energy = vasprun.fields['energies'][-1]
        grad = -vasprun.fields['gradient'][-1]
        hess = vasprun.fields['hessian'].reshape((len(numbers),3,len(numbers),3 ))
        masses = vasprun.fields['masses']
        rvecs = vasprun.fields['rvecs_init']
        pbc = [1,1,1]
    else: raise NotImplementedError
    return numbers, coords, energy, grad, hess, masses, rvecs, pbc


def make_yaff_ei(fn, charges, bcis=None, radii=None, scales=[1,1,1]):
    assert charges is not None or bcis is not None, 'Either charges or bcis should be parsed'
    f = open(fn, 'w')
    print('#Fixed charges', file=f)
    print('#---------------', file=f)
    print('', file=f)
    print('FIXQ:UNIT Q0 e', file=f)
    print('FIXQ:UNIT P e', file=f)
    print('FIXQ:UNIT R angstrom', file=f)
    print('FIXQ:SCALE 1 %3.1f' %(scales[0]), file=f)
    print('FIXQ:SCALE 2 %3.1f' %(scales[1]), file=f)
    print('FIXQ:SCALE 3 %3.1f' %(scales[2]), file=f)
    print('FIXQ:DIELECTRIC 1.0', file=f)
    print('', file=f)
    if charges is not None or radii is not None:
        print('# Atomic parameters', file=f)
        print('# ----------------------------------------------------', file=f)
        print('# KEY        label  Q_0A              R_A', file=f)
        print('# ----------------------------------------------------', file=f)
        if charges is not None:
            ffatypes = list(charges.keys())
        else:
            ffatypes = list(radii.keys())
        for ffatype in ffatypes:
            charge, radius = 0.0, 0.0
            if charges is not None: charge = charges[ffatype]
            if radii is not None: radius = radii[ffatype]
            print('FIXQ:ATOM %8s % 13.10f  %12.10f' %(ffatype, charge, radius/angstrom), file=f)
    if bcis is not None:
        print('# Bond parameters', file=f)
        print('# ----------------------------------------------------', file=f)
        print('# KEY         label0   label1           P_AB          ', file=f)
        print('# ----------------------------------------------------', file=f)
        for key, bci in bcis.items():
            ffatype1, ffatype2 = key.split('.')
            print('FIXQ:BOND  %8s  %8s  % 12.10f' %(ffatype1, ffatype2, bci), file=f)
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
            if label in list(constraints.keys()):
                for slave in slaves: constraints[label].append(slave)
            else:
                constraints[label] = slaves
    return constraints


def _check_charmm22(valence):
    """Print warnings for all kinds of energy terms not supported by CHARMM22.

       **Arguments**

       valence
            Instance of ValenceFF, which defines the force field.
    """
    # First report all types of terms not supported by CHARMM22.
    skipped = set([])
    supported_kinds = set(['BONDHARM', 'BENDAHARM', 'TORSION', 'TORSCHEBY1', 'TORSCHEBY2', 'TORSCHEBY3', 'TORSCHEBY4', 'TORSCHEBY6'])
    for term in valence.iter_terms():
        kind = term.basename[:term.basename.find('/')].upper()
        if kind not in supported_kinds and kind not in skipped:
            log.warning('Not writing {} term to CHARMM22 parameter file.'.format(kind))
            skipped.add(kind)


# Template based on the PROTEINS-ZINC-NUCLEID-ACIDS force field file by MacKerell.
# (par_all27_prot_na.prm)
charmm22_prm_template = '''\

BONDS
!
!V(bond) = Kb(b - b0)**2
!
!Kb: kcal/mole/A**2
!b0: A
!
!at1  at2        Kb         b0
!
{}


ANGLES
!
!V(angle) = Ktheta(Theta - Theta0)**2
!
!V(Urey-Bradley) = Kub(S - S0)**2
!
!Ktheta: kcal/mole/rad**2
!Theta0: degrees
!Kub: kcal/mole/A**2 (Urey-Bradley)
!S0: A
!
!at1  at2  at3    Ktheta     Theta0     Kub     S0
!
{}


DIHEDRALS
!
!V(dihedral) = Kchi(1 + cos(n(chi) - delta))
!
!Kchi: kcal/mole
!n: multiplicity
!delta: degrees
!
!at1  at2  at3  at4       Kchi  n    delta
!
{}


IMPROPER
!
!V(improper) = Kpsi(psi - psi0)**2
!
!Kpsi: kcal/mole/rad**2
!psi0: degrees
!note that the second column of numbers (0) is ignored
!
!at1  at2  at3  at4   Kpsi        psi0
!

'''


def _bonds_to_charmm22_prm(valence):
    """Construct a BONDS section of a CHARMM22 parameter file.

       **Arguments**

       valence
            Instance of ValenceFF, which defines the force field.
    """
    result = []
    for master in valence.iter_masters(label='BONDHARM'):
        if valence.is_negligible(master.index): continue
        ffatypes = master.basename.split('/')[1].split('.')
        K, q0 = valence.get_params(master.index)
        result.append('{:4s} {:4s} {:9.3f} {:10.4f}'.format(
            ffatypes[0], ffatypes[1], 0.5*K/(kcalmol/angstrom**2), q0/angstrom
        ))
    return '\n'.join(result)


def _angles_to_charmm22_prm(valence):
    """Construct a ANGLES section of a CHARMM22 parameter file.

       **Arguments**

       valence
            Instance of ValenceFF, which defines the force field.
    """
    result = []
    for master in valence.iter_masters(label='BENDAHARM'):
        if valence.is_negligible(master.index): continue
        ffatypes = master.basename.split('/')[1].split('.')
        K, q0 = valence.get_params(master.index)
        result.append('{:4s} {:4s} {:4s} {:9.3f} {:10.4f}'.format(
            ffatypes[0], ffatypes[1], ffatypes[2], 0.5*K/kcalmol, q0/deg
        ))
    return '\n'.join(result)


def _dihedrals_to_charmm22_prm(valence):
    """Construct a DIHEDRALS section of a CHARMM22 parameter file.

       **Arguments**

       valence
            Instance of ValenceFF, which defines the force field.
    """
    result = []
    for master in valence.iter_masters(label='TORS'):
        if valence.is_negligible(master.index): continue
        ffatypes = master.basename.split('/')[1].split('.')
        kind = master.basename.split('/')[0].lower()
        if kind=='torsion':
            m, K, q0 = valence.get_params(master.index)
        elif kind.startswith('torscheby'):
            m = int(kind[-1])
            K = valence.get_params(master.index, only='fc')
            if valence.get_params(master.index, only='sign')<0:
                q0 = 0.0
            else:
                q0 = 180*deg/m
        result.append('{:4s} {:4s} {:4s} {:4s} {:10.4f} {:2.0f} {:8.2f}'.format(
            ffatypes[0], ffatypes[1], ffatypes[2], ffatypes[3], 0.5*K/kcalmol, m, q0/deg+180
        ))
    return '\n'.join(result)


def dump_charmm22_prm(valence, fn):
    """Dump supported parameters in a CHARMM22 parameter file.

       **Arguments**

       valence
            Instance of ValenceFF, which defines the force field.

       fn
            The filename to write to.
    """
    # First report all types of terms not supported by CHARMM22.
    _check_charmm22(valence)

    # Dump supported parameters in prm file.
    with open(fn, 'w') as f:
        f.write(charmm22_prm_template.format(
            _bonds_to_charmm22_prm(valence),
            _angles_to_charmm22_prm(valence),
            _dihedrals_to_charmm22_prm(valence)))


charmm22_psf_template = '''\
PSF

       1 !NTITLE
 REMARKS Generated with QuickFF

{}

{}

{}

{}

       0 !NIMPHI: impropers


       0 !NDON: donors


       0 !NACC: acceptors


       0 !NNB


'''

def _atoms_to_charmm22_psf(system):
    result = ['{:8d} !NATOM'.format(system.natom)]
    for iatom in range(system.natom):
        ffatype = system.get_ffatype(iatom)
        if len(ffatype) > 4:
            log.warning('Atom type too long for CHARMM PSF file: {}'.format(ffatype))
        if system.charges is None:
            charge = 0.0
        else:
            charge = system.charges[iatom]
        result.append('{:8d} A    1    MOL  {:4} {:4} {:10.6f} {:13.4f}           0'.format(
            iatom+1, ffatype, ffatype, charge, system.masses[iatom]/amu))
    return '\n'.join(result)


def _ics_to_charmm22_psf(valence, label, maxwidth, header):
    result = []
    width = maxwidth
    nic = 0
    for term in valence.iter_terms(label=label):
        if valence.is_negligible(term.master): continue
        # Extract iatom0, iatom1
        iatoms = term.get_atoms()
        # Predict what the width will be and fix if needed
        width += len(iatoms)
        if width >= maxwidth:
            width = 0
        nic += 1
        # Convert to string and add to result
        s = ''.join(['{:8d}'.format(iatom+1) for iatom in iatoms])
        if width == 0:
            result.append(s)
        else:
            result[-1] += s
    result.insert(0, header.format(nic))
    return '\n'.join(result)


def _bonds_to_charmm22_psf(valence):
    return _ics_to_charmm22_psf(valence, 'BONDHARM', 8, '{:8d} !NBOND: bonds')


def _angles_to_charmm22_psf(valence):
    return _ics_to_charmm22_psf(valence, 'BENDAHARM', 9, '{:8d} !NTHETA: angles')


def _dihedrals_to_charmm22_psf(valence):
    return _ics_to_charmm22_psf(valence, 'TORS', 8, '{:8d} !NPHI: dihedrals')


def dump_charmm22_psf(system, valence, fn):
    """Dump supported internal coordinates in a CHARMM psf file.

       **Arguments**

       system
            Instance of yaff.System class

       valence
            Instance of ValenceFF, which defines the force field.

       fn
            The filename to write to.
    """
    # First report all types of terms not supported by CHARMM22.
    _check_charmm22(valence)

    # Dump supported internal coordinates into PSF file.
    with open(fn, 'w') as f:
        f.write(charmm22_psf_template.format(
            _atoms_to_charmm22_psf(system),
            _bonds_to_charmm22_psf(valence),
            _angles_to_charmm22_psf(valence),
            _dihedrals_to_charmm22_psf(valence)))


def _bondharm_to_yaff(valence):
    'construct a BONDHARM section of a yaff parameter file'
    prefix = 'BONDHARM'
    units = ParameterDefinition('UNIT', lines=['K kjmol/A**2', 'R0 A'])
    pars = ParameterDefinition('PARS')
    for i, master in enumerate(valence.iter_masters(label=prefix)):
        if valence.is_negligible(master.index): continue
        ffatypes = master.basename.split('/')[1].split('.')
        K, q0 = valence.get_params(master.index)
        pars.lines.append('%8s  %8s  %.10e  %.10e' %(
            ffatypes[0], ffatypes[1], K/(kjmol/angstrom**2), q0/angstrom
        ))
    return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})

def _bondfues_to_yaff(valence):
    'construct a BONDHARM section of a yaff parameter file'
    prefix = 'BONDFUES'
    units = ParameterDefinition('UNIT', lines=['K kjmol/A**2', 'R0 A'])
    pars = ParameterDefinition('PARS')
    for i, master in enumerate(valence.iter_masters(label=prefix)):
        if valence.is_negligible(master.index): continue
        ffatypes = master.basename.split('/')[1].split('.')
        K, q0 = valence.get_params(master.index)
        pars.lines.append('%8s  %8s  %.10e  %.10e' %(
            ffatypes[0], ffatypes[1], K/(kjmol/angstrom**2), q0/angstrom
        ))
    return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})

def _bondmm3_to_yaff(valence):
    'construct a BONDHARM section of a yaff parameter file'
    prefix = 'MM3QUART'
    units = ParameterDefinition('UNIT', lines=['K kjmol/A**2', 'R0 A'])
    pars = ParameterDefinition('PARS')
    for i, master in enumerate(valence.iter_masters(label='BondMM3')):
        if valence.is_negligible(master.index): continue
        ffatypes = master.basename.split('/')[1].split('.')
        K, q0 = valence.get_params(master.index)
        pars.lines.append('%8s  %8s  %.10e  %.10e' %(
            ffatypes[0], ffatypes[1], K/(kjmol/angstrom**2), q0/angstrom
        ))
    return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})

def _bendaharm_to_yaff(valence):
    'construct a BENDAHARM section of a yaff parameter file'
    prefix = 'BENDAHARM'
    units = ParameterDefinition('UNIT', lines=['K kjmol/rad**2', 'THETA0 deg'])
    pars = ParameterDefinition('PARS')
    for i, master in enumerate(valence.iter_masters(label=prefix)):
        if valence.is_negligible(master.index): continue
        ffatypes = master.basename.split('/')[1].split('.')
        K, q0 = valence.get_params(master.index)
        pars.lines.append('%8s  %8s  %8s  %.10e  %.10e' %(
            ffatypes[0], ffatypes[1], ffatypes[2], K/kjmol, q0/deg
        ))
    return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})

def _bendmm3_to_yaff(valence):
    'construct a BENDAHARM section of a yaff parameter file'
    prefix = 'MM3BENDA'
    units = ParameterDefinition('UNIT', lines=['K kjmol/rad**2', 'THETA0 deg'])
    pars = ParameterDefinition('PARS')
    for i, master in enumerate(valence.iter_masters(label='BendMM3')):
        if valence.is_negligible(master.index): continue
        ffatypes = master.basename.split('/')[1].split('.')
        K, q0 = valence.get_params(master.index)
        pars.lines.append('%8s  %8s  %8s  %.10e  %.10e' %(
            ffatypes[0], ffatypes[1], ffatypes[2], K/kjmol, q0/deg
        ))
    return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})

def _bendcharm_to_yaff(valence):
    'construct a BENDCHARM section of a yaff parameter file'
    prefix = 'BENDCHARM'
    units = ParameterDefinition('UNIT', lines=['K kjmol', 'COS0 au'])
    pars = ParameterDefinition('PARS')
    for i, master in enumerate(valence.iter_masters(label=prefix)):
        if valence.is_negligible(master.index): continue
        ffatypes = master.basename.split('/')[1].split('.')
        K, q0 = valence.get_params(master.index)
        pars.lines.append('%8s  %8s  %8s  %.10e  %.10e' %(
            ffatypes[0], ffatypes[1], ffatypes[2], K/kjmol, q0
        ))
    return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})

def _bendcheby_to_yaff(valence, m):
    'construct a BENDCHEBY section of a yaff parameter file'
    prefix = 'BENDCOS'
    units = ParameterDefinition('UNIT', lines=['A kjmol', 'PHI0 deg'])
    pars = ParameterDefinition('PARS')
    for i, master in enumerate(valence.iter_masters(label='BendCheby%i'%m)):
        if valence.is_negligible(master.index): continue
        ffatypes = master.basename.split('/')[1].split('.')
        K, sign = valence.get_params(master.index)
        if sign==-1: q0=0.0
        else: q0=np.pi/m
        pars.lines.append('%8s  %8s  %8s  %1i %.10e  %.10e' %(
            ffatypes[0], ffatypes[1], ffatypes[2], m, K/kjmol, q0/deg
        ))
    return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})

def _torsions_to_yaff(valence):
    'construct a TORSION section of a yaff parameter file'
    prefix = 'TORSION'
    units = ParameterDefinition('UNIT', lines=['A kjmol', 'PHI0 deg'])
    pars = ParameterDefinition('PARS')
    for i, master in enumerate(valence.iter_masters(label=prefix)):
        if valence.is_negligible(master.index): continue
        ffatypes = master.basename.split('/')[1].split('.')
        m, K, q0 = valence.get_params(master.index)
        pars.lines.append('%8s  %8s  %8s  %8s  %1i %.10e  %.10e' %(
            ffatypes[0], ffatypes[1],  ffatypes[2], ffatypes[3], m,
            K/kjmol, q0/deg
        ))
    return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})

def _torscheby_to_yaff(valence, m):
    'construct a TORSION section of a yaff parameter file for chebychev dihedrals'
    prefix = 'TORSION'
    units = ParameterDefinition('UNIT', lines=['A kjmol', 'PHI0 deg'])
    pars = ParameterDefinition('PARS')
    for i, master in enumerate(valence.iter_masters(label='TorsCheby%i' %m)):
        if valence.is_negligible(master.index): continue
        ffatypes = master.basename.split('/')[1].split('.')
        K, sign = valence.get_params(master.index)
        if sign<0: q0 = 0.0
        else: q0 = np.pi/m
        pars.lines.append('%8s  %8s  %8s  %8s  %1i %.10e  %.10e' %(
            ffatypes[0], ffatypes[1],  ffatypes[2], ffatypes[3], m,
            K/kjmol, q0/deg
        ))
    return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})

def _torsc2harm_to_yaff(valence):
    'construct a TORSC2HARM section of a yaff parameter file'
    prefix = 'TORSC2HARM'
    units = ParameterDefinition('UNIT', lines=['A kjmol', 'COS0 au'])
    pars = ParameterDefinition('PARS')
    for i, master in enumerate(valence.iter_masters(label=prefix)):
        if valence.is_negligible(master.index): continue
        ffatypes = master.basename.split('/')[1].split('.')
        a = valence.get_params(master.index)
        K = 0.5*a[3]
        cos0 = np.arccos(np.sqrt(-0.5*a[1]/a[3]))
        pars.lines.append('%8s  %8s  %8s  %8s  %.10e  %.10e' %(
            ffatypes[0], ffatypes[1],  ffatypes[2], ffatypes[3],
            K/kjmol, cos0
        ))
    return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})

def _dihedharm_to_yaff(valence):
    'construct a DIHEDHARM section of a yaff parameter file'
    prefix = 'DIHEDHARM'
    units = ParameterDefinition('UNIT', lines=['K kjmol/rad**2', 'PHI0 deg'])
    pars = ParameterDefinition('PARS')
    for i, master in enumerate(valence.iter_masters(label=prefix)):
        if valence.is_negligible(master.index): continue
        ffatypes = master.basename.split('/')[1].split('.')
        K, phi0 = valence.get_params(master.index)
        pars.lines.append('%8s  %8s  %8s  %8s  %.10e  %.10e' %(
            ffatypes[0], ffatypes[1],  ffatypes[2], ffatypes[3],
            K/kjmol, phi0/deg
        ))
    return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})

def _opdists_to_yaff(valence):
    'construct a OOPDIST section of a yaff parameter file'
    prefix = 'OOPDIST'
    units = ParameterDefinition('UNIT', lines=['K kjmol/A**2', 'D0 A'])
    pars = ParameterDefinition('PARS')
    for i, master in enumerate(valence.iter_masters(label=prefix)):
        if valence.is_negligible(master.index): continue
        if 'sqoopdist' in master.basename.lower(): continue
        ffatypes = master.basename.split('/')[1].split('.')
        K, q0 = valence.get_params(master.index)
        pars.lines.append('%8s  %8s  %8s  %8s  %.10e  %.10e' %(
            ffatypes[0], ffatypes[1], ffatypes[2], ffatypes[3],
            K/(kjmol/angstrom**2), q0/angstrom
        ))
    return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})

def _sqopdists_to_yaff(valence):
    'construct a SQOOPDIST section of a yaff parameter file'
    prefix = 'SQOOPDIST'
    units = ParameterDefinition('UNIT', lines=['K kjmol/A**4', 'D0 A**2'])
    pars = ParameterDefinition('PARS')
    for i, master in enumerate(valence.iter_masters(label=prefix)):
        if valence.is_negligible(master.index): continue
        ffatypes = master.basename.split('/')[1].split('.')
        K, q0 = valence.get_params(master.index)
        pars.lines.append('%8s  %8s  %8s  %8s  %.10e  %.10e' %(
            ffatypes[0], ffatypes[1], ffatypes[2], ffatypes[3],
            K/(kjmol/angstrom**4), q0/angstrom**2
        ))
    return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})

def _cross_to_yaff(valence):
    'construct a CROSS section of a yaff parameter file'
    prefix = 'CROSS'
    units = ParameterDefinition(
        'UNIT',
        lines=[
            'KSS kjmol/angstrom**2', 'KBS0 kjmol/(angstrom*rad)',
            'KBS1 kjmol/(angstrom*rad)', 'R0 angstrom', 'R1 angstrom',
            'THETA0 deg'
        ]
    )
    done = []
    pars = ParameterDefinition('PARS')
    for i, master in enumerate(valence.iter_masters(label=prefix+'/')):
        if valence.is_negligible(master.index): continue
        prefix, ffatypes, suffix = master.basename.split('/')
        label = prefix+'/'+ffatypes+'/'
        if label in done: continue
        bb, b0a, b1a = None, None, None
        for j, other in enumerate(valence.iter_masters(label=label)):
            if 'bb' in other.basename:
                bb = valence.get_params(other.index)
            elif 'b0a' in other.basename:
                b0a = valence.get_params(other.index)
            elif 'b1a' in other.basename:
                b1a = valence.get_params(other.index)
            else:
                raise ValueError('Invalid Cross term %s' %other.basename)
        Kss, Kbs0, Kbs1 = 0.0, 0.0, 0.0
        r0, r1, theta0 = 0.0, 0.0, 0.0
        if bb is not None: Kss, r0, r1 = bb
        if b0a is not None: Kbs0, r0, theta0 = b0a
        if b1a is not None: Kbs1, r1, theta0 = b1a
        ffatypes = ffatypes.split('.')
        if ffatypes[::-1]!=ffatypes:
            pars.lines.append(
                '%8s  %8s  %8s  % .10e  % .10e  % .10e  %.10e  %.10e  %.10e' %(
                    ffatypes[0], ffatypes[1], ffatypes[2],
                    Kss/(kjmol/angstrom**2), Kbs0/(kjmol/angstrom),
                    Kbs1/(kjmol/angstrom), r0/angstrom, r1/angstrom, theta0/deg,
            ))
        #average parameters for symmetric terms
        else:
            pars.lines.append(
            '%8s  %8s  %8s  % .10e  % .10e  % .10e  %.10e  %.10e  %.10e' %(
                ffatypes[2], ffatypes[1], ffatypes[0],
                Kss/(kjmol/angstrom**2), 0.5*(Kbs0+Kbs1)/(kjmol/angstrom),
                0.5*(Kbs0+Kbs1)/(kjmol/angstrom), r1/angstrom, r0/angstrom,
                theta0/deg,
            ))
        done.append(label)
    return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})

def _crosscbend_to_yaff(valence):
    'construct a CROSSCBEND section of a yaff parameter file'
    prefix = 'CROSSCBEND'
    units = ParameterDefinition(
        'UNIT',
        lines=[
            'KSS kjmol/angstrom**2', 'KBS0 kjmol/angstrom',
            'KBS1 kjmol/angstrom', 'R0 angstrom', 'R1 angstrom',
            'COS0 au'
        ]
    )
    done = []
    pars = ParameterDefinition('PARS')
    for i, master in enumerate(valence.iter_masters(label=prefix+'/')):
        if valence.is_negligible(master.index): continue
        prefix, ffatypes, suffix = master.basename.split('/')
        label = prefix+'/'+ffatypes+'/'
        if label in done: continue
        bb, b0a, b1a = None, None, None
        for j, other in enumerate(valence.iter_masters(label=label)):
            if 'bb' in other.basename:
                bb = valence.get_params(other.index)
            elif 'b0a' in other.basename:
                b0a = valence.get_params(other.index)
            elif 'b1a' in other.basename:
                b1a = valence.get_params(other.index)
            else:
                raise ValueError('Invalid Cross term %s' %other.basename)
        Kss, Kbs0, Kbs1 = 0.0, 0.0, 0.0
        r0, r1, cos0 = 0.0, 0.0, 0.0
        if bb is not None: Kss, r0, r1 = bb
        if b0a is not None: Kbs0, r0, cos0 = b0a
        if b1a is not None: Kbs1, r1, cos0 = b1a
        ffatypes = ffatypes.split('.')
        if ffatypes[::-1]!=ffatypes:
            pars.lines.append(
                '%8s  %8s  %8s  % .10e  % .10e  % .10e  %.10e  %.10e  %.10e' %(
                    ffatypes[0], ffatypes[1], ffatypes[2],
                    Kss/(kjmol/angstrom**2), Kbs0/(kjmol/angstrom),
                    Kbs1/(kjmol/angstrom), r0/angstrom, r1/angstrom, cos0,
            ))
        else:
            pars.lines.append(
                '%8s  %8s  %8s  % .10e  % .10e  % .10e  %.10e  %.10e  %.10e' %(
                    ffatypes[0], ffatypes[1], ffatypes[2],
                    Kss/(kjmol/angstrom**2), 0.5*(Kbs0+Kbs1)/(kjmol/angstrom),
                    0.5*(Kbs0+Kbs1)/(kjmol/angstrom), r0/angstrom, r1/angstrom, cos0,
            ))
        done.append(label)
    return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})

def _crossbonddihed_to_yaff(valence, m):
    'construct a CROSS section of a yaff parameter file'
    prefix = 'CROSSBONDDIH%i' %m
    prefix.rstrip('1')
    units = ParameterDefinition(
        'UNIT',
        lines=[
            'KSS kjmol/angstrom**2', 'KSD0 kjmol/angstrom',
            'KSD1 kjmol/angstrom', 'KSD2 kjmol/angstrom',
            'R0 angstrom', 'R1 angstrom', 'R2 angstrom', 'CPSI0 au'
        ]
    )
    done = []
    pars = ParameterDefinition('PARS')
    for i, master in enumerate(valence.iter_masters(label=prefix+'/')):
        if valence.is_negligible(master.index): continue
        prefix, ffatypes, suffix = master.basename.split('/')
        label = prefix+'/'+ffatypes+'/'
        if label in done: continue
        bb, b0d, b1d, b2d = None, None, None, None
        for j, other in enumerate(valence.iter_masters(label=label)):
            if 'bb' in other.basename:
                bb = valence.get_params(other.index)
            elif 'b0d' in other.basename:
                b0d = valence.get_params(other.index)
            elif 'b1d' in other.basename:
                b1d = valence.get_params(other.index)
            elif 'b2d' in other.basename:
                b2d = valence.get_params(other.index)
            else:
                raise ValueError('Invalid Cross term %s' %other.basename)
        Kss, Ksd0, Ksd1, Ksd2 = 0.0, 0.0, 0.0, 0.0
        r0, r1, r2, cpsi0 = 0.0, 0.0, 0.0, 0.0
        if bb is not None: Kss, r0, r2 = bb
        if b0d is not None: Ksd0, r0, cpsi0 = b0d
        if b1d is not None: Ksd1, r1, cpsi0 = b1d
        if b2d is not None: Ksd2, r2, cpsi0 = b2d
        ffatypes = ffatypes.split('.')
        if ffatypes[::-1]!=ffatypes:
            pars.lines.append(
                '%8s  %8s  %8s  %8s  % .10e  % .10e  % .10e  %.10e  %.10e  %.10e  %.10e  %.10e' %(
                    ffatypes[0], ffatypes[1], ffatypes[2], ffatypes[3],
                    Kss/(kjmol/angstrom**2), Ksd0/(kjmol/angstrom),
                    Ksd1/(kjmol/angstrom), Ksd2/(kjmol/angstrom),
                    r0/angstrom, r1/angstrom, r2/angstrom, cpsi0,
            ))
        else:
            pars.lines.append(
            '%8s  %8s  %8s  %8s  % .10e  % .10e  % .10e  %.10e  %.10e  %.10e  %.10e  %.10e' %(
                ffatypes[0], ffatypes[1], ffatypes[2], ffatypes[3],
                Kss/(kjmol/angstrom**2), 0.5*(Ksd0+Ksd2)/(kjmol/angstrom),
                Ksd1/(kjmol/angstrom), 0.5*(Ksd0+Ksd2)/(kjmol/angstrom),
                r0/angstrom, r1/angstrom, r2/angstrom, cpsi0,
            ))
        done.append(label)
    return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})

def _crossbenddihed_to_yaff(valence, m):
    'construct a CROSS section of a yaff parameter file'
    prefix = 'CROSSBENDDIH%i' %m
    prefix.rstrip('1')
    units = ParameterDefinition(
        'UNIT',
        lines=[
            'KBB kjmol/rad**2', 'KBD0 kjmol', 'KBD1 kjmol',
            'THETA0 deg', 'THETA1 deg', 'CPSI0 au',
        ]
    )
    done = []
    pars = ParameterDefinition('PARS')
    for i, master in enumerate(valence.iter_masters(label=prefix+'/')):
        if valence.is_negligible(master.index): continue
        prefix, ffatypes, suffix = master.basename.split('/')
        label = prefix+'/'+ffatypes+'/'
        if label in done: continue
        aa, a0d, a1d = None, None, None
        for j, other in enumerate(valence.iter_masters(label=label)):
            if 'aa' in other.basename:
                aa = valence.get_params(other.index)
            elif 'a0d' in other.basename:
                a0d = valence.get_params(other.index)
            elif 'a1d' in other.basename:
                a1d = valence.get_params(other.index)
            else:
                raise ValueError('Invalid Cross term %s' %other.basename)
        Kaa, Kad0, Kad1 = 0.0, 0.0, 0.0
        theta0, theta1, cpsi0 = 0.0, 0.0, 0.0
        if aa is not None: Kaa, theta0, theta1 = aa
        if a0d is not None: Kad0, theta0, cpsi0 = a0d
        if a1d is not None: Kad1, theta1, cpsi0 = a1d
        ffatypes = ffatypes.split('.')
        if ffatypes[::-1]!=ffatypes:
            pars.lines.append(
                '%8s  %8s  %8s  %8s  % .10e  % .10e  % .10e  %.10e  %.10e  %.10e' %(
                    ffatypes[0], ffatypes[1], ffatypes[2], ffatypes[3],
                    Kaa/kjmol, Kad0/kjmol, Kad1/kjmol,
                    theta0/deg, theta1/deg, cpsi0,
            ))
        else:
            pars.lines.append(
                '%8s  %8s  %8s  %8s  % .10e  % .10e  % .10e  %.10e  %.10e  %.10e' %(
                    ffatypes[0], ffatypes[1], ffatypes[2], ffatypes[3],
                    Kaa/kjmol, 0.5*(Kad0+Kad1)/kjmol, 0.5*(Kad0+Kad1)/kjmol,
                    theta0/deg, theta1/deg, cpsi0,
            ))
        done.append(label)
    return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})

def _crosscbenddihed_to_yaff(valence):
    'construct a CROSS section of a yaff parameter file'
    prefix = 'CROSSCBENDDIH'
    units = ParameterDefinition(
        'UNIT',
        lines=[
            'KBB kjmol', 'KBD0 kjmol', 'KBD1 kjmol',
            'CTHETA0 au', 'CTHETA1 au', 'CPSI0 au'
        ]
    )
    done = []
    pars = ParameterDefinition('PARS')
    for i, master in enumerate(valence.iter_masters(label=prefix+'/')):
        if valence.is_negligible(master.index): continue
        prefix, ffatypes, suffix = master.basename.split('/')
        label = prefix+'/'+ffatypes+'/'
        if label in done: continue
        aa, a0d, a1d = None, None, None
        for j, other in enumerate(valence.iter_masters(label=label)):
            if 'aa' in other.basename:
                aa = valence.get_params(other.index)
            elif 'a0d' in other.basename:
                a0d = valence.get_params(other.index)
            elif 'a1d' in other.basename:
                a1d = valence.get_params(other.index)
            else:
                raise ValueError('Invalid Cross term %s' %other.basename)
        Kaa, Kad0, Kad1 = 0.0, 0.0, 0.0
        ctheta0, ctheta1, cpsi0 = 0.0, 0.0, 0.0
        if aa is not None: Kaa, ctheta0, ctheta1 = aa
        if a0d is not None: Kad0, ctheta0, cpsi0 = a0d
        if a1d is not None: Kad1, ctheta1, cpsi0 = a1d
        ffatypes = ffatypes.split('.')
        if ffatypes[::-1]!=ffatypes:
            pars.lines.append(
                '%8s  %8s  %8s  %8s  % .10e  % .10e  % .10e  %.10e  %.10e  %.10e' %(
                    ffatypes[0], ffatypes[1], ffatypes[2], ffatypes[3],
                    Kaa/kjmol, Kad0/kjmol, Kad1/kjmol,
                    ctheta0, ctheta1, cpsi0,
            ))
        else:
            pars.lines.append(
                '%8s  %8s  %8s  %8s  % .10e  % .10e  % .10e  %.10e  %.10e  %.10e' %(
                    ffatypes[0], ffatypes[1], ffatypes[2], ffatypes[3],
                    Kaa/kjmol, 0.5*(Kad0+Kad1)/kjmol, 0.5*(Kad0+Kad1)/kjmol,
                    ctheta0, ctheta1, cpsi0,
            ))
        done.append(label)
    return ParameterSection(prefix, definitions={'UNIT': units, 'PARS': pars})


def dump_yaff(valence, fn):
    sections = [
        _bondharm_to_yaff(valence), _bondfues_to_yaff(valence),
        _bondmm3_to_yaff(valence),
        _bendaharm_to_yaff(valence), _bendmm3_to_yaff(valence),
        _bendcharm_to_yaff(valence), _bendcheby_to_yaff(valence, 1),
        _bendcheby_to_yaff(valence, 4),
        _torsions_to_yaff(valence), _torsc2harm_to_yaff(valence),
        _torscheby_to_yaff(valence,1), _torscheby_to_yaff(valence,2),
        _torscheby_to_yaff(valence,3), _torscheby_to_yaff(valence,4),
        _torscheby_to_yaff(valence,6), _dihedharm_to_yaff(valence),
        _opdists_to_yaff(valence), _sqopdists_to_yaff(valence),
        _cross_to_yaff(valence), _crosscbend_to_yaff(valence),
        _crossbonddihed_to_yaff(valence,1), _crossbonddihed_to_yaff(valence,2),
        _crossbonddihed_to_yaff(valence,3), _crossbonddihed_to_yaff(valence,4),
        _crossbonddihed_to_yaff(valence,6),
        _crossbenddihed_to_yaff(valence,1), _crossbenddihed_to_yaff(valence,2),
        _crossbenddihed_to_yaff(valence,3), _crossbenddihed_to_yaff(valence,4),
        _crossbenddihed_to_yaff(valence,6),
        _crosscbenddihed_to_yaff(valence),
    ]
    f = open(fn, 'w')
    for section in sections:
        if len(section['PARS'].lines)==0: continue
        print('# %s' %section.prefix, file=f)
        print('#-%s' %('-'*len(section.prefix)), file=f)
        for line in section['UNIT'].lines:
            print('%s:UNIT  %s' %(section.prefix, line), file=f)
        print('', file=f)
        for line in section['PARS'].lines:
            print('%s:PARS  %s' %(section.prefix, line), file=f)
        print('', file=f)
        print('', file=f)
    f.close()
