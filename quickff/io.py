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
from molmod.units import angstrom, electronvolt, amu, kcalmol, kjmol, deg
from molmod.io.fchk import FCHKFile

from yaff.pes.ext import PairPotEI
from yaff.pes.ff import ForcePartPair

from quickff.log import log


__all__ = ['VASPRun', 'read_abinitio', 'make_yaff_ei', 'dump_charmm22_prm',
           'dump_charmm22_psf']


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
        energy = vasprun.fields['energies'][0]
        grad = -vasprun.fields['gradient'][0]
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


def _check_charmm22(valence):
    """Print warnings for all kinds of energy terms not supported by CHARMM22.

       **Arguments**

       valence
            Instance of ValenceFF, which defines the force field.
    """
    # First report all types of terms not supported by CHARMM22.
    skipped = set([])
    supported_kinds = set(['BONDHARM', 'BENDAHARM', 'TORSION'])
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
        ffatypes = master.basename.split('/')[1].split('.')
        K, q0 = valence.get_params(master.index)
        if K < 1e-6*kjmol/angstrom**2: continue
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
        ffatypes = master.basename.split('/')[1].split('.')
        K, q0 = valence.get_params(master.index)
        if K < 1e-6*kjmol: continue
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
    for master in valence.iter_masters(label='TORSION'):
        ffatypes = master.basename.split('/')[1].split('.')
        m, K, q0 = valence.get_params(master.index)
        if K < 1e-6*kjmol: continue
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
    for iatom in xrange(system.natom):
        ffatype = system.get_ffatype(iatom)
        if len(ffatype) > 4:
            log.warning('Atom type too long for CHARMM PSF file: {}'.format(ffatype))
        result.append('{:8d} A    1    MOL  {:4} {:4} {:10.6f} {:13.4f}           0'.format(
            iatom+1, ffatype, ffatype, system.charges[iatom], system.masses[iatom]/amu))
    return '\n'.join(result)


def _ics_to_charmm22_psf(valence, label, maxwidth, header, K_index, K_threshold):
    result = []
    width = maxwidth
    nic = 0
    for term in valence.iter_terms(label=label):
        K = valence.get_params(term.master)[K_index]
        if K < K_threshold: continue
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
    return _ics_to_charmm22_psf(valence, 'BONDHARM', 8, '{:8d} !NBOND: bonds', 0, 1e-6*kjmol/angstrom**2)


def _angles_to_charmm22_psf(valence):
    return _ics_to_charmm22_psf(valence, 'BENDAHARM', 9, '{:8d} !NTHETA: angles', 0, 1e-6*kjmol)


def _dihedrals_to_charmm22_psf(valence):
    return _ics_to_charmm22_psf(valence, 'TORSION', 8, '{:8d} !NPHI: dihedrals', 1, 1e-6*kjmol)


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
