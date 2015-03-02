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

from molmod.periodic import periodic as pt
from molmod.molecules import Molecule
from molmod.molecular_graphs import MolecularGraph
from molmod.io.chk import load_chk, dump_chk
from molmod.io.fchk import FCHKFile
from molmod.io.psf import PSFFile
from molmod.units import deg, angstrom, kjmol
from molmod.ic import bond_length, bend_angle, dihed_angle, opbend_dist

import numpy as np

from quickff.refdata import ReferenceData
from quickff.ic import IC
from quickff.tools import find_opdist_patterns
from quickff.uff import get_uff_sigmas_epsilons

__all__ = ['System']

class System(object):
    '''
        A class for storing all system properties such as atom numbers, force
        field atom types, charges, vdW parameters, bonds, bends, dihedrals,
        out-of-plane patterns, neighborlists and the ab initio reference data.
    '''
    def __init__(self, numbers, ref, ffatypes=None, charges=None, radii=None,
                    epsilons=None, sigmas=None, bonds=None, bends=None, 
                    diheds=None, opdists=None, nlist=None):
        '''
           **Arguments:**

           numbers
                A numpy array (N) with atomic numbers.

           ref
                An instance of ReferenceData containing the reference data
                from which the force field will be estimated.

           **Optional Arguments**

           ffatypes
                A numpy array (N) with force field atom types.

           charges
                A numpy array (N) with atomic charges. By default, all charges
                are zero.

           radii
                A numpy array (N) with atomic radii determining the width of the
                charge distribution. By default, all radii are zero.

           epsilons
                A numpy array (N) with vdw epsilons values. By default, all
                epsilons are zero.
           
           sigmas
                A numpy array (N) with vdw sigma values. By default, all sigmas
                are zero.

           bonds
                A numpy array (B,2) with atom indexes (counting starts from
                zero) to define topological bond patterns. If no bonds are
                given, they are estimated from the geometry in ref.coords.

           bends
                A numpy array (A,3) with atom indexes (counting starts from
                zero) to define topological bend patterns. If no bends are
                given, they are derived from the bonds.

           diheds
                A numpy array (D,4) with atom indexes (counting starts from
                zero) to define topological dihedrall patterns. If no dihedrals
                are given, they are derived from the bonds.

           opdists
                A numpy array (O,4) with atom indexes (counting starts from
                zero) to define topological out of plane patterns. The first
                three atoms define the plane, while the last atom defines
                the center atom. If no opdists are given, they are derived
                from the bonds.

           nlist
                A dictionnairy containing the neighbors for every atom. If nlist
                is not given, it is derived from the bonds.
        '''
        self.numbers = numbers
        self.ffatypes = ffatypes
        self.charges = charges
        self.radii = radii
        self.epsilons = epsilons
        self.sigmas = sigmas
        self.bonds = bonds
        self.bends = bends
        self.diheds = diheds
        self.opdists = opdists
        self.nlist = nlist
        self.ref = ref
        self.ics = {}
        #Set attributes to zero if they are not defined
        if self.charges is None:
            self.charges = np.zeros(self.natoms, float)
        if self.radii is None:
            self.radii = np.zeros(self.natoms, float)
        if self.epsilons is None:
            self.epsilons = np.zeros(self.natoms, float)
        if self.sigmas is None:
            self.sigmas = np.zeros(self.natoms, float)
        #perfrom some checks
        self.ref._check(natoms=self.natoms)
        self._check_topology()
        if self.ffatypes is not None:
            self._average_charges_ffatypes()

    def _get_natoms(self):
        'Get the number of atoms in the system'
        return len(self.numbers)

    natoms = property(_get_natoms)

    @classmethod
    def from_files(cls, fns, ei_path=None, vdw_path=None):
        '''
           A method to construct a System instance from input files. If
           the input files do not contain the topology, it is estimated
           from the geometry.

           **Arguments:**

           fns
                A list of file names. Files further on in the list will
                overwrite earlier files if there is overlap. txt files
                can only be used in hipart format to define charges. HDF5
                files can only be used in Horton format to define charges.
                The set of charges that will be extracted from the HDF5
                file is dependant on the charge scheme defined in the
                kwargs.

           **Optional Arguments:**

           ei_path
                A string defining the path in the HDF5 file which contains a
                dataset `EI_PATH/charges` from which the atomic charges will be
                extracted.

           vdw_path
                A string defining the path in the HDF5 file which contains 2
                datasets, `VDW_PATH/epsilons` and `VDW_PATH/sigmas` from which the
                atomic vdW parameters will be extracted.
        '''
        #initialise
        numbers = None
        charges = None
        radii = None
        epsilons = None
        sigmas = None
        ffatypes = None
        bonds = None
        bends = None
        diheds = None
        opdists = None
        nlist = None
        ref = ReferenceData()
        #Read all input files
        for fn in fns:
            extension = fn.split('.')[-1]
            if extension in ['chk', 'mfs']:
                sample = load_chk(fn)
                for key, value in sample.iteritems():
                    if key in ['numbers']:
                        numbers = value
                    elif key in ['charges', 'ac', 'aq']:
                        charges = value
                    elif key in ['radii']:
                        radii = value
                    elif key in ['epsilons']:
                        epsilons = value
                    elif key in ['sigmas']:
                        sigmas = value
                    elif key in ['ffatypes', 'atypes']:
                        ffatypes = value
                    elif key in ['bonds']:
                        bonds = value
                    elif key in ['bends', 'angles']:
                        bends = value
                    elif key in ['diheds', 'dihedrals']:
                        diheds = value
                    elif key in ['opdists']:
                        opdists = value
                    elif key in ['coords', 'coordinates', 'pos']:
                        ref.update(coords=value.reshape([-1, 3]))
                    elif key in ['gradient', 'grad', 'gpos', 'forces']:
                        ref.update(grad=value.reshape([-1, 3]))
                    elif key in ['hessian', 'hess']:
                        natoms = int(np.sqrt(value.size/9))
                        ref.update(hess=value.reshape([natoms, 3, natoms, 3]))
            elif extension in ['fchk']:
                fchk = FCHKFile(fn)
                numbers = fchk.fields.get('Atomic numbers')
                ref.update(coords=fchk.fields.get('Current cartesian coordinates').reshape([len(numbers), 3]))
                ref.update(grad=fchk.fields.get('Cartesian Gradient').reshape([len(numbers), 3]))
                ref.update(hess=fchk.get_hessian().reshape([len(numbers), 3, len(numbers), 3]))
            elif extension in ['psf']:
                psf = PSFFile(fn)
                numbers = np.array(psf.numbers)
                ffatypes = np.array(psf.atom_types)
                charges = np.array(psf.charges)
                bonds = np.array(psf.bonds)
                bends = np.array(psf.bends)
                diheds = np.array(psf.dihedrals)
                nlist = psf.get_molecular_graph().neighbors
            elif extension in ['h5']:
                import h5py
                f = h5py.File(fn, 'r')
                if ei_path is None and vdw_path is None:
                    print '\nWARNING: a HDF5 file was given but both ei_path and'  +\
                          'vdw_path are None, HDF5 file is ignored.\n'
                if ei_path is not None:
                    charges = f['%s/charges' % ei_path][:]
                    if '%s/radii' %ei_path in f:
                        radii = f['%s/radii' %ei_path][:]
                if vdw_path is not None:
                    epsilons = f['%s/epsilons' % vdw_path][:]
                    sigmas = f['%s/sigmas' % vdw_path][:]
            else:
                raise IOError('Unsupported file format for %s' %fn)
        return cls(numbers, ref, ffatypes=ffatypes, charges=charges,
                   radii=radii, epsilons=epsilons, sigmas=sigmas, bonds=bonds, 
                   bends=bends, diheds=diheds, opdists=opdists, nlist=nlist)

    def read_uff_vdw(self):
        '''
            Read Lennart-Jones vdW parameters from the UFF force field.
        '''
        self.sigmas, self.epsilons = get_uff_sigmas_epsilons(self.numbers)

    def _check_topology(self):
        '''
           Check wether all topology information (bonds, bends, dihedrals,
           out-of-plane patterns and neighborlist) is present and complete if
           necessary.
        '''
        assert self.numbers is not None
        if self.bonds is None:
            molecule = Molecule(self.numbers, coordinates=self.ref.coords)
            graph = MolecularGraph.from_geometry(molecule)
            psf = PSFFile()
            psf.add_molecule(molecule)
            self.bonds = np.array(psf.bonds)
            self.bends = np.array(psf.bends)
            self.diheds = np.array(psf.dihedrals)
            self.opdists = find_opdist_patterns(graph)
            self.nlist = graph.neighbors
        else:
            graph = MolecularGraph(self.bonds, self.numbers)
            psf = PSFFile()
            psf.add_molecular_graph(graph)
            if self.bends is None:
                self.bends = np.array(psf.bends)
            if self.diheds is None:
                self.diheds = np.array(psf.dihedrals)
            if self.opdists is None:
                self.opdists = find_opdist_patterns(graph)
            if self.nlist is None:
                self.nlist = graph.neighbors

    def guess_ffatypes(self, level):
        '''
           A method to guess atom types. This will overwrite ffatypes
           that are already defined in the system.

           **Arguments:**

           level
                A string used for guessing atom types:

                    * low:     based on atomic number
                    * medium:  based on atomic number and number of neighbors
                    * high:    based on atomic number, number of neighbors and atomic number of neighbors
                    * highest: based on index in the molecule
        '''
        if self.ffatypes is not None:
            raise ValueError('Atom types are already defined in the system.')
        if level == 'low':
            atypes = np.array([pt[number].symbol for number in self.numbers])
        elif level == 'medium':
            atypes = []
            for index, number in enumerate(self.numbers):
                nind = self.nlist[index]
                sym = pt[self.numbers[index]].symbol.upper()
                atype = '%s%i' % (sym, len(nind))
                atypes.append(atype)
        elif level == 'high':
            atypes = []
            for index, number in enumerate(self.numbers):
                nind = self.nlist[index]
                nsym = sorted([
                    pt[self.numbers[neigh]].symbol.lower() for neigh in nind
                ])
                sym = pt[self.numbers[index]].symbol.upper()
                if len(nsym)==1:
                    atype = '%s1_%s' % (sym, nsym[0])
                elif len(nsym)==2:
                    atype = '%s2_%s%s' % (sym, nsym[0], nsym[1])
                else:
                    atype = '%s%i' % (sym, len(nind))
                    num_c = sum([1.0 for sym in nsym if sym == 'c'])
                    num_n = sum([1.0 for sym in nsym if sym == 'n'])
                    num_o = sum([1.0 for sym in nsym if sym == 'o'])
                    if num_c > 0: atype += '_c%i' % num_c
                    if num_n > 0: atype += '_n%i' % num_n
                    if num_o > 0: atype += '_o%i' % num_o
                atypes.append(atype)
        elif level == 'highest':
            atypes = np.array([
                '%s%i' % (pt[n].symbol, i) for i, n in enumerate(self.numbers)
            ])
        else:
            raise ValueError('Invalid level, recieved %s' % level)
        self.ffatypes = np.array(atypes)
        self._average_charges_ffatypes()

    def _average_charges_ffatypes(self):
        'Take the average of the atomic charges and radii over the force field atom types'
        charges = {}
        radii = {}
        for atype, charge, radius in zip(self.ffatypes, self.charges, self.radii):
            if atype not in charges:
                charges[atype] = [charge]
                radii[atype] = [radius]
            else:
                charges[atype].append(charge)
                radii[atype].append(radius)
        for i, atype in enumerate(self.ffatypes):
            self.charges[i] = sum(charges[atype])/len(charges[atype])
            self.radii[i] = sum(radii[atype])/len(radii[atype])

    def determine_ics_from_topology(self):
        '''
            Method to generate IC instances corresponding to all ics
            present in the bonds, bends, dihedrals and opdists attributes.
        '''
        self.ics = {}
        def sort_ffatypes(ffatypes, kind=''):
            if kind != 'opdist':
                if ffatypes[0] < ffatypes[-1]:
                    return ffatypes
                elif ffatypes[0]==ffatypes[-1]:
                    if len(ffatypes)==4:
                        if ffatypes[1]<ffatypes[2]:
                            return ffatypes
                        else:
                            return ffatypes[::-1]
                    else:
                        return ffatypes
                else:
                    return ffatypes[::-1]
            else:
                result = sorted(ffatypes[:3])
                result.append(ffatypes[3])
                return result
        #Find bonds
        number = {}
        for bond in self.bonds:
            name = 'bond/'+'.'.join(sort_ffatypes([self.ffatypes[at] for at in bond]))
            if name not in number.keys():
                number[name] = 0
            else:
                number[name] += 1
            ic = IC(
                name+str(number[name]), bond, bond_length,
                qunit='A', kunit='kjmol/A**2'
            )
            if name in self.ics.keys():
                self.ics[name].append(ic)
            else:
                self.ics[name] = [ic]
        #Find bends
        number = {}
        for bend in self.bends:
            name = 'angle/'+'.'.join(sort_ffatypes([self.ffatypes[at] for at in bend]))
            if name not in number.keys():
                number[name] = 0
            else:
                number[name] += 1
            ic = IC(
                name+str(number[name]), bend, bend_angle,
                qunit='deg', kunit='kjmol/rad**2'
            )
            if name in self.ics.keys():
                self.ics[name].append(ic)
            else:
                self.ics[name] = [ic]
        #Find dihedrals
        number = {}
        for dihed in self.diheds:
            if bend_angle(self.ref.coords[dihed[0:3]], deriv=0)[0] > 175*deg \
            or bend_angle(self.ref.coords[dihed[1:4]], deriv=0)[0] > 175*deg:
                continue
            name = 'dihed/'+'.'.join(sort_ffatypes([self.ffatypes[at] for at in dihed]))
            if name not in number.keys():
                number[name] = 0
            else:
                number[name] += 1
            ic = IC(
                name+str(number[name]), dihed, dihed_angle,
                qunit='deg', kunit='kjmol'
            )
            if name in self.ics.keys():
                self.ics[name].append(ic)
            else:
                self.ics[name] = [ic]
        #Find out-of-plane distances
        number = {}
        for opdist in self.opdists:
            if bend_angle(self.ref.coords[opdist[0:3]], deriv=0)[0] > 175*deg \
            or bend_angle(self.ref.coords[opdist[0:3]], deriv=0)[0] < 5*deg:
                continue
            name = 'opdist/'+'.'.join(sort_ffatypes([self.ffatypes[at] for at in opdist], kind='opdist'))
            if name not in number.keys(): number[name] = 0
            else: number[name] += 1
            ic = IC(
                name+str(number[name]), opdist, opbend_dist,
                qunit='A', kunit='kjmol/A**2'
            )
            if name in self.ics.keys():
                self.ics[name].append(ic)
            else:
                self.ics[name] = [ic]

    def print_atom_info(self):
        '''
            Print information about the atoms in the system to the screen
        '''
        print '    ----------------------------------------------------------------------------------------------'
        print '      index   symbol   ffatype       charge [e]     radius [A]    epsilon [kJ/mol]    sigma [A]   '
        print '    ----------------------------------------------------------------------------------------------'
        for index, (number, ffatype, charge, radius, epsilon, sigma) \
        in enumerate(zip(self.numbers, self.ffatypes, self.charges, self.radii, self.epsilons, self.sigmas)):
            index_fmt   = str(index)        + ' '*(5 - len(str(index)))
            symbol_fmt  = pt[number].symbol + ' '*(6 - len(pt[number].symbol))
            ffatype_fmt = ffatype           + ' '*(10 - len(ffatype))
            print '      %s   %s   %s      % 6.3f        %6.3f          %6.3f          %6.3f' % (
                index_fmt, symbol_fmt, ffatype_fmt, charge, radius/angstrom, epsilon/kjmol, sigma/angstrom
            )
        print '    ----------------------------------------------------------------------------------------------'

    def print_ic_info(self):
        '''
            Print information about the internal coordinate in the system
            to the screen
        '''
        print '    -----------------------------------------------------'
        print '                            icname     indices           '
        print '    -----------------------------------------------------'
        for icname, ics in self.ics.iteritems():
            for ic in ics:
                print '    %30s    ' % icname + ' '.join(
                    ['%3i' % index for index in ic.indexes]
                )
        print '    -----------------------------------------------------'

    def dump(self, fn):
        '''
            Dump system to a MolMod checkpoint file.
        '''
        sample = {}
        sample['numbers']   = self.numbers
        sample['charges']   = self.charges
        sample['radii']     = self.radii
        sample['epsilons']  = self.epsilons
        sample['sigmas']    = self.sigmas 
        sample['ffatypes']  = self.ffatypes
        sample['masses']    = np.array([pt[Z].mass for Z in self.numbers])
        sample['bonds']     = self.bonds
        sample['bends']     = self.bends
        sample['dihedrals'] = self.diheds
        sample['opdists']   = self.opdists
        sample['pos']       = self.ref.coords
        sample['grad']      = self.ref.grad
        sample['hess']      = self.ref.hess
        dump_chk(fn, sample)

    def dump_charges_yaff(self, fn, ei, mode='w'):
        '''
            Write or append charges to a file in Yaff format.
        '''
        f = open(fn, mode)
        print >> f, ''
        print >> f, '# Fixed charges'
        print >> f, '# ============='
        print >> f, ''
        print >> f, "# Mathematical form: q_A = q_0A + sum'_B p_BA"
        print >> f, '# where q0_A is the reference charge of atom A. It is mostly zero, sometimes a'
        print >> f, '# non-zero integer. The total charge of a fragment is the sum of all reference'
        print >> f, '# charges. The parameter p_BA is the charge transfered from B to A. Such charge'
        print >> f, '# transfers are only carried out over bonds in the FF topology.'
        print >> f, '# The charge on an atom is modeled as a Gaussian distribution. The spread on the'
        print >> f, '# Gaussian is called the radius R. When the radius is set to zero, point charges'
        print >> f, '# will be used instead of smeared charges.'
        print >> f, ''
        print >> f, 'FIXQ:UNIT Q0 e'
        print >> f, 'FIXQ:UNIT P e'
        print >> f, 'FIXQ:UNIT R angstrom'
        print >> f, 'FIXQ:SCALE 1 %3.1f' % ei.scales[0]
        print >> f, 'FIXQ:SCALE 2 %3.1f' % ei.scales[1]
        print >> f, 'FIXQ:SCALE 3 %3.1f' % ei.scales[2]
        print >> f, 'FIXQ:DIELECTRIC 1.0'
        print >> f, ''
        print >> f, '# Atomic parameters'
        print >> f, '# ----------------------------------------------------'
        print >> f, '# KEY        label  Q_0A              R_A'
        print >> f, '# ----------------------------------------------------'
        added = []
        for atype, q, radius in zip(self.ffatypes, self.charges, self.radii):
            if atype not in added:
                print >> f, 'FIXQ:ATOM %8s % 13.10f  %12.10f' % (atype, q, radius/angstrom)
                added.append(atype)
        f.close()

    def dump_vdw_yaff(self, fn, vdw, mode='w'):
        '''
            Write or append van der Waals parameters to a file in Yaff format.
        '''
        if vdw.pot.kind.startswith('MM3Buckingham'):
            pot_id = 'MM3'
        elif vdw.pot.kind.startswith('LennartJones'):
            pot_id = 'LJ'
        else:
            raise ValueError('VDWPart in model has unsupported pot_kind: %s' %vdw.kind)
        f = open(fn, mode)
        print >> f, '# van der Waals'
        print >> f, '#=============='
        print >> f, '# The following mathemetical form is supported:'
        print >> f, '#  - MM3:   EPSILON*(1.84e5*exp(-12*r/SIGMA)-2.25*(SIGMA/r)^6)'
        print >> f, '#  - LJ:    4.0*EPSILON*((SIGMA/r)^12 - (SIGMA/r)^6)'
        print >> f, '#  - PAULI: A*exp(-B*r)'
        print >> f, ''
        print >> f, '%s:UNIT SIGMA angstrom' %(pot_id)
        print >> f, '%s:UNIT EPSILON kjmol' %(pot_id)
        print >> f, '%s:SCALE 1 %4.2f' % (pot_id, vdw.scales[0])
        print >> f, '%s:SCALE 2 %4.2f' % (pot_id, vdw.scales[1])
        print >> f, '%s:SCALE 3 %4.2f' % (pot_id, vdw.scales[2])
        print >> f, ''
        print >> f, '# -------------------------------------------'
        print >> f, '# KEY    ffatype  SIGMA         EPSILON'
        print >> f, '# -------------------------------------------'
        added = []
        for atype, sigma, epsilon in zip(self.ffatypes, self.sigmas, self.epsilons):
            if atype not in added:
                print >> f, '%s:PARS %8s % .10f % .10f' %(pot_id, atype, sigma/angstrom, epsilon/kjmol)
                added.append(atype)
        f.close()
