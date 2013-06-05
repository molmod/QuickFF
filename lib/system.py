#! /usr/bin/env python

from molmod.periodic import periodic as pt
from molmod.molecules import Molecule
from molmod.molecular_graphs import MolecularGraph
from molmod.io.psf import PSFFile
from molmod.io.chk import load_chk, dump_chk
from molmod.io.fchk import FCHKFile
from molmod.io.psf import PSFFile
from molmod.units import *
from molmod.ic import *

import numpy as np

from refdata import ReferenceData
from ic import IC
from atypes import assign_atypes

__all__ = ['System']

class System(object):
    def __init__(self, numbers, ffatypes, charges, ref, bonds, bends, diheds, nlist):
        '''
           **Arguments:**

           numbers
                A numpy array (N) with atomic numbers.

           ffatypes
                A numpy array (N) with force field atom types of a string
                specifying how to guess the atom types [high, medium or low].
 
           charges
                A numpy array (N) with atomic charges. The default charges
                are all zero.
           
           ref
                An instance of ReferenceData containing the reference data
                from which the force field will be estimated.
                     
           bonds
                A numpy array (B,2) with atom indexes (counting starts from
                zero) to define topological bond patterns. If no bonds are 
                given, they are estimated from the geometry in ref.coords.

           bends
                A numpy array (A,3) with atom indexes (counting starts from
                zero) to define topological bend patterns.

           diheds
                A numpy array (D,4) with atom indexes (counting starts from
                zero) to define topological dihedrall patterns.
           
           nlist
                A dictionnairy containing the neighbors for every every atom.
        '''
        self.numbers = numbers
        self.ffatypes = ffatypes
        self.charges = charges
        self.bonds = bonds
        self.bends = bends
        self.diheds = diheds
        self.ref = ref
        self.nlist = nlist
        self.ref.check()
        self.check_topology()
    
    def _get_natoms(self):
        return len(self.numbers)
    
    natoms = property(_get_natoms)
    
    @classmethod
    def from_files(cls, fns):
        '''
           A method to construct a System instance from input files. If
           the input files do not contain the topology, it is estimated
           from the geometry.
        
           **Arguments:**

           fns
                A list of file names. Files further on in the list will
                overwrite earlier files if there is overlap. txt files
                can only be used in hipart format to define charges.
        '''
        #initialise
        numbers = None
        charges = None
        ffatypes = None
        bonds = None
        bends = None
        diheds = None
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
                    elif key in ['ffatypes', 'atypes']:
                        ffatypes = value
                    elif key in ['bonds']:
                        bonds = value
                    elif key in ['bends', 'angles']:
                        bends = value
                    elif key in ['diheds', 'dihedrals']:
                        diheds = value
                    elif key in ['coords', 'coordinates', 'pos']:
                        ref.update(coords=value)
                    elif key in ['gradient', 'grad', 'gpos', 'forces']:
                        ref.update(grad=value)
                    elif key in ['hessian', 'hess']:
                        ref.update(hess=value)
                    else:
                        print 'WARNING: Skipped key %s in sample %s' %(key, fn)
            elif extension in ['fchk']:
                fchk = FCHKFile(fn)
                numbers = fchk.fields.get('Atomic numbers')
                ref.update(coords=fchk.fields.get('Current cartesian coordinates').reshape([len(numbers),3]))
                ref.update(grad=fchk.fields.get('Cartesian Gradient').reshape([len(numbers),3]))
                ref.update(hess=fchk.get_hessian().reshape([len(numbers),3, len(numbers),3]))
            elif extension in ['psf']:
                psf = PSFFile(fn)
                numbers = np.array(psf.numbers)
                ffatypes = np.array(psf.atom_types)
                charges = np.array(psf.charges)
                bonds = np.array(psf.bonds)
                bends = np.array(psf.bends)
                diheds = np.array(psf.dihedrals)
                nlist = psf.get_molecular_graph().neighbors
            elif extension in ['txt']:
                from hipart.io import load_atom_scalars
                charges = load_atom_scalars(fn)
        #Finalise topology and estimate atom types if necessary
        if charges is None:
            charges = np.zeros(len(numbers), float)
        return cls(numbers, ffatypes, charges, ref, bonds, bends, diheds, nlist)
    
    def check_topology(self):
        '''
           Check wether all topology information (bonds, bends, diheds
           and nlist) is present and complete if necessary.
        '''   
        if self.bonds is None:
            molecule = Molecule(self.numbers, coordinates=self.ref.coords)
            graph = MolecularGraph.from_geometry(molecule)
            psf = PSFFile()
            psf.add_molecule(molecule)
            self.bonds = np.array(psf.bonds)
            self.bends = np.array(psf.bends)
            self.diheds = np.array(psf.dihedrals)
            self.nlist = graph.neighbors
        if self.bends is None or self.diheds is None or self.nlist is None:
            graph = MolecularGraph(self.bonds, self.numbers)
            psf = PSFFile()
            psf.add_molecular_graph(graph)
            self.bends = np.array(psf.bends)
            self.dihedrals = np.array(psf.dihedrals)
            self.nlist = graph.neighbors
    
    def guess_atypes(self, atypes_level):
        '''
           A method to guess atom types. This will overwrite ffatypes
           that are already defined in the system.
        
           **Arguments:**
           
           atypes_level
                A string used for guessing atom types based on atomic
                number (low), local topology (medium) or atomic index
                in the molecule (high).
        '''
        self.ffatypes = assign_atypes(self, atypes_level)

    def determine_ics_from_topology(self):
        '''
            Method to generate IC instances corresponding to all ics
            present in the bonds, bends and dihedrals attributes.
        '''
        self.ics = {}
        number = {}
        def sort_ffatypes(ffatypes):
            if ffatypes[0]<=ffatypes[-1]: return ffatypes
            else: return ffatypes[::-1]
        #Find bonds
        for bond in self.bonds:
            name = 'bond/'+'.'.join(sort_ffatypes([self.ffatypes[at] for at in bond]))
            if name not in number.keys(): number[name] = 0
            else: number[name] += 1
            ic = IC(name+str(number[name]), bond, bond_length, qunit='A', kunit='kjmol/A**2')
            if name in self.ics.keys():
                self.ics[name].append(ic)
            else:
                self.ics[name] = [ic]
        #Find bends
        for bend in self.bends:
            name = 'angle/'+'.'.join(sort_ffatypes([self.ffatypes[at] for at in bend]))
            if name not in number.keys(): number[name] = 0
            else: number[name] += 1
            ic = IC(name+str(number[name]), bend, bend_angle, qunit='deg', kunit='kjmol/rad**2')
            if name in self.ics.keys():
                self.ics[name].append(ic)
            else:
                self.ics[name] = [ic]
        #Find dihedrals
        for dihed in self.diheds:
            name = 'dihed/'+'.'.join(sort_ffatypes([self.ffatypes[at] for at in dihed]))
            if name not in number.keys(): number[name] = 0
            else: number[name] += 1
            ic = IC(name+str(number[name]), dihed, dihed_angle, qunit='deg', kunit='kjmol')
            if name in self.ics.keys():
                self.ics[name].append(ic)
            else:
                self.ics[name] = [ic]

    def determine_ics_from_names(self, icnames):
        '''
            Method to generate IC instances corresponding to the names
            given as argument. For example, bond/Cph.Hph will search for
            all ICs corresponding to a bond between atoms of type Cph and 
            Hph.
            
            **Arguments**
            
            icnames
                A list containing the names of ics. Format of icnames is:
                kind/atype1.atype2...atypeN with kind=bond, bend or dihed
        '''
        self.ics = {}
        for icname in icnames:
            match = []
            atypes = icname.split('/')[1].split('.')
            ickind = icname.split('/')[0]
            if ickind in ['bond']:
                assert len(atypes)==2
                for bond in self.bonds:
                    if (self.ffatypes[bond]==atypes).all() or (self.ffatypes[bond]==atypes[::-1]).all():
                        match.append(IC(icname+str(len(match)), bond, bond_length, qunit='A', kunit='kjmol/A**2'))
            elif ickind in ['angle']:
                assert len(atypes)==3
                for bend in self.bends:
                    if (self.ffatypes[bend]==atypes).all() or (self.ffatypes[bend]==atypes[::-1]).all():
                        match.append(IC(icname+str(len(match)), bend, bend_angle, qunit='deg', kunit='kjmol/rad**2'))
            elif ickind in ['dihed']:
                assert len(atypes)==4
                for dihed in self.dihedrals:
                    if (self.ffatypes[dihed]==atypes).all() or (self.ffatypes[dihed]==atypes[::-1]).all():
                        match.append(IC(icname+str(len(match)), dihed, dihed_angle, qunit='deg', kunit='kjmol'))
            else:
                raise ValueError('Invalid icname, recieved %s' %icname)
            if len(match)==0:
                raise ValueError('Found no match for icname %s' %icname)
            self.ics[icname] = match

    def print_atom_info(self):
        print '    ---------------------------------'
        print '      index   symbol  ffatype  charge'
        print '    ---------------------------------'
        for index, (number, ffatype, charge) in enumerate(zip(self.numbers, self.ffatypes, self.charges)):
            print '        %3i      %3s %8s  % 6.3f' %(index, pt[number].symbol, ffatype, charge)
        print '    ---------------------------------'
        print
    
    def print_ic_info(self):
        print '    -----------------------------------------------------'
        print '                            icname     indices           '
        print '    -----------------------------------------------------'
        for icname, ics in self.ics.iteritems():
            for ic in ics:
                print '    %30s    ' %icname + ' '.join(['%3i' %index for index in ic.indexes])
        print '    -----------------------------------------------------'
        print

    def dump(self, fn):
        sample = {}
        sample['numbers']  = self.numbers
        sample['charges']  = self.charges
        sample['ffatypes'] = self.ffatypes
        sample['masses']   = np.array([pt[number].mass for number in self.numbers])
        sample['bonds']    = self.bonds
        sample['pos']      = self.ref.coords
        sample['grad']     = self.ref.grad
        sample['hess']     = self.ref.hess
        dump_chk(fn, sample)

    def dump_charges_yaff(self, fn, eirule, mode='w'):
        if eirule==-1:
            return
        elif eirule==0:
            scale1 = 1.0
            scale2 = 1.0
            scale3 = 1.0
        elif eirule==1:
            scale1 = 0.0
            scale2 = 1.0
            scale3 = 1.0
        elif eirule==2:
            scale1 = 0.0
            scale2 = 0.0
            scale3 = 1.0
        elif eirule==3:
            scale1 = 0.0
            scale2 = 0.0
            scale3 = 0.0
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
        print >> f, 'FIXQ:SCALE 1 %3.1f' %scale1
        print >> f, 'FIXQ:SCALE 2 %3.1f' %scale2
        print >> f, 'FIXQ:SCALE 3 %3.1f' %scale3
        print >> f, 'FIXQ:DIELECTRIC 1.0'
        print >> f, ''
        print >> f, '# Atomic parameters'
        print >> f, '# ----------------------------------------------------'
        print >> f, '# KEY        label  Q_0A              R_A'
        print >> f, '# ----------------------------------------------------'
        for atype, charge in zip(self.ffatypes, self.charges):
            print >> f, 'FIXQ:ATOM %8s % .10f  0.0000000000e+00' %(atype, charge)
        f.close()
