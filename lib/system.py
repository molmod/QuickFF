#! /usr/bin/env python

from molmod.io.chk import load_chk
from molmod.io.fchk import FCHKFile
from molmod.io.psf import PSFFile
from molmod.ic import *
from molmod.periodic import periodic as pt
from molmod.molecular_graphs import MolecularGraph
from molmod.unit_cells import UnitCell
from molmod.molecules import Molecule
from molmod.units import *
from copy import deepcopy
import numpy as np, sys

from model import *
from ic import *
from tools import *
from atypes import assign_atypes

__all__=['System']

class System(object):
    def __init__(self, fn_chk, fn_psf=None, guess_atypes_level=None, charges=None):
        self.load_sample(fn_chk)
        if fn_psf is not None:
            self.load_topology(fn_psf)
        self.check_topology()
        if guess_atypes_level is not None:
            self.guess_atypes(guess_atypes_level)
        if charges is not None:
            self.load_charges(charges)
        self.print_info()
    
    def load_sample(self, fn_chk):
        self.fn_chk = fn_chk
        if fn_chk.endswith('.chk') or fn_chk.endswith('.mfs'):
            print 'SYSTEM LOAD : loading system sample from MolMod checkpoint file %s' %fn_chk
            self.sample = load_chk(fn_chk)
            self.Natoms = len(self.sample['numbers'])
        elif fn_chk.endswith('.fchk'):
            print 'SYSTEM LOAD : loading system sample from Gaussian fchk file %s' %fn_chk
            fchk = FCHKFile(fn_chk)
            self.sample = {}
            self.sample['title'] = 'Sample title'
            self.sample['numbers'] = fchk.fields.get('Atomic numbers')
            self.Natoms = len(self.sample['numbers'])
            self.sample['coordinates'] = fchk.fields.get('Current cartesian coordinates').reshape([self.Natoms, 3])
            self.sample['gradient'] = fchk.fields.get('Cartesian Gradient').reshape([self.Natoms, 3])
            self.sample['hessian']  = fchk.get_hessian().reshape([self.Natoms, 3, self.Natoms, 3])
        else:
            raise ValueError('Invalid file extension, recieved %s' %fn_chk)
        if 'colors' not in self.sample.keys():
            self.sample['colors'] = ['#7777CC' for i in xrange(self.Natoms)]
    
    def load_topology(self, fn_psf):
        print 'SYSTEM TOPOL: loading topology from %s' %fn_psf
        self.fn_psf = fn_psf
        psf = PSFFile(fn_psf)
        self.sample['bonds'] = psf.bonds
        if len(psf.bends)>0: self.sample['bends'] = psf.bends
        if len(psf.dihedrals)>0: self.sample['dihedrals'] = psf.dihedrals
        self.sample['ffatypes'] = psf.atom_types
        self.neighbor_list = psf.get_molecular_graph().neighbors
        if 'ac' not in self.sample:
            print 'SYSTEM EI   : the charges are taken from %s and the radii are set to 1e-4 A' %fn_psf
            self.sample['ac']   = psf.charges
            self.sample['radii'] = 1e-4*angstrom*np.ones([self.Natoms], float)

    def check_topology(self):
        """Check wether all topology information is present and complete topology
           if necessary
        """    
        if 'bonds' not in self.sample.keys():
            print "SYSTEM TOPOL: creating topology from the geometry because 'bonds' was missing"
            unit_cell = UnitCell(np.diag([50.0, 50.0, 50.0]),active=np.array([False, False, False]))
            molecule = Molecule(self.sample['numbers'], self.sample['coordinates'], unit_cell=unit_cell)
            graph = MolecularGraph.from_geometry(molecule)
            psf = PSFFile()
            psf.add_molecular_graph(graph)
            self.sample['bonds'] = psf.bonds
            if len(psf.bends)>0: self.sample['bends'] = psf.bends
            if len(psf.dihedrals)>0: self.sample['dihedrals'] = psf.dihedrals
            self.neighbor_list = graph.neighbors
        elif (len(self.sample['bonds'])>1 and 'bends' not in self.sample.keys()) \
          or (len(self.sample['bends'])>1 and 'dihedrals' not in self.sample.keys()):
            print "SYSTEM TOPOL: estimating bends, dihedrals and neighbor_list from bonds"
            graph = MolecularGraph(self.sample['bonds'], self.sample['numbers'])
            psf = PSFFile()
            psf.add_molecular_graph(graph)
            if len(psf.bends)>0: self.sample['bends'] = psf.bends
            if len(psf.dihedrals)>0: self.sample['dihedrals'] = psf.dihedrals
            self.neighbor_list = graph.neighbors
        elif not hasattr(self, 'neighbor_list'):
            print "SYSTEM TOPOL: neighbors estimated from bonds"
            graph = MolecularGraph(self.sample['bonds'], self.sample['numbers'])
            self.neighbor_list = graph.neighbors
    
    def guess_atypes(self, guess_atypes_level):
        print 'SYSTEM ATYPE: guessing atom types at %s level' %guess_atypes_level
        self.sample['ffatypes'] = assign_atypes(self.sample['bonds'], self.sample['numbers'], guess_atypes_level)
        
    def load_charges(self, charges):
        if charges.endswith('.txt'):
            print 'SYSTEM CHARG: the charges are set to the hipart charge file %s and the radii are set to 1e-4 A' %charges
            from hipart.io import load_atom_scalars
            self.sample['ac'] = load_atom_scalars(charges)
            self.sample['mc'] = sum(self.sample['ac'])
            self.sample['radii'] = 1e-4*angstrom*np.ones([self.Natoms], float)
        else:
            print 'SYSTEM CHARG: the charges are set according to the command line input and the radii are set to 1e-4 A'
            self.sample['ac'] = np.array([float(x) for x in charges.split(',')])
            self.sample['mc'] = sum(self.sample['ac'])
            self.sample['radii'] = 1e-4*angstrom*np.ones([self.Natoms], float)
    
    def define_models(self, eirule=-1, eikind=None):
        #Setting up total energy model
        print 'SYSTEM MODEL: Projecting translational and rotional dofs out of the hessian'
        hess = self.sample['hessian'].reshape([3*self.Natoms, 3*self.Natoms])
        VTx, VTy, VTz = global_translation(self.sample['coordinates'])
        VRx, VRy, VRz = global_rotation(self.sample['coordinates'])
        U, S, Vt = np.linalg.svd( np.array([VTx, VTy, VTz, VRx, VRy, VRz]).transpose() )
        P = np.dot(U[:,6:], U[:,6:].T)
        hess_proj = np.dot(P, np.dot(hess, P)).reshape([self.Natoms, 3, self.Natoms, 3])
        print 'SYSTEM MODEL: Setting up Harmonic model for total energy'
        self.totmodel = HarmonicModel(self.sample['coordinates'], self.sample['gradient'], hess_proj, name='Harmonic Total Energy')
        #Setting up electrostatic model
        if eirule==3 and not has_15_bonded(self):
            print 'SYSTEM MODEL: electrostatics switched off, because no 15 bonded pairs and exclude rule was %i' %eirule
            self.eirule = -1
            self.eikind = 'Zero'
        elif eirule==2 and 'dihedrals' not in self.sample.keys():
            print 'SYSTEM MODEL: electrostatics switched off, because no 14 bonded pairs and exclude rule was %i' %eirule
            self.eirule = -1
            self.eikind = 'Zero'
        elif eirule==1 and 'bends' not in self.sample.keys():
            print 'SYSTEM MODEL: electrostatics switched off, because no 13 bonded pairs and exclude rule was %i' %eirule
            self.eirule = -1
            self.eikind = 'Zero'
        else:
            if eikind is None:
                if eirule==-1:
                    eikind='Zero'
                else:
                    eikind='Harmonic'
            print 'SYSTEM MODEL: electrostatics switched on with kind %s and exlusion rule %i' %(eikind, eirule)
            self.eikind = eikind
            self.eirule = eirule
        if 'ac' not in self.sample.keys() and self.eikind!='Zero':
            raise ValueError("Charges need to be specified when electrostatics are switched on")
        if self.eikind=='Zero':
            self.eimodel = ZeroModel()
        else:
            exclude = []
            if self.eirule>0:
                for bond in self.sample['bonds']:
                    exclude.append([bond[0], bond[1]])
            if self.eirule>1 and 'bends' in self.sample.keys():
                for bend in self.sample['bends']:
                    exclude.append([bend[0], bend[2]])
            if self.eirule>2 and 'dihedrals' in self.sample.keys():
                for dihed in self.sample['dihedrals']:
                    exclude.append([dihed[0], dihed[3]])
            if self.eikind=='Harmonic':
                forces_ei, hess_ei = electrostatics(self.sample, exclude_pairs=exclude)
                self.eimodel = HarmonicModel(self.sample['coordinates'], forces_ei, hess_ei, name='Harmonic Electrostatic Energy')
            elif self.eikind=='Coulomb':
                self.eimodel = CoulombModel(self.sample['coordinates'], self.sample['ac'], name='Coulomb Electrostatic Energy', exclude_pairs=exclude)
            else:
                raise ValueError('Invalid electrostatic model kind, received %s' %self.eikind)

    def print_info(self):
        if 'ac' in self.sample.keys():
            print 'SYSTEM INFOR: The atom types, charges and charge radii are respectively'
            print
            for i, (t, q, r) in enumerate(zip(self.sample['ffatypes'], self.sample['ac'], self.sample['radii'])):
                print '                %3i  %8s % 6.3f %.4f' %(i, t, q, r)
            print
        else:
            print 'SYSTEM INFOR: The atom types are'
            print
            for i, t in enumerate(self.sample['ffatypes']):
                print '                %3i  %8s' %(i, t)
            print
    
    def dump_sample_qff(self, fn):
        if self.fn_chk==fn:
            raise ValueError('I refuse to overwrite input system file %s' %self.fn_chk)
        from molmod.io.chk import dump_chk
        print 'SYSTEM DUMPS: dumping system sample to QuickFF checkpoint file %s' %fn
        dump_chk(fn, self.sample)
    
    def dump_sample_yaff(self, fn):
        if self.fn_chk==fn:
            raise ValueError('I refuse to overwrite input system file %s' %self.fn_chk)
        print 'SYSTEM DUMPS: dumping system sample to Yaff checkpoint file %s' %fn
        yaff = {}
        yaff['bonds']    = self.sample['bonds']
        yaff['charges']  = self.sample['ac']
        yaff['ffatypes'] = self.sample['ffatypes']
        yaff['masses']   = np.array([pt[number].mass for number in self.sample['numbers']])
        yaff['numbers']  = self.sample['numbers']
        yaff['pos']      = self.sample['coordinates']
        from molmod.io.chk import dump_chk
        dump_chk(fn, yaff)

    def icnames_from_topology(self):
        print 'SYSTEM ICNAM: determining ic names from topolgy'
        icnames = []
        atypes = self.sample['ffatypes']
        for bond in self.sample['bonds']:
            name01 = 'bond/%s.%s' %(atypes[bond[0]], atypes[bond[1]])
            name10 = 'bond/%s.%s' %(atypes[bond[1]], atypes[bond[0]])
            if not (name01 in icnames or name10 in icnames):
                icnames.append(name01)
        if 'bends' in self.sample.keys():
            for bend in self.sample['bends']:
                name01 = 'angle/%s.%s.%s' %(atypes[bend[0]], atypes[bend[1]], atypes[bend[2]])
                name10 = 'angle/%s.%s.%s' %(atypes[bend[2]], atypes[bend[1]], atypes[bend[0]])
                if not (name01 in icnames or name10 in icnames):
                    icnames.append(name01)
        if 'dihedrals' in self.sample.keys():
            for dihedral in self.sample['dihedrals']:
                name01 = 'dihedral/%s.%s.%s.%s' %(atypes[dihedral[0]], atypes[dihedral[1]], atypes[dihedral[2]], atypes[dihedral[3]])
                name10 = 'dihedral/%s.%s.%s.%s' %(atypes[dihedral[3]], atypes[dihedral[2]], atypes[dihedral[1]], atypes[dihedral[0]])
                if not (name01 in icnames or name10 in icnames):
                    icnames.append(name01)
        return icnames
    
    def find_ic_patterns(self, icnames, units=None):
        if icnames==['all']: icnames = self.icnames_from_topology()
        print 'SYSTEM ICPAT: determining ic patterns'
        atypes = self.sample['ffatypes']
        self.icnames = icnames
        self.ics = {}
        if units is not None:
            assert isinstance(units, dict)
            self.units = units
        else:
            self.units = {}
        for icname in icnames:
            match = []
            values = []
            ictypes = icname.split('/')[1].split('.')
            ickind = icname.split('/')[0]
            if ickind in ['dist', 'bond']:
                assert len(ictypes)==2
                if units is None: self.units[icname] = {'q': 'A', 'k': 'kjmol/A**2'}
                for bond in self.sample['bonds']:
                    if (atypes[bond[0]]==ictypes[0] and atypes[bond[1]]==ictypes[1]) \
                    or (atypes[bond[0]]==ictypes[1] and atypes[bond[1]]==ictypes[0]):
                        match.append(IC(bond, bond_length, name=icname+str(len(match))))
            elif ickind in ['angle', 'bend']:
                assert len(ictypes)==3
                if units is None: self.units[icname] = {'q': 'deg', 'k': 'kjmol/rad**2'}
                for bend in self.sample['bends']:
                    if (atypes[bend[0]]==ictypes[0] and atypes[bend[1]]==ictypes[1] and atypes[bend[2]]==ictypes[2]) \
                    or (atypes[bend[0]]==ictypes[2] and atypes[bend[1]]==ictypes[1] and atypes[bend[2]]==ictypes[0]):
                        match.append(IC(bend, bend_angle, name=icname+str(len(match))))
            elif ickind in ['dihed', 'dihedral', 'torsion']:
                assert len(ictypes)==4
                if units is None: self.units[icname] = {'q': 'deg', 'k': 'kjmol'}
                for dihed in self.sample['dihedrals']:
                    if (atypes[dihed[0]]==ictypes[0] and atypes[dihed[1]]==ictypes[1] and atypes[dihed[2]]==ictypes[2] and atypes[dihed[3]]==ictypes[3]) \
                    or (atypes[dihed[0]]==ictypes[3] and atypes[dihed[1]]==ictypes[2] and atypes[dihed[2]]==ictypes[1] and atypes[dihed[3]]==ictypes[0]):
                        match.append(IC(dihed, dihed_angle, name=icname+str(len(match))))
            else:
                raise ValueError('Recieved invalid ic kind: %s' %ickind)
            if len(match)==0:
                print 'SYSTEM ICPAT: No match found for %s' %icname
            self.ics[icname] = match

    def get_neighbors(self, indices, depth=1):
        if depth==-1: return range(self.Natoms)
        neighbors = deepcopy(indices)
        edge = deepcopy(indices)
        current = 0
        while current<depth:
            new_edge = []
            for index in edge:
                for neighbor in self.neighbor_list[index]:
                    if neighbor not in neighbors:
                        neighbors.append(neighbor)
                        new_edge.append(neighbor)
            edge = new_edge
            current += 1
        return neighbors
