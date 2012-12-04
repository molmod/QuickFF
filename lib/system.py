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

__all__=['System']

class System(object):
    def __init__(self, name, fn_chk, fn_psf=None, eikind='Zero', eirule=-1, charges=None):
        print 'SYSTEM INIT : initializing system %s' %name
        self.name = name
        self.load_sample(fn_chk)
        self.get_topology(fn_psf)
        self.define_ei_model(eikind, eirule, charges)
        
        print 'SYSTEM PROJ : Projecting translational and rotional dofs out of the hessian'
        hess = self.sample['hessian'].reshape([3*self.Natoms, 3*self.Natoms])
        VTx, VTy, VTz = global_translation(self.sample['coordinates'])
        VRx, VRy, VRz = global_rotation(self.sample['coordinates'])
        U, S, Vt = np.linalg.svd( np.array([VTx, VTy, VTz, VRx, VRy, VRz]).transpose() )
        P = np.dot(U[:,6:], U[:,6:].T)
        hess_proj = np.dot(P, np.dot(hess, P)).reshape([self.Natoms, 3, self.Natoms, 3])
        self.totmodel = HarmonicModel(self.sample['coordinates'], self.sample['gradient'], hess_proj, name='Harmonic Total Energy')
    
    def load_sample(self, fn_chk):
        self.fn_chk = fn_chk
        if fn_chk.endswith('.chk'):
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
            self.sample['ffatypes'] = [pt[i].symbol for i in self.sample['numbers']]
        else:
            raise ValueError('Invalid file extension, recieved %s' %fn_chk)
    
    def get_topology(self, fn_psf=None, unit_cell=UnitCell(np.diag([50.0, 50.0, 50.0]),active=np.array([False, False, False])) ):
        self.fn_psf = fn_psf
        if fn_psf is not None:
            print 'SYSTEM TOPO : loading topology from %s' %fn_psf
            psf = PSFFile(fn_psf)
            self.sample['bonds'] = psf.bonds
            if len(psf.bends)>0: self.sample['bends'] = psf.bends
            if len(psf.dihedrals)>0: self.sample['dihedrals'] = psf.dihedrals
            self.sample['ffatypes'] = psf.atom_types
            self.neighbor_list = psf.get_molecular_graph().neighbors
            if 'charges' not in self.sample:
                print 'SYSTEM EI   : the charges are taken from %s and the radii are set to 0.01 A' %fn_psf
                self.sample['charges'] = psf.charges
                self.sample['radii'] = 0.01*angstrom*np.ones([self.Natoms], float)
        elif 'bonds' not in self.sample.keys():
            print "SYSTEM TOPO : creating topology from the geometry because 'bonds' was missing in %s" %self.fn_chk 
            molecule = Molecule(self.sample['numbers'], self.sample['coordinates'], unit_cell=unit_cell)
            graph = MolecularGraph.from_geometry(molecule)
            psf = PSFFile()
            psf.add_molecular_graph(graph)
            self.sample['bonds'] = psf.bonds
            if len(psf.bends)>0: self.sample['bends'] = psf.bends
            if len(psf.dihedrals)>0: self.sample['dihedrals'] = psf.dihedrals
            self.neighbor_list = graph.neighbors
        elif 'neighbor_list' not in self.sample.keys():
            print "SYSTEM TOPO : creating topology from the geometry because 'neighbors' was missing in %s" %self.fn_chk 
            molecule = Molecule(self.sample['numbers'], self.sample['coordinates'], unit_cell=unit_cell)
            graph = MolecularGraph.from_geometry(molecule)
            psf = PSFFile()
            psf.add_molecular_graph(graph)
            self.sample['bonds'] = psf.bonds
            if len(psf.bends)>0: self.sample['bends'] = psf.bends
            if len(psf.dihedrals)>0: self.sample['dihedrals'] = psf.dihedrals
            self.neighbor_list = graph.neighbors
        else:
            self.neighbor_list = self.sample['neighbor_list']
            print "SYSTEM TOPO : topology was loaded from %s" %self.fn_chk
    
    def define_ei_model(self, eikind, eirule, charges):
        self.eikind = eikind
        self.eirule = eirule
        print 'SYSTEM EI   : Constructing electrostatic model of kind %s and with exlusion rule = %i' %(eikind, eirule)
        if charges is not None:
            print 'SYSTEM EI   : the charges are set according to the command line input and the radii are set to 0.01 A'
            self.sample['charges'] = np.array([float(x) for x in charges.split(',')])
            self.sample['radii'] = 0.01*angstrom*np.ones([self.Natoms], float)
        elif charges not in self.sample.keys():
            raise ValueError("No charges present in sample, please specify them on the command line using the options '--charges'")
        if eikind!='Zero':
            exclude = []
            if eirule>0:
                for bond in self.sample['bonds']:
                    exclude.append([bond[0], bond[1]])
            if eirule>1 and 'bends' in self.sample.keys():
                for bend in self.sample['bends']:
                    exclude.append([bend[0], bend[2]])
            if eirule>2 and 'dihedrals' in self.sample.keys():
                for dihed in self.sample['dihedrals']:
                    exclude.append([dihed[0], dihed[3]])
        if eikind=='Zero':
            self.eimodel = ZeroModel()
        elif eikind=='Harmonic':
            forces_ei, hess_ei = electrostatics(self.sample, exclude_pairs=exclude)
            self.eimodel = HarmonicModel(self.sample['coordinates'], forces_ei, hess_ei, name='Harmonic Electrostatic Energy')
        elif eikind=='Coulomb':
            self.eimodel = CoulombModel(self.sample['coordinates'], self.sample['charges'], name='Coulomb Electrostatic Energy', exclude_pairs=exclude)
        else:
            raise ValueError('Invalid model kind in System.define_ei_model, received %s' %eikind)

    def dump_sample(self, fn):
        if not self.fn_chk==fn:
            from molmod.io.chk import dump_chk
            print 'SYSTEM DUMP : dumping system sample to MolMod checkpoint file %s' %fn
            dump_chk(fn, self.sample)

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
            print 'SYSTEM ICPAT: processing %s' %icname
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
    
    def test_ics(self, epsilon=1e-4, ntests=50, threshold=1e-5):
        print 'SYSTEM TEST : testing ic gradient and hessian'
        for icname, ics in self.ics.iteritems():
            for ic in ics:
                ic.test(self.sample['coordinates'], epsilon=epsilon, ntests=ntests, threshold=threshold)
        print 'SYSTEM TEST : testing non-bond pairs gradient and hessian'
        for i, itype in enumerate(self.sample['ffatypes']):
            for j, jtype in enumerate(self.sample['ffatypes'][:i]):
                if (i,j) in self.sample['bonds']: continue
                name = 'pair/%i(%s)-%i(%s)' %(i,itype,j,jtype)
                ic = IC([i,j], bond_length, name=name)
                ic.test(self.sample['coordinates'], epsilon=epsilon, ntests=ntests, threshold=threshold)

    def get_neighbors(self, indices, depth=1):
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
