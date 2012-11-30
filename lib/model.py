#! /usr/bin/env python

from molmod.units import *
from molmod.ic import bond_length
from molmod.periodic import periodic as pt
import numpy as np, sys

from ic import *

__all__=['ZeroModel', 'HarmonicModel', 'CoulombModel', 'electrostatics']


class ZeroModel(object):
    def __init__(self, name='ZeroModel'):
        self.name = name
    
    def get_energy(self, coords):
        return 0.0
    
    def get_forces(self, coords):
        return np.zeros(coords.shape, float)
    
    def get_hessian(self, coords):
        return np.zeros([3*len(coords),3*len(coords)], float)


class HarmonicModel(object):
    def __init__(self, coords0, forces, hess, reference=0.0, name='Harmonic Model', ridge=1e-10):
        self.coords0 = coords0
        self.forces = forces
        self.hess = hess
        self.reference = reference
        self.ridge = ridge
        self.diagonalize_hess()
    
    def get_energy(self, coords):
        energy = self.reference
        dx = (coords - self.coords0).reshape((-1,1))
        energy += np.dot(self.forces.T, dx)[0]
        energy += 0.5*np.dot(dx.T, np.dot(self.hess, dx))[0,0]
        return energy
    
    def get_forces(self, coords):
        forces = self.forces
        dx = (coords - self.coords0).reshape((-1,1))
        return forces + np.dot(self.hess, dx)[0]
    
    def get_hessian(self, coords):
        return self.hess
    
    def diagonalize_hess(self):
        self.evals, self.evecs = np.linalg.eigh(self.hess)
        self.ievals = np.zeros(len(self.evals), float)
        for i, eigval in enumerate(self.evals):
            if abs(eigval)>self.ridge:
                self.ievals[i] = 1.0/eigval
        self.ihess = np.dot(self.evecs, np.dot(np.diag(self.ievals), self.evecs.T))        
    
    def get_constrained_hess(self, free_indices, spring=10.0*kjmol/angstrom**2):
        D = spring*np.identity(len(self.hess))
        for i in free_indices:
            D[i,i] = 0.0
        return self.hess + D
    
    def get_constrained_ihess(self, free_indices, spring=10.0*kjmol/angstrom**2):
        evals, evecs = np.linalg.eigh(self.get_constrained_hess(free_indices, spring=spring))
        ievals = np.zeros(len(evals), float)
        for i, eigval in enumerate(evals):
            if abs(eigval)>self.ridge:
                ievals[i] = 1.0/eigval
        return np.dot(evecs, np.dot(np.diag(ievals), evecs.T))
    
    def print_hess(self):
        from tools import global_translation, global_rotation, calc_angles
        red = "\033[31m%s\033[0m"
        green = "\033[32m%s\033[0m"
        print '====================================================================='
        print 'Printing hessian of %s' %self.name
        print '---------------------------------------------------------------------'
        VTx, VTy, VTz = global_translation(self.coords0)
        VRx, VRy, VRz = global_rotation(self.coords0)
        print ' Eigenvalues  |     Tx       Ty       Tz       Rx       Ry       Rz  '
        print '--------------+------------------------------------------------------'
        for i, eigval in enumerate(self.evals):
            ATx, ATy, ATz, ARx, ARy, ARz = calc_angles(self.evecs.T[i], [VTx, VTy, VTz, VRx, VRy, VRz])
            if abs(eigval)>self.ridge:
                print green %( '% .6e |  %7.3f  %7.3f  %7.3f  %7.3f  %7.3f  %7.3f' %(eigv, ATx/deg, ATy/deg, ATz/deg, ARx/deg, ARy/deg, ARz/deg) )
            else:
                print red %( '% .6e |  %7.3f  %7.3f  %7.3f  %7.3f  %7.3f  %7.3f' %(eigv, ATx/deg, ATy/deg, ATz/deg, ARx/deg, ARy/deg, ARz/deg) )
        print '====================================================================='



class CoulombModel(object):
    def __init__(self, coords, charges, name='Coulomb Model', exclude_pairs=[]):
        self.coords = coords
        self.charges = charges
        self.name = name
        self.exclude_pairs = exclude_pairs
        self.shift = self.get_energy(self.coords, shift=False)
    
    def get_energy(self, coords, shift=True):
        energy = 0.0
        if shift: energy -= self.shift
        for i, qi in enumerate(self.charges):
            for j, qj in enumerate(self.charges):
                if j<=i: continue
                if [i,j] in self.exclude_pairs or [j,i] in self.exclude_pairs: continue
                bond = IC([i, j], bond_length)
                energy += qi*qj/bond.value(coords)
        return energy



def electrostatics(sample, exclude_pairs=[], exclude_types=[]):
    qs = sample['charges']
    atypes = sample['ffatypes']
    forces = np.zeros(3*len(qs), float)
    hess = np.zeros([3*len(qs), 3*len(qs)], float)
    for i in xrange(len(qs)):
        if atypes[i] in exclude_types: continue
        for j in xrange(i):
            if atypes[j] in exclude_types: continue
            if [i,j] in exclude_pairs or [j,i] in exclude_pairs: continue
            bond = IC([i, j], bond_length)
            r = bond.value(sample['coordinates'])
            qgrad = bond.grad(sample['coordinates'])
            hess += qs[i]*qs[j]/(r**2)*(2.0/r*np.outer(qgrad, qgrad) - bond.hess(sample['coordinates']))
            forces += -qs[i]*qs[j]/(r**2)*qgrad
    return forces, hess
