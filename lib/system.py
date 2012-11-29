#! /usr/bin/env python

from molmod.io.chk import load_chk
from molmod.io.psf import PSFFile
from molmod.io.xyz import XYZWriter
from molmod.units import *
from molmod.ic import *
from molmod.periodic import periodic as pt
from pytool.math import statistics
from pytool.fit import fitpar
from pytool.physics import global_translation, global_rotation, calc_angles
from copy import deepcopy
import numpy as np, sys

sys.path.append('/home/louis/Documents/Doctoraat/Code/mil53/model/metals/covalent/est/lib')
from model import *
from ic import *
from evaluators import *

__all__=['System', 'FFTable']

class System(object):
    def __init__(self, name, fn_chk, fn_psf=None, ei_exclude=0):
        print '%s: initializing' %name
        self.name = name
        self.fn_chk = fn_chk
        self.fn_psf = fn_psf
        self.ei_exclude = ei_exclude
        self.sample = load_chk(fn_chk)
        if self.fn_psf is not None:
            print '%s: Getting topology info' %self.name
            psf = PSFFile(self.fn_psf)
            self.sample['bends'] = psf.bends
            self.sample['dihedrals'] = psf.dihedrals
            self.sample['neighbors'] = psf.get_molecular_graph().neighbors
        print '%s: Projecting translational and rotional dofs out of the hessian' %self.name
        self.Nat = len(self.sample['coordinates'])
        hess = self.sample['hessian'].reshape([3*self.Nat, 3*self.Nat])
        VTx, VTy, VTz = global_translation(self.sample['coordinates'])
        VRx, VRy, VRz = global_rotation(self.sample['coordinates'])
        U, S, Vt = np.linalg.svd( np.array([VTx, VTy, VTz, VRx, VRy, VRz]).transpose() )
        oc = U[:,6:]
        hess_proj = np.dot(oc.transpose(), np.dot(hess, oc))
        hess_proj = np.dot(oc, np.dot(hess_proj, oc.transpose()))
        self.totharm = HarmonicModel(self.sample['coordinates'], np.zeros(3*self.Nat, float), hess_proj, name='Harmonic Total Energy')
        print '%s: Constructing electrostatic model' %self.name
        exclude = []
        if self.ei_exclude>0:
            for bond in self.sample['bonds']:
                exclude.append([bond[0], bond[1]])
        if self.ei_exclude>1:
            for bend in self.sample['bends']:
                exclude.append([bend[0], bend[2]])
        if self.ei_exclude>2:
            for dihed in self.sample['dihedrals']:
                exclude.append([dihed[0], dihed[3]])
        forces_ei, hess_ei = electrostatics(self.sample, exclude_pairs=exclude)
        self.eiharm = HarmonicModel(self.sample['coordinates'], forces_ei, hess_ei, name='Harmonic Electrostatic Energy')
        self.eicoul = CoulombModel(self.sample['coordinates'], self.sample['charges'], name='Coulomb Electrostatic Energy')

    def find_ic_patterns(self, icnames):
        print '%s: finding ic patterns' %self.name
        atypes = self.sample['ffatypes']
        self.ics = {}
        self.qunit = {}
        self.kunit = {}
        for icname, (kunit, qunit) in icnames.iteritems():
            match = []
            values = []
            if '%s' in icname: icname = icname %( (self.name.split('/')[0].upper(),)*icname.count('%') )
            self.qunit[icname] = qunit
            self.kunit[icname] = kunit
            ictypes = icname.split('/')[1].split('.')
            ickind = icname.split('/')[0]
            if ickind=='bond':
                assert len(ictypes)==2
                for bond in self.sample['bonds']:
                    if (atypes[bond[0]]==ictypes[0] and atypes[bond[1]]==ictypes[1]) or (atypes[bond[0]]==ictypes[1] and atypes[bond[1]]==ictypes[0]):
                        match.append(IC(bond, bond_length, name=icname+str(len(match))))
            elif ickind=='angle':
                assert len(ictypes)==3
                for bend in self.sample['bends']:
                    if (atypes[bend[0]]==ictypes[0] and atypes[bend[1]]==ictypes[1] and atypes[bend[2]]==ictypes[2]) \
                    or (atypes[bend[0]]==ictypes[2] and atypes[bend[1]]==ictypes[1] and atypes[bend[2]]==ictypes[0]):
                        match.append(IC(bend, bend_angle, name=icname+str(len(match))))
            elif ickind=='dihed':
                assert len(ictypes)==4
                for dihed in self.sample['dihedrals']:
                    if (atypes[dihed[0]]==ictypes[0] and atypes[dihed[1]]==ictypes[1] and atypes[dihed[2]]==ictypes[2] and atypes[dihed[3]]==ictypes[3]) \
                    or (atypes[dihed[0]]==ictypes[3] and atypes[dihed[1]]==ictypes[2] and atypes[dihed[2]]==ictypes[1] and atypes[dihed[3]]==ictypes[0]):
                        match.append(IC(dihed, dihed_angle, name=icname+str(len(match))))
            else:
                raise ValueError('Recieved invalid ic kind: %s' %ickind)
            if len(match)==0:
                print '  No match found for %s' %(icname)
            self.ics[icname] = match

    def test_ics(self, epsilon=1e-4, ntests=50, threshold=1e-5):
        print '%s: testing ic gradient and hessian implementations' %self.name
        print ' - testing ics'
        for icname, ics in self.ics.iteritems():
            for ic in ics:
                ic.test(self.sample['coordinates'], epsilon=epsilon, ntests=ntests, threshold=threshold)
        print ' - testing non-bond pairs'
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
                for neighbor in self.sample['neighbors'][index]:
                    if neighbor not in neighbors:
                        neighbors.append(neighbor)
                        new_edge.append(neighbor)
            edge = new_edge
            current += 1
        return neighbors
    
    def geometry_perturbation(self, icname, coupling=None, free_depth=None, spring=10.0*kjmol/angstrom**2):
        """
            Get the perturbation on the geometry resulting from perturbing along
            the ics belonging to icname. Coupling is a numpy array with the 
            coefficients of the linear combination describing the coupling.
            If free_depth is None, the hessian remains unbiased for calculating 
            the perturbed geometry. If free_depth is e.g. equal to 3, all atoms 
            that are separated by more then 3 atoms from the atoms in the ic 
            under consideration are fixed with a spring (with strength defined
            in spring).
            
            The following strings are also supported for coupling:
                ``symmetric``   A symmetric coupling 
                                i.e. coupling = np.ones(N)/np.sqrt(N)
                                with N the number of ics for the given icname
        """
        ics = self.ics[icname]
        if coupling=='symmetric':
            coupling = np.ones(len(ics), float)/len(ics)
        elif isinstance(coupling, np.ndarray):
            assert len(ics)==len(coupling)
        elif coupling is not None:
            raise NotImplementedError('Unsupported coupling: %s' %(str(coupling)))
        qgrads = []
        for i, ic in enumerate(ics):
            qgrads.append(ic.grad(self.sample['coordinates']))
        if coupling is not None:
            qgrad_coupled = np.zeros(3*self.Nat, float)
            for i, qgrad in enumerate(qgrads):
                qgrad_coupled += qgrad*coupling[i]
            for i in xrange(len(ics)):
                qgrads[i] = qgrad_coupled
        vs = []
        for qgrad in qgrads:
            if free_depth is None:
                ihess = self.totharm.ihess
            else:
                indices = [i for i in xrange(self.Nat) if np.linalg.norm(qgrad.reshape([self.Nat, 3])[i])>1e-3]
                free_indices = self.get_neighbors(indices, depth=free_depth)
                ihess = self.totharm.get_constrained_ihess(free_indices, spring=spring)
            v = np.dot(ihess, qgrad)
            kl = np.dot(qgrad.T, np.dot(ihess, qgrad))
            vs.append(v.reshape((-1,3))/kl)
        return vs
    
    def analyze_trajectory(self, v, evaluators=[], fn_xyz=None, amplitude=0.25*angstrom, steps=101):
        symbols = np.array([pt[n].symbol for n in self.sample['numbers']])
        if fn_xyz is not None: traj = XYZWriter(file(fn_xyz, 'w'), symbols)
        values = [[] for i in xrange(len(evaluators))]
        for n in xrange(steps):
            pre = amplitude*np.sin(-np.pi/2+np.pi*n/(steps-1))
            coords = self.sample['coordinates'] + pre*v
            for i, evaluator in enumerate(evaluators):
                values[i].append(evaluator(self, coords))
            if fn_xyz is not None: traj.dump('frame %i' %n, coords)
        if fn_xyz is not None: del(traj)
        return values



class Array(object):
    def __init__(self, data, unit=1.0):
        if isinstance(data, np.ndarray):
            self.data = data
        else:
            self.data = np.array(data)
        self.unit = unit
        self.mean, self.std, self.num = statistics(self.data)
    
    def __len__(self):
        return len(self.data)
    
    def html(self, fmt):
        if self.num==0: return ''
        if self.std==None: self.std = 0.0
        return fmt %(self.mean/parse_unit(self.unit), self.std/parse_unit(self.unit), self.num)



class FFTable(object):    
    def __init__(self, system, coupling=None, free_depth=None, spring=10.0*kjmol/angstrom**2):
        self.coupling = coupling
        self.free_depth = free_depth
        self.spring = spring
        self.estimate(system)
    
    def estimate(self, system):
        self.name = system.name
        print '%s: estimating force constants [coupling=%s, free_depth=%s and spring=%.3e au]' %(self.name, str(self.coupling), str(self.free_depth), self.spring)
        self.k     = {}
        self.k_cov = {}
        self.q     = {}
        self.q_cov = {}
        for icname, ics in system.ics.iteritems():
            k_arr     = []
            k_cov_arr = []
            q_arr     = []
            q_cov_arr = []
            vs = system.geometry_perturbation(icname, self.coupling, free_depth=self.free_depth, spring=self.spring)
            if   icname.split('/')[0]=='bond' : amplitude_factor = 0.1
            elif icname.split('/')[0]=='angle': amplitude_factor = 0.01
            elif icname.split('/')[0]=='dihed': amplitude_factor = 0.01
            else: raise ValueError('Invalid ictype: %s' %icname.split('/')[0])
            for iv, v in enumerate(vs):
                amplitude = amplitude_factor*ics[iv].value(system.sample['coordinates'])
                values = system.analyze_trajectory(
                    v, evaluators=[ic_evaluator(icname, iv), energy_evaluator('totharm'), energy_evaluator('eiharm')], amplitude=amplitude
                )
                qs = np.array(values[0])
                tot = np.array(values[1])
                ei = np.array(values[2])
                if self.coupling is not None:
                    tot /= len(ics)
                    ei /= len(ics)
                pars = fitpar(qs, tot, rcond=1e-6)
                k = 2*pars[0]
                q = -pars[1]/k
                pars_cov = fitpar(qs, tot-ei, rcond=1e-6)
                k_cov = 2*pars_cov[0]
                q_cov = -pars_cov[1]/k_cov
                covfit = pars_cov[0]*qs**2 + pars_cov[1]*qs + pars_cov[2]
                #rmsd = np.sqrt( ((tot - ei - covfit)**2).sum()/len(totharm) )
                k_arr.append(k)
                k_cov_arr.append(k_cov)
                q_arr.append(q)
                q_cov_arr.append(q_cov)
            self.k[icname]     = Array(k_arr    , unit=system.kunit[icname])
            self.k_cov[icname] = Array(k_cov_arr, unit=system.kunit[icname])
            self.q[icname]     = Array(q_arr    , unit=system.qunit[icname])
            self.q_cov[icname] = Array(q_cov_arr, unit=system.qunit[icname])

    def html(self, icname, attrs, suffix=None):
        """
            Returns strings containing the statistical information of all attributes 
            in <attrs> for the given icname. Allowed attributes: q, q_cov, k, k_cov
        """
        fmts  = {
            'q'     :   ('q<sub> </sub> = ' ,'% 8.3f &#177 %7.3f (%i)'),
            'q_cov' :   ('q<sub>c</sub> = ' ,'% 8.3f &#177 %7.3f (%i)'),
            'k'     :   ('k<sub> </sub> = ' ,'% 8.1f &#177 %5.1f (%i)'),
            'k_cov' :   ('k<sub>c</sub> = ' ,'% 8.1f &#177 %5.1f (%i)'),
        }
        if '%s' in icname: icname = icname %( (self.name.split('/')[0].upper(),)*icname.count('%') )
        result = ''
        for attr in attrs:
            if len(attrs)>1: result += fmts[attr][0]
            result += getattr(self, attr)[icname].html(fmts[attr][1])
            if len(result)>0:
                if suffix is not None: result += suffix
                result += '<br>'
        return result
