#! /usr/bin/env python

import numpy as np
from fctable import *
from tools import fitpar
from evaluators import *

__all__ = ['estimate', 'geometry_perturbation', 'analyze_trajectory']

def estimate(system, coupling=None, free_depth=0, spring=10.0*kjmol/angstrom**2):
    print '   ESTIMATOR: estimating ff pars'
    fctab = FCTable(system.icnames, system.units)
    for icname, ics in system.ics.iteritems():
        kdata = []
        qdata = []
        vs = geometry_perturbation(system, ics, coupling=coupling, free_depth=free_depth, spring=spring)
        if   icname.split('/')[0]=='bond' : amplitude_factor = 0.1
        elif icname.split('/')[0]=='angle': amplitude_factor = 0.01
        elif icname.split('/')[0]=='dihed': amplitude_factor = 0.01
        else: raise ValueError('Invalid ictype: %s' %icname.split('/')[0])
        for iv, v in enumerate(vs):
            amplitude = amplitude_factor*ics[iv].value(system.sample['coordinates'])
            values = analyze_trajectory(
                system, v, evaluators=[ic_evaluator(icname, iv), energy_evaluator('totmodel'), energy_evaluator('eimodel')], amplitude=amplitude
            )
            qs = np.array(values[0])
            tot = np.array(values[1])
            ei = np.array(values[2])
            if coupling is not None:
                tot /= len(ics)
                ei /= len(ics)
            pars = fitpar(qs, tot-ei, rcond=1e-6)
            k = 2*pars[0]
            q = -pars[1]/k
            fit = pars[0]*qs**2 + pars[1]*qs + pars[2]
            kdata.append(k)
            qdata.append(q)
        fctab.add(icname, kdata, qdata)
    return fctab


def geometry_perturbation(system, ics, coupling=None, free_depth=0, spring=10.0*kjmol/angstrom**2):
    """
        Get the perturbation on the geometry resulting from perturbing along
        the ics. Coupling is a numpy array with the 
        coefficients of the linear combination describing the coupling.
        If free_depth is None, the hessian remains unbiased for calculating 
        the perturbed geometry. If free_depth is e.g. equal to 3, all atoms 
        that are separated by more then 3 atoms from the atoms in the ic 
        under consideration are fixed with a spring to their original position
        (with strength defined in spring).
        
        The following strings are also supported for coupling:
            ``symmetric``   A symmetric coupling 
                            i.e. coupling = np.ones(N)/np.sqrt(N)
                            with N the number of ics for the given icname
    """
    if coupling=='symmetric':
        coupling = np.ones(len(ics), float)/len(ics)
    elif isinstance(coupling, np.ndarray):
        assert len(ics)==len(coupling)
    elif coupling is not None:
        raise NotImplementedError('Unsupported coupling: %s' %(str(coupling)))
    qgrads = []
    for i, ic in enumerate(ics):
        qgrads.append(ic.grad(system.sample['coordinates']))
    if coupling is not None:
        qgrad_coupled = np.zeros(3*system.Nat, float)
        for i, qgrad in enumerate(qgrads):
            qgrad_coupled += qgrad*coupling[i]
        for i in xrange(len(ics)):
            qgrads[i] = qgrad_coupled
    vs = []
    for qgrad in qgrads:
        if free_depth==0:
            ihess = system.totmodel.ihess
        else:
            indices = [i for i in xrange(system.Nat) if np.linalg.norm(qgrad.reshape([system.Nat, 3])[i])>1e-3]
            free_indices = system._get_neighbors(indices, depth=free_depth)
            ihess = system.totmodel.get_constrained_ihess(free_indices, spring=spring)
        v = np.dot(ihess, qgrad)
        kl = np.dot(qgrad.T, np.dot(ihess, qgrad))
        vs.append(v.reshape((-1,3))/kl)
    return vs

def analyze_trajectory(system, v, evaluators=[], fn_xyz=None, amplitude=0.25*angstrom, steps=101):
    symbols = np.array([pt[n].symbol for n in system.sample['numbers']])
    if fn_xyz is not None: traj = XYZWriter(file(fn_xyz, 'w'), symbols)
    values = [[] for i in xrange(len(evaluators))]
    for n in xrange(steps):
        pre = amplitude*np.sin(-np.pi/2+np.pi*n/(steps-1))
        coords = system.sample['coordinates'] + pre*v
        for i, evaluator in enumerate(evaluators):
            values[i].append(evaluator(self, coords))
        if fn_xyz is not None: traj.dump('frame %i' %n, coords)
    if fn_xyz is not None: del(traj)
    return values
