#! /usr/bin/env python

from molmod.units import *
from molmod.periodic import periodic as pt
from molmod.io.xyz import XYZWriter
import numpy as np, matplotlib.pyplot as pp

from fftable import FFTable
from tools import fitpar, add_plot
from evaluators import *
from ffit import FFitProgram

__all__ = ['estimate', 'calculate_perturbation', 'analyze_perturbation', 'plot_perturbation']


def estimate(system, coupling=None, free_depth=0, spring=10.0*kjmol/angstrom**2):
    print 'PERTUR ESTIM: calculating harmonic ff pars directly from hessian'
    fctab = FFTable(system.icnames, system.units)
    for icname, ics in system.ics.iteritems():
        kdata = []
        qdata = []
        vs = calculate_perturbation(system, icname, coupling=coupling, free_depth=free_depth, spring=spring)
        if   icname.split('/')[0] in ['dist', 'bond'] : amplitude_factor = 0.1
        elif icname.split('/')[0] in ['angle', 'bend']: amplitude_factor = 0.01
        elif icname.split('/')[0] in ['dihedral', 'dihed', 'torsion']: amplitude_factor = 0.01
        else: raise ValueError('Invalid ictype: %s' %icname.split('/')[0])
        for iv, v in enumerate(vs):
            amplitude = amplitude_factor*ics[iv].value(system.sample['coordinates'])
            values = analyze_perturbation(
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



def calculate_perturbation(system, icname, coupling=None, free_depth=0, spring=10.0*kjmol/angstrom**2):
    """
        Calculate the perturbation on the geometry resulting from perturbing
        along the ics of <icname>. Coupling is a numpy array with the 
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
    print 'PERTUR CALC : calculation perturbation vector for %s with coupling=%s, free_depth=%i and spring=%.3f kjmol/A^2' %(
        icname, coupling, free_depth, spring/(kjmol/angstrom**2)
    )
    ics = system.ics[icname]
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
        qgrad_coupled = np.zeros(3*system.Natoms, float)
        for i, qgrad in enumerate(qgrads):
            qgrad_coupled += qgrad*coupling[i]
        for i in xrange(len(ics)):
            qgrads[i] = qgrad_coupled
    vs = []
    for qgrad in qgrads:
        if free_depth==0:
            ihess = system.totmodel.ihess.reshape([3*system.Natoms, 3*system.Natoms])
        else:
            indices = [i for i in xrange(system.Natoms) if np.linalg.norm(qgrad.reshape([system.Natoms, 3])[i])>1e-3]
            free_indices = system.get_neighbors(indices, depth=free_depth)
            ihess = system.totmodel.get_constrained_ihess(free_indices, spring=spring).reshape([3*system.Natoms, 3*system.Natoms])
        v = np.dot(ihess, qgrad)
        kl = np.dot(qgrad.T, np.dot(ihess, qgrad))
        vs.append(v.reshape((-1,3))/kl)
    return vs


def analyze_perturbation(system, v, evaluators=[], fn_xyz=None, amplitude=0.5, steps=101):
    if fn_xyz is not None:
        symbols = np.array([pt[n].symbol for n in system.sample['numbers']])
        traj = XYZWriter(file(fn_xyz, 'w'), symbols)
    values = [[] for i in xrange(len(evaluators))]
    for n in xrange(steps):
        pre = amplitude*np.sin(-np.pi/2+np.pi*n/(steps-1))
        coords = system.sample['coordinates'] + pre*v
        for i, evaluator in enumerate(evaluators):
            values[i].append(evaluator(system, coords))
        if fn_xyz is not None: traj.dump('frame %i' %n, coords)
    if fn_xyz is not None: del(traj)
    return values


def plot_perturbation(system, icname, relative_amplitude=0.1, coupling=None, free_depth=0, spring=10*kjmol/angstrom**2, qunit='au', eunit='kjmol', kunit=None):
    if kunit is None: kunit = '%s/%s**2' %(eunit, qunit)
    vs = calculate_perturbation(system, icname, coupling=coupling, free_depth=free_depth, spring=spring)
    pp.clf()
    fig, axs = pp.subplots(len(vs), 1)
    if len(vs)==1: axs = np.array([axs])
    
    for i, v in enumerate(vs):
        fn_xyz = 'traj-%s-%i.xyz' %(icname.split('/')[1], i)
        amplitude = system.ics[icname][i].value(system.sample['coordinates'])*relative_amplitude
        values = analyze_perturbation(
            system, v, 
            evaluators=[ic_evaluator(icname, i), energy_evaluator('totmodel'), energy_evaluator('eimodel')],
            fn_xyz=fn_xyz, amplitude=amplitude, steps=101
        )
        qs = np.array(values[0])
        total = np.array(values[1])
        ei = np.array(values[2])
        if coupling is not None: total /= len(vs)
        pars = fitpar(qs, total-ei, rcond=1e-6)
        k = 2*pars[0]
        q0 = -pars[1]/k
        fit = pars[0]*qs**2 + pars[1]*qs + pars[2]
        rmsd = np.sqrt( ((total - ei - fit)**2).sum()/len(total) )
        title = "perturbation %s:\n" %('-'.join(['%s[%i]' %(system.sample['ffatypes'][atindex],atindex) for atindex in system.ics[icname][i].indexes])) \
              + "---------------------------------------------\n" \
              + "Fit:   k  = % 4.0f %s   q0 = % 6.2f %s\n" %(k/parse_unit(kunit), kunit, q0/parse_unit(qunit), qunit) \
              + "---------------------------------------------\n" \
              + "RMS(Total - EI - Fit) = %.3e" %(rmsd/parse_unit(eunit))
        curves = [
            [qs, total, 'b-', 'AI Total Energy'       ],
            [qs, ei   , 'r-', 'Electrostatic Energy'  ],
            [qs, fit  , 'g-', 'Fitted Covalent Energy'],
        ]
        add_plot(axs[i], curves, title=title, xunit=qunit, yunit=eunit)

    fig.set_size_inches([8, 8*len(vs)])
    fig.tight_layout()
    pp.savefig('energy-%s.pdf' %(icname.split('/')[1]))
