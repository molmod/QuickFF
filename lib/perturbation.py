#! /usr/bin/env python

from molmod.units import *
from molmod.periodic import periodic as pt
from molmod.io.xyz import XYZWriter
import numpy as np, matplotlib.pyplot as pp
from scipy.optimize import minimize

from fftable import FFTable
from tools import fitpar
from evaluators import *

__all__ = ['PerturbationTheory', 'RelaxedGeometryPT']


class PerturbationTheory(object):
    def __init__(self):
        self.description = 'BasePerturbationTheory (NOT FOR DIRECT USE)'
    
    def estimate(self, system):
        """
            Estimate force constant and rest value of all ics in the system
            directly from the hessian using perturbation theory.
        """
        print 'PERTUR ESTIM: estimate all pars using %s' %self.description
        fctab = FFTable(system.icnames, system.units)
        for icname, ics in system.ics.iteritems():
            kdata = []
            qdata = []
            trajectories = self.perturbation_trajectory(system, icname)
            for itraj, trajectory in enumerate(trajectories):
                evaluators = [ic_evaluator(icname, itraj), energy_evaluator('totmodel'), energy_evaluator('eimodel')]
                values = self.analyze(system, trajectory, evaluators=evaluators)
                qs = np.array(values[0])
                tot = np.array(values[1])
                ei = np.array(values[2])
                pars = fitpar(qs, tot-ei, rcond=1e-6)
                k = 2*pars[0]
                q = -pars[1]/k
                fit = pars[0]*qs**2 + pars[1]*qs + pars[2]
                kdata.append(k)
                qdata.append(q)
            fctab.add(icname, kdata, qdata)
        print
        return fctab

    def perturbation_trajectory(self, system, icname):
        """
            Calculate the perturbation on the geometry when perturbing
            along <icname>. This is calculated in a derived class.
        """
        raise NotImplementedError


    def analyze(self, system, trajectory, evaluators=[], fn_xyz=None):
        if fn_xyz is not None:
            symbols = np.array([pt[n].symbol for n in system.sample['numbers']])
            xyzwriter = XYZWriter(file(fn_xyz, 'w'), symbols)
        values = [[] for i in xrange(len(evaluators))]
        for dx in trajectory:
            coords = system.sample['coordinates'] + dx
            for i, evaluator in enumerate(evaluators):
                values[i].append(evaluator(system, coords))
            if fn_xyz is not None: xyzwriter.dump('frame %i' %n, coords)
        if fn_xyz is not None: del(xyzwriter)
        return values


    def plot_single(self, system, icname,start=0.9, end=1.1, steps=101, qunit='au', eunit='kjmol', kunit=None):
        print 'PERTUR SINGL: Estimate single par using %s' %self.description
        if kunit is None: kunit = '%s/%s**2' %(eunit, qunit)
        trajectories = self.perturbation_trajectory(system, icname)
        fig, axs = pp.subplots(len(trajectories), 1)
        if len(trajectories)==1: axs = np.array([axs])
        
        for itraj, trajectory in enumerate(trajectories):
            fn_xyz = 'traj-%s-%i.xyz' %(icname.split('/')[1], itraj)
            evaluators = [ic_evaluator(icname, itraj), energy_evaluator('totmodel'), energy_evaluator('eimodel')]
            values = self.analyze(system, trajectory, evaluators=evaluators, fn_xyz=fn_xyz)
            qs    = np.array(values[0])
            total = np.array(values[1])
            ei    = np.array(values[2])
            pars_ref = fitpar(qs, total, rcond=1e-6)
            k_ref    = 2*pars_ref[0]
            q0_ref   = -pars_ref[1]/k_ref
            pars_ei = fitpar(qs, ei, rcond=1e-6)
            k_ei    = 2*pars_ei[0]
            if pars_ei[1]==0:
                q0_ei = 0.0
            else:
                q0_ei = -pars_ei[1]/k_ei
            pars = fitpar(qs, total-ei, rcond=1e-6)
            k    = 2*pars[0]
            q0   = -pars[1]/k
            fit  = pars[0]*qs**2 + pars[1]*qs + pars[2]
            rmsd = np.sqrt( ((total - ei - fit)**2).sum()/len(total) )
            title = "perturbation %s:\n" %('-'.join(['%s[%i]' %(system.sample['ffatypes'][atindex],atindex) for atindex in system.ics[icname][itraj].indexes])) \
                  + "---------------------------------------------\n" \
                  + "Fit Tot: k  = % 4.0f %s   q0 = % 6.2f %s\n" %(k_ref/parse_unit(kunit), kunit, q0_ref/parse_unit(qunit), qunit) \
                  + "Fit EI : k  = % 4.0f %s   q0 = % 6.2f %s\n" %(k_ei/parse_unit(kunit), kunit, q0_ei/parse_unit(qunit), qunit) \
                  + "Fit Cov: k  = % 4.0f %s   q0 = % 6.2f %s\n" %(k/parse_unit(kunit), kunit, q0/parse_unit(qunit), qunit) \
                  + "---------------------------------------------\n" \
                  + "RMS(Total - EI - Fit) = %.3e %s" %(rmsd/parse_unit(eunit), eunit)
            curves = [
                [qs, total, 'b-', 'AI Total Energy'       ],
                [qs, ei   , 'r-', 'Electrostatic Energy'  ],
                [qs, fit  , 'g-', 'Fitted Covalent Energy'],
            ]
            add_plot(axs[itraj], curves, title=title, xunit=qunit, yunit=eunit)

        fig.set_size_inches([8, 8*len(trajectories)])
        fig.tight_layout()
        pp.savefig('energy-%s.pdf' %(icname.split('/')[1]))


class RelaxedGeometryPT(PerturbationTheory):
    def __init__(self, energy_penalty=1.0, ic_penalty=1.0, dq_rel=0.1, qsteps = 101,
                 bond_thresshold=1e-6*angstrom, bend_thresshold=1e-4*deg, dihed_thresshold=1e-2*deg):
        self.energy_penalty = energy_penalty
        self.ic_penalty = ic_penalty
        self.dq_rel = dq_rel
        self.qsteps = qsteps
        self.Dr = bond_thresshold
        self.Dt = bend_thresshold
        self.Dp = dihed_thresshold
        self.description = 'Relaxed Geometry Perturbation Theory with: \n'    + \
                           '               - an energy cost with weight of %.3e \n' %(self.energy_penalty)          + \
                           '               - an ic deviation cost with weight of %.3e \n' %(self.ic_penalty)        + \
                           '                    Err(r)     = %.3e A   \n' %(self.Dr/angstrom)                       + \
                           '                    Err(theta) = %.3e deg \n' %(self.Dt/deg)                            + \
                           '                    Err(psi)   = %.3e deg   ' %(self.Dp/deg)
    
    def calc_qdm(self, system, ic, coords=None):
        """
            Calculate the Internal Coordinate Deveation Matrix (QDM).
            If sandwiched between a geometry perturbation vector, this
            matrix represents the weighted sum of the deviations of the
            internal coordinates, except for the ic given in args, 
            from their non-perturbed values.
            
        """
        if coords is None: coords = system.sample['coordinates']
        QDM = np.zeros([3*system.Natoms, 3*system.Natoms], float)
        for icname2 in system.icnames:
            for ic2 in system.ics[icname2]:
                if ic2.name==ic.name:
                    continue
                if icname2.split('/')[0] in ['bond', 'dist']: sigma = self.Dr
                elif icname2.split('/')[0] in ['bend', 'angle']: sigma = self.Dt
                elif icname2.split('/')[0] in ['dihedral', 'dihed', 'torsion']: sigma = self.Dp
                else: raise ValueError('Invalid icname, recieved %s' %icname2)
                QDM += np.outer(ic2.grad(coords), ic2.grad(coords))/(sigma**2)
        return QDM
    
    
    def perturbation_trajectory(self, system, icname):
        """
            Calculate the perturbation on the geometry when perturbing
            along <icname> with magnitudes given in icrange. 
            
              icrange = [1.0 - self.dq_rel  ,  1.0 + self.dq_rel]*ic0 
              with a total of self.qsteps
            
            The perturbation on the geometry is calculated by minimizing
            a cost function containing the energy change and the deviation
            of all ics except the one under consideration. This minimization
            is performed under the constraint of a value q (in icrange)
            for the ic under consideration.
        """
        qunit = parse_unit(system.units[icname]['q'])
        trajectories = []
        for iic, ic in enumerate(system.ics[icname]):
            print 'PERTUR TRAJC: %s' %ic.name
            qarray = ( 1.0-self.dq_rel + 2.0*self.dq_rel/(self.qsteps-1)*np.array(range(self.qsteps),float) )*ic.value(system.sample['coordinates'])
            trajectory = []
            #Define cost function (and its derivative) that needs to be minimized
            QDM = self.calc_qdm(system, ic)
            def chi(dx):
                X = self.ic_penalty*0.5*np.dot(dx.T, np.dot(QDM, dx)) \
                  + self.energy_penalty*system.totmodel.get_energy(dx.reshape((-1, 3)))
                return X
            def dchi_dx(dx):
                X = self.ic_penalty*0.5*np.dot(QDM, dx) \
                  + self.energy_penalty*system.totmodel.get_gradient(dx.reshape((-1, 3))).reshape(3*system.Natoms)
                return X 
            for iq, q in enumerate(qarray):
                #Define the constraint under which the cost function needs to be minimized
                constraints = ({
                    'type': 'eq',
                    'fun' : lambda dx: ic.value(system.sample['coordinates'] + dx.reshape((-1, 3))) - q,
                    'jac' : lambda dx: ic.grad(system.sample['coordinates'] + dx.reshape((-1, 3))),
                })
                result = minimize(chi, np.zeros(3*system.Natoms, float), method='SLSQP', jac=dchi_dx, constraints=constraints, tol=1e-9)
                trajectory.append(result.x.reshape((-1, 3)))
            assert len(trajectory)==self.qsteps
            trajectories.append(trajectory)
        return trajectories
    
   


def add_plot(ax, curves, title='', xunit='au', yunit='au', xlabel='IC [%s]', ylabel='Energy [%s]'):
    for x, y, style, label in curves:
        ax.plot(x/parse_unit(xunit), y/parse_unit(yunit), style, label=label)
    if '%s' in xlabel: xlabel = xlabel %xunit
    if '%s' in ylabel: ylabel = ylabel %yunit
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid()
    ax.set_title(title)
    ax.legend(loc='lower left')
