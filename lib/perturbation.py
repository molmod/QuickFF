#! /usr/bin/env python

from molmod.units import *
from molmod.periodic import periodic as pt
from molmod.io.xyz import XYZWriter
import numpy as np, matplotlib.pyplot as pp
from scipy.optimize import minimize

from fftable import FFTable
from tools import *
from evaluators import *

__all__ = ['PerturbationTheory', 'RelaxedGeometryPT']


class PerturbationTheory(object):
    def __init__(self, skip_diheds):
        self.skip_diheds = skip_diheds
        self.description = 'BasePerturbationTheory (NOT FOR DIRECT USE)'
    
    def estimate(self, system):
        """
            Estimate force constant and rest value of all ics in the system
            directly from the hessian using perturbation theory.
        """
        print 'PERTUR ESTIM: estimate all pars using %s' %self.description
        fctab = FFTable(system.icnames, system.units)
        for icname, ics in system.ics.iteritems():
            if self.skip_diheds and icname.split('/')[0] in ['dihedral', 'dihed', 'torsion']:
                fctab.add(icname, [0.0]*len(ics), [0.0]*len(ics))
                continue
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


    def plot_single(self, system, icname, qunit='au', eunit='kjmol', kunit=None):
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
    def __init__(self, skip_diheds=True, energy_penalty=1.0, ic_penalty=1.0, dq_rel=0.1, qsteps = 101,
                 bond_thresshold=1.0*angstrom, bend_thresshold=1.0*deg, dihed_thresshold=1.0*deg):
        PerturbationTheory.__init__(self, skip_diheds)
        self.energy_penalty = energy_penalty
        self.ic_penalty = ic_penalty
        self.dq_rel = dq_rel
        self.qsteps = qsteps
        self.Dr = bond_thresshold
        self.Dt = bend_thresshold
        self.Dp = dihed_thresshold
        self.description = """Relaxed Geometry Perturbation Theory with:
               - an energy cost with weight of %.3e
               - an ic deviation cost with weight of %.3e
                    Err(r)     = %.3e A
                    Err(theta) = %.3e deg
                    Err(psi)   = %.3e deg   
        """ %(self.energy_penalty, self.ic_penalty, self.Dr/angstrom, self.Dt/deg, self.Dp/deg)
    
    def get_qdm(self, system, ic):
        """
            Calculate the Internal Coordinate Deveation Matrix (QDM).
            If sandwiched between a geometry perturbation vector, this
            represents the weighted sum of the deviations of the
            internal coordinates, except for the ic given in args, 
            from their non-perturbed values.
            
        """
        QDM = np.zeros([3*system.Natoms, 3*system.Natoms], float)
        indexes = system.get_neighbors(list(ic.indexes), depth=-1)
        for icname2 in system.icnames:
            for ic2 in system.ics[icname2]:
                if ic2.name==ic.name: continue
                include=True
                for i in ic2.indexes:
                    if not i in indexes:
                        include=False
                        break
                if not include: continue
                if icname2.split('/')[0] in ['bond', 'dist']: sigma = self.Dr
                elif icname2.split('/')[0] in ['bend', 'angle']: sigma = self.Dt
                elif icname2.split('/')[0] in ['dihedral', 'dihed', 'torsion']: sigma = self.Dp
                else: raise ValueError('Invalid icname, recieved %s' %icname2)
                print 'DEBUG: included %s in QDM with sigma = %.3e au' %(ic2.name, sigma)
                QDM += np.outer(ic2.grad(system.sample['coordinates']), ic2.grad(system.sample['coordinates']))/(sigma**2)
        VTx, VTy, VTz = global_translation(system.sample['coordinates'])
        VRx, VRy, VRz = global_rotation(system.sample['coordinates'])
        U, S, Vt = np.linalg.svd( np.array([VTx, VTy, VTz, VRx, VRy, VRz]).transpose() )
        P = np.dot(U[:,6:], U[:,6:].T)
        QDM = np.dot(P, np.dot(QDM, P))
        evals, evecs = np.linalg.eigh(QDM)
        print 'DEBUG: eigenvalues of QDM:', evals
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
        coords0 = system.sample['coordinates']
        for iic, ic in enumerate(system.ics[icname]):
            print 'PERTUR TRAJC: %s' %ic.name
            qarray = ( 1.0-self.dq_rel + 2.0*self.dq_rel/(self.qsteps-1)*np.array(range(self.qsteps),float) )*ic.value(system.sample['coordinates'])
            trajectory = []
            #Define cost function (and its derivative) that needs to be minimized
            QDM = self.get_qdm(system, ic)
            hess = system.totmodel.hess.reshape([3*system.Natoms, 3*system.Natoms])
            def chi(dx):
                ic = 0.5*np.dot(dx.T, np.dot(QDM, dx))
                ener = 0.5*np.dot(dx.T, np.dot(hess, dx))
                return self.ic_penalty*ic + self.energy_penalty*ener
            def dchi_dx(dx):
                ic = np.dot(QDM, dx)
                ener = np.dot(hess, dx)
                return self.ic_penalty*ic + self.energy_penalty*ener
            for iq, q in enumerate(qarray):
                #Define the constraint under which the cost function needs to be minimized
                constraints = (
                    {'type': 'eq',
                     'fun' : lambda dx: ic.value(coords0 + dx.reshape((-1, 3))) - q,
                     'jac' : lambda dx: ic.grad(coords0 + dx.reshape((-1, 3)))},
                )
                dx0 = np.zeros(3*system.Natoms, float)
                result = minimize(chi, dx0, method='SLSQP', jac=dchi_dx, constraints=constraints, tol=1e-9)
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
