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
    def __init__(self, system, skip_diheds):
        self.system = system
        self.skip_diheds = skip_diheds
        self.description = 'BasePerturbationTheory (NOT FOR DIRECT USE)'
    
    def estimate(self):
        """
            Estimate force constant and rest value of all ics in the system
            directly from the hessian using perturbation theory.
        """
        print 'PERTUR ESTIM: estimate all pars using %s' %self.description
        fctab = FFTable(self.system.icnames)
        for icname, ics in self.system.ics.iteritems():
            unit = {'q': ics[0].qunit, 'k': ics[0].kunit}
            if self.skip_diheds and icname.split('/')[0] in ['dihedral', 'dihed', 'torsion']:
                fctab.add(icname, [0.0]*len(ics), [0.0]*len(ics))
                continue
            kdata = []
            qdata = []
            trajectories = self.perturbation_trajectory(icname)
            for itraj, trajectory in enumerate(trajectories):
                evaluators = [ic_evaluator(name=icname, i=itraj), energy_evaluator('totmodel'), energy_evaluator('eimodel')]
                qs, tot, ei = self.analyze(trajectory, evaluators=evaluators)
                pars = fitpar(qs, tot-ei, rcond=1e-6)
                k = 2*pars[0]
                q = -pars[1]/k
                kdata.append(k)
                qdata.append(q)
            fctab.add(icname, kdata, qdata, unit=unit)
        print
        return fctab

    def perturbation_trajectory(self, icname):
        """
            Calculate the perturbation on the geometry when perturbing
            along <icname>. This is calculated in a derived class.
        """
        raise NotImplementedError


    def analyze(self, trajectory, evaluators=[], fn_xyz=None):
        if fn_xyz is not None:
            symbols = np.array([pt[n].symbol for n in self.system.sample['numbers']])
            xyzwriter = XYZWriter(file(fn_xyz, 'w'), symbols)
        values = [[] for i in xrange(len(evaluators))]
        for dx in trajectory:
            coords = self.system.sample['coordinates'] + dx
            for i, evaluator in enumerate(evaluators):
                values[i].append(evaluator(self.system, coords))
            if fn_xyz is not None: xyzwriter.dump('frame %i' %n, coords)
        if fn_xyz is not None: del(xyzwriter)
        return np.array(values)

    
    def plot_icname(self, icname, eunit='kjmol'):
        print 'PERTUR SINGL: Calculating perturbation trajectory for %s' %self.description
        trajectories = self.perturbation_trajectory(icname)
        N = len(trajectories)
        
        for iic, ic in enumerate(self.system.ics[icname]):
            self.plot_ic_constrained([N,5,5*iic+1], ic, trajectories[iic])
            QDM, other_ics = self.get_qdm(ic, return_other_ics=True)
            bonds = [other_ic for other_ic in other_ics if other_ic.name.split('/')[0] in ['bond', 'dist']]
            if len(bonds)>0: self.plot_ic_other_ics([N,5,5*iic+2], ic, bonds, trajectories[iic])
            bends = [other_ic for other_ic in other_ics if other_ic.name.split('/')[0] in ['bend', 'angle']]
            if len(bends)>0: self.plot_ic_other_ics([N,5,5*iic+3], ic, bends, trajectories[iic])
            diheds = [other_ic for other_ic in other_ics if other_ic.name.split('/')[0] in ['dihed', 'dihedral', 'torsion']]
            if len(diheds)>0: self.plot_ic_other_ics([N,5,5*iic+4], ic, diheds, trajectories[iic])
            self.plot_ic_energies([N,5,5*iic+5], ic, trajectories[iic], eunit=eunit)
        
        fig.set_size_inches([8*5, 8*N])
        fig.tight_layout()
        pp.savefig( '%s.pdf' %(icname.replace('/', '-')) )


    def plot_ic_energies(self, plc, ic, trajectory, eunit='kjmol'):
        fn_xyz = '%s.xyz' %(ic.name.replace('/','-'))
        evaluators = [ic_evaluator(ic=ic), energy_evaluator('totmodel'), energy_evaluator('eimodel')]
        qs, tot, ei = self.analyze(trajectory, evaluators=evaluators, fn_xyz=fn_xyz)
        pars_tot = fitpar(qs, tot, rcond=1e-6)
        k_tot    = 2*pars_tot[0]
        q0_tot   = -pars_tot[1]/k_tot
        pars_ei = fitpar(qs, ei, rcond=1e-6)
        k_ei    = 2*pars_ei[0]
        q0_ei   = -pars_ei[1]/k_ei
        pars = fitpar(qs, tot-ei, rcond=1e-6)
        k    = 2*pars[0]
        q0   = -pars[1]/k
        fit  = pars[0]*qs**2 + pars[1]*qs + pars[2]
        rmsd = np.sqrt( ((tot - ei - fit)**2).sum()/len(tot) )
        title = "Fit Tot: k  = % 4.0f %s   q0 = % 6.2f %s\n" %(k_tot/parse_unit(ic.kunit), ic.kunit, q0_tot/parse_unit(ic.qunit), ic.qunit) \
              + "Fit EI : k  = % 4.0f %s   q0 = % 6.2f %s\n" %(k_ei/parse_unit(ic.kunit), ic.kunit, q0_ei/parse_unit(ic.qunit), ic.qunit) \
              + "Fit Cov: k  = % 4.0f %s   q0 = % 6.2f %s\n" %(k/parse_unit(ic.kunit), ic.kunit, q0/parse_unit(ic.qunit), ic.qunit) \
              + "------------------------------------------------------------\n" \
              + "RMS(Tot - EI - Fit) = %.3e %s" %(rmsd/parse_unit(eunit), eunit)
        curves = [
            [qs, tot, 'b-', 'AI Total Energy'       ],
            [qs, ei , 'r-', 'Electrostatic Energy'  ],
            [qs, fit, 'g-', 'Fitted Covalent Energy'],
        ]
        ax = pp.subplot(plc[0], plc[1], plc[2])
        add_plot(ax, curves, title=title, xunit=ic.qunit, yunit=eunit)




class RelaxedGeometryPT(PerturbationTheory):
    def __init__(self, system, skip_diheds=True, energy_penalty=1.0, ic_penalty=1.0, dq_rel=0.1, qsteps = 101,
                 bond_thresshold=0.01*angstrom, bend_thresshold=0.1*deg, dihed_thresshold=1.0*deg):
        PerturbationTheory.__init__(self, system, skip_diheds)
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
    
    def get_qdm(self, ic, return_other_ics=False):
        """
            Calculate the Internal Coordinate Deveation Matrix (QDM).
            If sandwiched between a geometry perturbation vector, this
            represents the weighted sum of the deviations of the
            internal coordinates, except for the ic given in args, 
            from their non-perturbed values.
            
        """
        QDM = np.zeros([3*self.system.Natoms, 3*self.system.Natoms], float)
        indexes = self.system.get_neighbors(list(ic.indexes), depth=-1)
        other_ics = []
        for icname2 in self.system.icnames:
            for ic2 in self.system.ics[icname2]:
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
                other_ics.append(ic2)
                QDM += np.outer(ic2.grad(self.system.sample['coordinates']), ic2.grad(self.system.sample['coordinates']))/(sigma**2)
        VTx, VTy, VTz = global_translation(self.system.sample['coordinates'])
        VRx, VRy, VRz = global_rotation(self.system.sample['coordinates'])
        U, S, Vt = np.linalg.svd( np.array([VTx, VTy, VTz, VRx, VRy, VRz]).transpose() )
        P = np.dot(U[:,6:], U[:,6:].T)
        QDM = np.dot(P.T, np.dot(QDM, P))
        evals, evecs = np.linalg.eigh(QDM)
        if return_other_ics:
            return QDM, other_ics
        else:
            return QDM
                    
    
    def perturbation_trajectory(self, icname):
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
        trajectories = []
        coords0 = self.system.sample['coordinates']
        for iic, ic in enumerate(self.system.ics[icname]):
            print 'PERTUR TRAJC: %s' %ic.name
            qarray = ( 1.0-self.dq_rel + 2.0*self.dq_rel/(self.qsteps-1)*np.array(range(self.qsteps),float) )*ic.value(coords0)
            trajectory = []
            #Define cost function (and its derivative) that needs to be minimized
            QDM = self.get_qdm(ic)
            hess = 1e-3*self.system.totmodel.hess.reshape([3*self.system.Natoms, 3*self.system.Natoms])
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
                dx0 = np.zeros(3*self.system.Natoms, float)
                result = minimize(chi, dx0, method='SLSQP', jac=dchi_dx, constraints=constraints, tol=1e-9)
                trajectory.append(result.x.reshape((-1, 3)))
            assert len(trajectory)==self.qsteps
            trajectories.append(trajectory)
        return trajectories
    
   
    def plot_icname(self, icname, eunit='kjmol'):
        trajectories = self.perturbation_trajectory(icname)
        print 'PERTUR SINGL: Plot analysis of perturbation trajectory for %s' %(icname)
        print
        print '              %s' %self.description
        N = len(trajectories)
        
        for iic, ic in enumerate(self.system.ics[icname]):
            self.plot_ic_costs([N,6,iic], ic, trajectories[iic])
            self.plot_ic_constrained([N,6,6*iic+2], ic, trajectories[iic])
            QDM, other_ics = self.get_qdm(ic, return_other_ics=True)
            bonds = [other_ic for other_ic in other_ics if other_ic.name.split('/')[0] in ['bond', 'dist']]
            self.plot_ic_other_ics([N,6,6*iic+3], ic, bonds, trajectories[iic], title='Free Bond Lengths', ylabel='Bond length [%s]', yunit='A')
            bends = [other_ic for other_ic in other_ics if other_ic.name.split('/')[0] in ['bend', 'angle']]
            self.plot_ic_other_ics([N,6,6*iic+4], ic, bends, trajectories[iic], title='Free Bending Angles', ylabel='Bending Angle [%s]', yunit='deg')
            diheds = [other_ic for other_ic in other_ics if other_ic.name.split('/')[0] in ['dihed', 'dihedral', 'torsion']]
            self.plot_ic_other_ics([N,6,6*iic+5], ic, diheds, trajectories[iic], title='Free Dihedral Angles', ylabel='Dihedral Angle [%s]', yunit='deg')
            self.plot_ic_energies([N,6,6*iic+6], ic, trajectories[iic], eunit=eunit)
        
        fig = pp.gcf()
        fig.set_size_inches([8*7, 8*N])
        fig.tight_layout()
        pp.savefig( '%s.pdf' %(icname.replace('/', '-')) )
    
    
    def plot_ic_costs(self, plc, ic, trajectory, eunit='kjmol'):
        'plc = plot configuration (array of form [Nrows, Ncols, current fig index])'
        QDM = self.get_qdm(ic)
        hess = self.system.totmodel.hess.reshape([3*self.system.Natoms, 3*self.system.Natoms])
        q = np.zeros(len(trajectory), float)
        icdev  = np.zeros(len(trajectory), float)
        energy = np.zeros(len(trajectory), float)
        total  = np.zeros(len(trajectory), float)
        for i, dx in enumerate(trajectory):
            q[i]      = ic.value(self.system.sample['coordinates'] + dx)
            icdev[i]  = 0.5*np.dot(dx.reshape([3*self.system.Natoms]).T, np.dot(QDM, dx.reshape([3*self.system.Natoms])))
            energy[i] = 0.5*np.dot(dx.reshape([3*self.system.Natoms]).T, np.dot(hess, dx.reshape([3*self.system.Natoms])))
        curves = [
            [q, icdev , 'b-', '_nolegend_'],
            [q, energy, 'g-', '_nolegend_'],
        ]
        ax = pp.subplot(2*plc[0], plc[1], 2*plc[2]*plc[1]+1)
        title = "Relative Deviation Cost"
        add_plot(ax, [curves[0]], title=title, xunit=ic.qunit, yunit='1.0', xlabel='%s [%%s]' %(ic.name.split('/')[0]), ylabel='Deviation [-]')
        ax = pp.subplot(2*plc[0], plc[1], (2*plc[2]+1)*plc[1]+1)
        add_plot(ax, [curves[1]], title='Energy Cost', xunit=ic.qunit, yunit=eunit, xlabel='%s [%%s]' %(ic.name.split('/')[0]), ylabel='Energy [%s]')
    
    
    def plot_ic_constrained(self, plc, ic, trajectory):
        q = np.zeros(len(trajectory), float)
        for i, dx in enumerate(trajectory):
            q[i]      = ic.value(self.system.sample['coordinates'] + dx)
        curves = [[np.array(range(len(trajectory))), q , 'b-', '_nolegend_']]
        ax = pp.subplot(plc[0], plc[1], plc[2])
        add_plot(ax, curves, title='Constrained IC', xunit='1.0', yunit=ic.qunit, xlabel='Step [-]', ylabel='%s [%%s]' %(ic.name.split('/')[0]))
    
    
    def plot_ic_other_ics(self, plc, ic, other_ics, trajectory, title='Other ICs', ylabel='IC [%%s]', yunit='au'):
        styles = ['b-', 'r-', 'g-', 'c-', 'm-', 'y-', 'b--', 'r--', 'g--', 'c--', 'm--', 'y--']
        curves = []
        for iother, other_ic in enumerate(other_ics):
            q  = np.zeros(len(trajectory), float)
            q2 = np.zeros(len(trajectory), float)
            for i, dx in enumerate(trajectory):
                q[i]  = ic.value(self.system.sample['coordinates'] + dx)
                q2[i] = other_ic.value(self.system.sample['coordinates'] + dx)
            curves.append([q, q2 , styles[iother], other_ic.name.split('/')[1]])
        ax = pp.subplot(plc[0], plc[1], plc[2])
        add_plot(ax, curves, title=title, xunit=ic.qunit, yunit=yunit, xlabel='%s [%%s]' %(ic.name.split('/')[0]), ylabel=ylabel)




def add_plot(ax, curves, title='', xunit='au', yunit='au', xlabel='IC [%s]', ylabel='Energy [%s]'):
    for x, y, style, label in curves:
        ax.plot(x/parse_unit(xunit), y/parse_unit(yunit), style, label=label)
    if '%s' in xlabel: xlabel = xlabel %xunit
    if '%s' in ylabel: ylabel = ylabel %yunit
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid()
    ax.set_title(title)
    ax.legend(loc='best')
