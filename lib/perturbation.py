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
    def __init__(self, system, skip_diheds, print_descr=True):
        self.system = system
        self.skip_diheds = skip_diheds
        self.description = "PERTUR THEOR: BasePerturbationTheory (NOT FOR DIRECT USE)"
        if print_descr: print self.description
    
    def estimate(self):
        """
            Estimate force constant and rest value of all ics in the system
            directly from the hessian using perturbation theory.
        """
        print 'PERTUR ESTIM: estimate all pars'
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
        for idx, dx in enumerate(trajectory):
            coords = self.system.sample['coordinates'] + dx
            for i, evaluator in enumerate(evaluators):
                values[i].append(evaluator(self.system, coords))
            if fn_xyz is not None: xyzwriter.dump('frame %i' %idx, coords)
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
        title = "Energy Contributions"
        curves = [
            [qs, tot, 'b-', 'AI Total Energy'       ],
            [qs, ei , 'r-', 'Electrostatic Energy'  ],
            [qs, fit, 'g-', 'Fitted Covalent Energy'],
        ]
        ax = pp.subplot(plc[0], plc[1], plc[2])
        add_plot(ax, curves, title=title, xunit=ic.qunit, yunit=eunit)




class RelaxedGeometryPT(PerturbationTheory):
    def __init__(self, system, skip_diheds=True, energy_penalty=1.0, strain_penalty=1.0, dq_rel=0.05, qsteps = 51,
                 bond_thresshold=1.0, bend_thresshold=1.0, dihed_thresshold=1.0):
        PerturbationTheory.__init__(self, system, skip_diheds, print_descr=False)
        self.energy_penalty = energy_penalty
        self.strain_penalty = strain_penalty
        self.dq_rel = dq_rel
        self.qsteps = qsteps
        self.Dr = bond_thresshold
        self.Dt = bend_thresshold
        self.Dp = dihed_thresshold
        self.description = """PERTUR THEOR: Relaxed Geometry Perturbation Theory with:
        
               - an energy cost with weight of %.3e
               - an strain cost with weight of %.3e
                    Err(r)     = %.3e
                    Err(theta) = %.3e
                    Err(psi)   = %.3e  
        """ %(self.energy_penalty, self.strain_penalty, self.Dr, self.Dt, self.Dp)
        print self.description
    
    def strain_matrix(self, ic=None):
        """
            Calculate the Strain Matrix.
            If sandwiched between a geometry perturbation vector, this
            represents the weighted sum of the deviations of the
            internal coordinates, except for the ic given in args, 
            from their non-perturbed values.
            
        """
        strain = np.zeros([3*self.system.Natoms, 3*self.system.Natoms], float)
        coords0 = self.system.sample['coordinates']
        for icname, ics in self.system.ics.iteritems():
            if   icname.split('/')[0] in ['dist' , 'bond']:                 sigma = self.Dr
            elif icname.split('/')[0] in ['angle', 'bend']:                 sigma = self.Dt
            elif icname.split('/')[0] in ['dihed', 'dihedral', 'torsion']:  sigma = self.Dp
            else: raise ValueError('Illegal ic kind, recieved %s' %other.name)
            Gq = np.array([ic0.grad(coords0) for ic0 in ics if ic is None or ic0.name!=ic.name])
            if len(Gq)>0:
                U, S, Vt = np.linalg.svd(Gq)
                svals = np.array([s for s in S if s > 1e-6])
                rank = len(svals)
                V  = Vt.T[:,:rank]
                Vo = Vt.T[:,rank:]
                strain += np.dot(V, V.T) + 0.01*np.dot(Vo, Vo.T)/(3*self.system.Natoms)
        return strain
    
    
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
        for ic in self.system.ics[icname]:
            print 'PERTUR TRAJC: %s' %ic.name
            qarray = ( 1.0-self.dq_rel + 2.0*self.dq_rel/(self.qsteps-1)*np.array(range(self.qsteps),float) )*ic.value(coords0)
            trajectory = []
            #Define cost function (and its derivative) that needs to be minimized
            H = self.system.totmodel.hess.reshape([3*self.system.Natoms, 3*self.system.Natoms])
            S = self.strain_matrix(ic)
            def chi(dx):
                strain = 0.5*np.dot(dx.T, np.dot(S, dx))
                energy = 0.5*np.dot(dx.T, np.dot(H, dx))
                return self.strain_penalty*strain + self.energy_penalty*energy
            for iq, q in enumerate(qarray):
                #Define the constraint under which the cost function needs to be minimized
                constraints = (
                    {'type': 'eq',
                     'fun' : lambda dx: ic.value(coords0 + dx.reshape((-1, 3))) - q,
                     'jac' : lambda dx: ic.grad(coords0 + dx.reshape((-1, 3)))},
                )
                result = minimize(chi, np.zeros([3*self.system.Natoms], float), method='SLSQP', constraints=constraints, tol=1e-15)
                trajectory.append(result.x.reshape([-1, 3]))
            assert len(trajectory)==self.qsteps
            trajectories.append(trajectory)
        return trajectories
    
   
    def plot_icname(self, icname, eunit='kjmol'):
        trajectories = self.perturbation_trajectory(icname)
        print 'PERTUR PLOT : Plot analysis of perturbation trajectory for %s' %(icname)
        N = len(trajectories)
        
        for iic, (ic, trajectory) in enumerate(zip(self.system.ics[icname], trajectories)):
            self.plot_ic_costs([N,6,iic], ic, trajectory)
            self.plot_ic_constrained([N,6,6*iic+2], ic, trajectory)
            other_ics = self.system.get_ics(exclude_ic=[ic.name])
            bonds = [other_ic for other_ic in other_ics if other_ic.name.split('/')[0] in ['bond', 'dist']]
            self.plot_ic_other_ics([N,6,6*iic+3], ic, bonds, trajectory, title='Free Bond Lengths', ylabel='Bond length [%s]', yunit='A')
            bends = [other_ic for other_ic in other_ics if other_ic.name.split('/')[0] in ['bend', 'angle']]
            self.plot_ic_other_ics([N,6,6*iic+4], ic, bends, trajectory, title='Free Bending Angles', ylabel='Bending Angle [%s]', yunit='deg')
            diheds = [other_ic for other_ic in other_ics if other_ic.name.split('/')[0] in ['dihed', 'dihedral', 'torsion']]
            self.plot_ic_other_ics([N,6,6*iic+5], ic, diheds, trajectory, title='Free Dihedral Angles', ylabel='Dihedral Angle [%s]', yunit='deg')
            self.plot_ic_energies([N,6,6*iic+6], ic, trajectory, eunit=eunit)
        
        fig = pp.gcf()
        fig.set_size_inches([10*6, 10*N])
        pp.tight_layout()
        
        pp.savefig( '%s.pdf' %(icname.replace('/', '-')))
    
    
    def plot_ic_costs(self, plc, ic, trajectory, eunit='kjmol'):
        'plc = plot configuration (array of form [Nrows, Ncols, current fig index])'
        H = self.system.totmodel.hess.reshape([3*self.system.Natoms, 3*self.system.Natoms])
        S = self.strain_matrix(ic)
        q = np.zeros(len(trajectory), float)
        strain = np.zeros(len(trajectory), float)
        energy = np.zeros(len(trajectory), float)
        total  = np.zeros(len(trajectory), float)
        for i, dx in enumerate(trajectory):
            q[i]      = ic.value(self.system.sample['coordinates'] + dx)
            strain[i] = 0.5*np.dot(dx.reshape([3*self.system.Natoms]).T, np.dot(S, dx.reshape([3*self.system.Natoms])))
            energy[i] = 0.5*np.dot(dx.reshape([3*self.system.Natoms]).T, np.dot(H, dx.reshape([3*self.system.Natoms])))
        curves = [
            [q, strain, 'b-', '_nolegend_'],
            [q, energy, 'g-', '_nolegend_'],
        ]
        ax1 = pp.subplot(2*plc[0], plc[1], 2*plc[2]*plc[1]+1)
        add_plot(ax1, [curves[0]], title='Strain Cost', xunit=ic.qunit, yunit='1.0', xlabel=None, ylabel='Strain [-]')
        ax2 = pp.subplot(2*plc[0], plc[1], (2*plc[2]+1)*plc[1]+1, sharex=ax1)
        add_plot(ax2, [curves[1]], title='Energy Cost', xunit=ic.qunit, yunit=eunit, xlabel='%s [%%s]' %(ic.name.split('/')[0]), ylabel='Energy [%s]')
    
    
    def plot_ic_constrained(self, plc, ic, trajectory):
        q = np.zeros(len(trajectory), float)
        for i, dx in enumerate(trajectory):
            q[i]      = ic.value(self.system.sample['coordinates'] + dx)
        curves = [[np.array(range(len(trajectory))), q , 'b-', '_nolegend_']]
        ax = pp.subplot(plc[0], plc[1], plc[2])
        add_plot(ax, curves, title='Constrained IC - %s' %ic.name, xunit='1.0', yunit=ic.qunit, xlabel='Step [-]', ylabel='%s [%%s]' %(ic.name.split('/')[0]))
    
    
    def plot_ic_other_ics(self, plc, ic, other_ics, trajectory, title='Other ICs', ylabel='IC [%%s]', yunit='au'):
        styles = ['b-', 'r-', 'g-', 'c-', 'm-', 'y-', 'b--', 'r--', 'g--', 'c--', 'm--', 'y--']*100
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
    if xlabel is not None:
        if '%s' in xlabel: xlabel = xlabel %xunit
        ax.set_xlabel(xlabel, fontsize=20)
    if ylabel is not None:
        if '%s' in ylabel: ylabel = ylabel %yunit
        ax.set_ylabel(ylabel, fontsize=20)
    ax.grid()
    ax.set_title(title, fontsize=24)
    ax.legend(loc='best')
