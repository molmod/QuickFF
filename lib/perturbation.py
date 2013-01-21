#! /usr/bin/env python

from molmod.units import *
from molmod.periodic import periodic as pt
from molmod.io.xyz import XYZWriter
import numpy as np, matplotlib.pyplot as pp
from scipy.optimize import minimize

from fftable import FFTable
from tools import fitpar
from evaluators import *

__all__ = ['PerturbationTheory', 'RelaxedGeometryPT', 'MinimalDeviationPT']


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
            trajectories = self.perturbation_trajectory(system, icname, start=0.9, end=1.1, steps=101)
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

    def perturbation_trajectory(self, system, icname, start=0.9, end=1.1, steps=101):
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
        trajectories = self.perturbation_trajectory(system, icname, start=0.9, end=1.1, steps=101)
        pp.clf()
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
    def __init__(self, coupling=None, free_depth=0, spring=10.0*kjmol/angstrom**2):
        """
            Coupling is a numpy array with the coefficients of the linear 
            combination describing the coupling.
            
            The following strings are also supported for coupling:
                ``symmetric``   A symmetric coupling 
                                i.e. coupling = np.ones(N)/np.sqrt(N)
                                with N the number of ics for the given icname

            If free_depth is None, the hessian remains unbiased for calculating 
            the perturbed geometry. If free_depth is e.g. equal to 3, all atoms 
            that are separated by more then 3 atoms from the atoms in the ic 
            under consideration are fixed with a spring to their original position
            (with strength defined in spring).
        """
        self.free_depth = free_depth
        self.spring = spring
        self.coupling = coupling
        self.description = 'Relaxed Geometry Perturbation Theory with \n' + \
                           '              coupling=%10s    free_depth=%1i    spring=%.3f kjmol/A^2' %(
            self.coupling, self.free_depth, self.spring/(kjmol/angstrom**2)
        )
    
    
    def perturbation_trajectory(self, system, icname, start=0.8, end=1.2, steps=101):
        #TODO: Method has to be updated so it calculates the complete trajectory at once accoring to given icrange.
        """
            Calculate the perturbation on the geometry when perturbing
            along <icname>. The perturbation trajectory is calculated with
            relaxed geometry and a value of icname equal to one.
        """
        raise NotImplementedError('Method has to be updated so it calculates the complete trajectory at once accoring to given icrange.')
        
        print '                %40s' %(icname)
        ics = system.ics[icname]
        if self.coupling=='symmetric':
            coupling = np.ones(len(ics), float)/len(ics)
        elif isinstance(self.coupling, np.ndarray):
            assert len(ics)==len(self.coupling)
            coupling = self.coupling.copy()
        elif self.coupling is None:
            coupling = None
        else:
            raise NotImplementedError('Unsupported coupling: %s' %(str(self.coupling)))
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
            if self.free_depth==0:
                ihess = system.totmodel.ihess.reshape([3*system.Natoms, 3*system.Natoms])
            else:
                indices = [i for i in xrange(system.Natoms) if np.linalg.norm(qgrad.reshape([system.Natoms, 3])[i])>1e-3]
                free_indices = system.get_neighbors(indices, depth=self.free_depth)
                ihess = system.totmodel.get_constrained_ihess(free_indices, spring=self.spring).reshape([3*system.Natoms, 3*system.Natoms])
            v = np.dot(ihess, qgrad)
            kl = np.dot(qgrad.T, np.dot(ihess, qgrad))*np.sqrt(len(qgrads))
            vs.append(v.reshape((-1,3))/kl)
        return vs


class MinimalDeviationPT(PerturbationTheory):
    def __init__(self, depth=-1, bond_err=1e-6*angstrom, bend_err=1e-4*deg, dihed_err=1e-2*deg):
        self.depth = depth
        self.bond_err = bond_err
        self.bend_err = bend_err
        self.dihed_err = dihed_err
        self.description = 'Minimal Deviation Perturbation Theory with \n' + \
                           '              depth=%1i    bond_error = %.3e A   bend_error = %.3e deg   dihed_error = %.3e deg' %(
            self.depth, self.bond_err/angstrom, self.bend_err/deg, self.dihed_err/deg
        )
    
    def calc_qdm(self, system, ic, coords=None):
        """
            Calculate the Internal Coordinate Deveation Matrix (QDM).
            If sandwiched between a geometry perturbation vector, this
            matrix represents the weighted sum of the deviations of the
            internal coordinates, except for the ic given in args, 
            from their non-perturbed values.
            
        """
        if coords is None: coords = system.sample['coordinates']
        neighbors = system.get_neighbors(ic.indexes, depth=self.depth)
        ics_t = []
        QDM = np.zeros([3*system.Natoms, 3*system.Natoms], float)
        for icname2 in system.icnames:
            for ic2 in system.ics[icname2]:
                if ic2.name==ic.name:
                    continue
                include = True
                for i in ic2.indexes:
                    if i not in neighbors:
                        include = False
                        break
                if include:
                    if icname2.split('/')[0] in ['bond', 'dist']: sigma = self.bond_err
                    elif icname2.split('/')[0] in ['bend', 'angle']: sigma = self.bend_err
                    elif icname2.split('/')[0] in ['dihedral', 'dihed', 'torsion']: sigma = self.dihed_err
                    else: raise ValueError('Invalid icname, recieved %s' %icname2)
                    QDM += np.outer(ic2.grad(coords), ic2.grad(coords))/(sigma**2)
                    ics_t.append(ic2)
        return QDM, ics_t
    
    
    def perturbation_trajectory(self, system, icname, start=0.9, end=1.1, steps=101, Ncycles=3):
        """
            Calculate the perturbation on the geometry when perturbing
            along <icname> with magnitudes given in icrange. 
            
              icrange = [start*ic0, end*ic0] with a total of <steps> steps
            
            For a given
            dq from icrange, the perturbed geometry has to have a value of
            the given ic equal to dq and a minimal value for the other ics.
        """
        qunit = parse_unit(system.units[icname]['q'])
        coords0 = system.sample['coordinates']
        ics = system.ics[icname]
        trajectories = [[[] for i in xrange(steps)] for j in xrange(len(ics))]
        for iic, ic in enumerate(ics):
            print 'PERTUR TRAJC: %s' %ic.name
            QDM, ics_t = self.calc_qdm(system, ic)
            qarray = ( start + (end-start)/(steps-1)*np.array(range(steps),float) )*ic.value(system.sample['coordinates'])
            for iq, q in enumerate(qarray):
                dx0 = np.zeros(3*system.Natoms, float)
                def chi(dx):
                   return 0.5*np.dot(dx.T, np.dot(QDM, dx))
                def dchi_dx(dx):
                    return np.dot(QDM, dx)
                constraints = ({
                    'type': 'eq',
                    'fun' : lambda dx: ic.value(coords0 + dx.reshape((-1, 3))) - q,
                    'jac' : lambda dx: ic.grad(coords0 + dx.reshape((-1, 3))),
                })
                result = minimize(chi, dx0, method='SLSQP', jac=dchi_dx, constraints=constraints, tol=1e-9)
                trajectories[iic][iq] = result.x.reshape((-1, 3))
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
