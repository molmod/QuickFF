#! /usr/bin/env python

from molmod.units import *
from molmod.ic import dihed_cos, dihed_angle

import numpy as np

from fftable import DataArray, FFTable
from perturbation import RelaxedGeoPertTheory
from cost import HessianFCCost
from terms import HarmonicTerm, CosineTerm

__all__ = ['Program']

class Program(object):
    def __init__(self, system, model):
        '''
            The central class to manage the entire program.

            **Arguments**

            system
                An instance of the System class and contains all the
                system information

            model
                An instance of the Model class and contains all the info
                to define the total PES and its electrostatic contribution.
        '''
        self.system = system
        self.model = model
        self.pt = RelaxedGeoPertTheory(system, model)
        self.cost = HessianFCCost(self.system, self.model)

    def estimate_from_pt(self, skip_dihedrals=True):
        '''
            First Step: calculate harmonic force field parameters from perturbation
            trajectories.

            **Optional Arguments**

            skip_dihedrals
                If set to True, the dihedral ff parameters will not
                be calculated.
        '''
        ff = FFTable()
        maxlength = max([len(icname) for icname in self.system.ics.keys()]) + 2
        for icname, ics in sorted(self.system.ics.iteritems()):
            if skip_dihedrals and icname.startswith('dihed'):
                continue
            ks  = DataArray(unit=ics[0].kunit)
            q0s = DataArray(unit=ics[0].qunit)
            for ic in ics:
                k, q0 = self.pt.estimate(ic)
                ks.append(k)
                q0s.append(q0)
            ff.add(icname, ks, q0s)
            descr = icname + ' '*(maxlength-len(icname))
            print '    %s   K = %s    q0 = %s' %(descr, ks.string(), q0s.string())
        self.model.val.update_fftable(ff)
        return ff

    def determine_dihedral_potentials(self, marge2=15*deg, marge3=15*deg):
        '''
            Determine the potential of every dihedral based on the values of
            the dihedral angles in the geometry. First try if a cosine potential
            of the form 0.5*K*[1 - cos(m(psi-psi0))] works well with m=2,3 and
            psi0 = 0,pi/m. Else use a term harmonic in cos(psi).
        '''
        maxlength = max([len(icname) for icname in self.system.ics.keys()]) + 2
        deleted_diheds = []
        for icname, ics in self.system.ics.iteritems():
            if not icname.startswith('dihed'):
                continue
            ms = []
            rvs = []
            descr = icname + ' '*(maxlength-len(icname))
            for ic in ics:
                psi0 = abs(ic.value(self.system.ref.coords))
                n1 = len(self.system.nlist[ic.indexes[1]])
                n2 = len(self.system.nlist[ic.indexes[2]])
                if psi0>=0 and psi0<=max([marge2,marge3]):
                    if 4 in [n1, n2]: #use m=3 if at least one of the central atoms has 4 neighbors
                        ms.append(3)
                        rvs.append(0.0)
                    elif 3 in [n1, n2]: #use m=2 if at least one of the central atoms has 3 neighbors
                        ms.append(2)
                        rvs.append(0.0)
                    else: #use m=1 else (system has exactly 4 atoms)
                        ms.append(1)
                        rvs.append(0.0)
                elif psi0>=60*deg-marge3 and psi0<=60*deg+marge3:
                    ms.append(3)
                    rvs.append(60.0*deg)
                elif psi0>=90*deg-marge2 and psi0<=90*deg+marge2:
                    ms.append(2)
                    rvs.append(90.0*deg)
                elif psi0>=120*deg-marge3 and psi0<=120*deg+marge3:
                    ms.append(3)
                    rvs.append(0.0*deg)
                elif psi0>=180*deg-marge2 and psi0<=180*deg+marge2:
                    if 4 in [n1, n2]: #use m=3 if at least one of the central atoms has 4 neighbors
                        ms.append(3)
                        rvs.append(60.0*deg)
                    elif 3 in [n1, n2]: #use m=2 if at least one of the central atoms has 3 neighbors
                        ms.append(2)
                        rvs.append(0.0)
                    else: #use m=1 else (system has exactly 4 atoms)
                        ms.append(1)
                        rvs.append(180.0*deg)
                else:
                    ms.append(-1)
                    rvs.append(np.cos(psi0))
            m = DataArray(ms, unit='au')
            if m.mean==-1 or m.std>0.0:
                print '    %s   WARNING: could not determine clear trent in dihedral angles, dihedral is ignored in force field !!!' %descr
                deleted_diheds.append(icname)
            else:
                rv = DataArray(rvs, unit='deg')
                print '    %s   0.5*K*[1 - cos(%i(psi - psi0))]    psi0 = %s' %(descr, m.mean, rv.string())
                for i, ic in enumerate(ics):
                    ic.icf = dihed_angle
                    ic.qunit = 'deg'
                    self.model.val.vterms[icname][i] = CosineTerm(ic, self.system.ref.coords, 0.0, rv.mean, m.mean)
        for icname in deleted_diheds:
            del self.system.ics[icname]
            del self.model.val.vterms[icname]


    def refine_cost(self, fixed=[]):
        '''
            Second Step: refine the force constants using a Hessian least
                         squares cost function.
        '''
        fcs = self.cost.estimate(fixed=fixed)
        self.model.val.update_fcs(fcs)
        self.model.val.get_fftable().print_screen()

    def run(self):
        print 'System information:\n'
        self.system.print_atom_info()
        print '\nDetermine dihedral potentials\n'
        self.determine_dihedral_potentials()
        print '\nEstimating ff pars from perturbation trajectories for bonds, bends and opdists\n'
        self.estimate_from_pt(skip_dihedrals=True)
        print '\nRefining force constants using a Hessian LSQ cost\n'
        self.refine_cost(fixed=[])
        return self.model.val.get_fftable()

    def plot_pt(self, icname):
        '''
            Generate and plot the perturbation trajectories for all ics with a
            name compatible with icname.
        '''
        for ic in self.system.ics[icname]:
            name = ic.name.replace('/', '-')
            trajectory = self.pt.generate(ic)
            self.pt.plot(ic, trajectory, 'energies-'+name+'.png')
            self.pt.write(trajectory, 'trajectory-'+name+'.xyz')
