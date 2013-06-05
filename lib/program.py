#! /usr/bin/env python

from molmod.units import *

import numpy as np

from fftable import DataArray, FFTable
from perturbation import RelaxedGeoPertTheory
from cost import HessianFCCost

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

    def estimate_pt(self, skip_dihedrals=True):
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

    def refine_cost(self, tol=1e-9):
        '''
            Second Step: refine the force constants using a Hessian least
                         squares cost function.

            **Optional Arguments**

            tol
                a float defining the convergence tolerence on the force
                constants.
        '''
        fcs = self.cost.estimate(tol)
        self.model.val.update_fcs(fcs)
        ff = self.model.val.get_fftable()
        return ff

    def run(self):
        print 'System information:'
        print
        self.system.print_atom_info()
        print
        print 'Estimating ff pars from perturbation trajectories'
        print
        fftab1 = self.estimate_pt()
        print
        print 'Refining force constants using a Hessian LSQ cost'
        print
        fftab2 = self.refine_cost()
        return fftab2

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
