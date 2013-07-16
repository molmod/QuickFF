#! /usr/bin/env python

from quickff.fftable import DataArray, FFTable
from quickff.perturbation import RelaxedGeoPertTheory
from quickff.cost import HessianFCCost

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
        self.pert_theory = RelaxedGeoPertTheory(system, model)
        self.cost = HessianFCCost(system, model)

    def estimate_from_pt(self, skip_dihedrals=True):
        '''
            First Step: calculate harmonic force field parameters from
            perturbation trajectories.

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
                k, q0 = self.pert_theory.estimate(ic)
                ks.append(k)
                q0s.append(q0)
            ff.add(icname, ks, q0s)
            descr = icname + ' '*(maxlength-len(icname))
            print '    %s   K = %s    q0 = %s' % (
                descr, ks.string(), q0s.string()
            )
        self.model.val.update_fftable(ff)
        return ff

    def refine_cost(self, fixed=None):
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
        self.model.val.determine_dihedral_potentials(self.system)
        print '\nEstimating all pars for bonds, bends and opdists\n'
        self.estimate_from_pt(skip_dihedrals=True)
        print '\nRefining force constants using a Hessian LSQ cost\n'
        self.refine_cost()
        return self.model.val.get_fftable()

    def plot_pt(self, icname):
        '''
            Generate and plot the perturbation trajectories for all ics with a
            name compatible with icname.
        '''
        for ic in self.system.ics[icname]:
            name = ic.name.replace('/', '-')
            trajectory = self.pert_theory.generate(ic)
            self.pert_theory.plot(ic, trajectory, 'energies-'+name+'.png')
            self.pert_theory.write(trajectory, 'trajectory-'+name+'.xyz')
