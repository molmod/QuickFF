from quickff.fftable import DataArray, FFTable
from quickff.perturbation import RelaxedGeoPertTheory
from quickff.cost import HessianFCCost

import os, sys, getpass, datetime

__all__ = ['Program']

header = r'''
____________/\\\________________________________________________________/\\\\\\\\\\\\\\\__/\\\\\\\\\\\\\\\______________
__________/\\\\/\\\\_______________________________________/\\\_________\/\\\///////////__\/\\\///////////______________
_________/\\\//\////\\\_________________/\\\_______________\/\\\_________\/\\\_____________\/\\\________________________
_________/\\\______\//\\\__/\\\____/\\\_\///______/\\\\\\\\_\/\\\\\\\\____\/\\\\\\\\\\\_____\/\\\\\\\\\\\_______________
_________\//\\\______/\\\__\/\\\___\/\\\__/\\\___/\\\//////__\/\\\////\\\__\/\\\///////______\/\\\///////_______________
___________\///\\\\/\\\\/___\/\\\___\/\\\_\/\\\__/\\\_________\/\\\\\\\\/___\/\\\_____________\/\\\_____________________
______________\////\\\//_____\/\\\___\/\\\_\/\\\_\//\\\________\/\\\///\\\___\/\\\_____________\/\\\____________________
__________________\///\\\\\\__\//\\\\\\\\\__\/\\\__\///\\\\\\\\_\/\\\_\///\\\_\/\\\_____________\/\\\___________________
_____________________\//////____\/////////___\///_____\////////__\///____\///__\///______________\///___________________

    Welcom to QuickFF 1.0 - a Python package to quickly derive force fields from ab initio input data

                                                Written by
                    Louis Vanduyfhuys(1)*, Steven Vandenbrande(1) and Toon Verstraelen(1)

(1) Center for Molecular Modeling, Ghent University Belgium.
* mailto: Louis.Vanduyfhuys@UGent.be
'''

footer = r'''
__/\\\__________________________________________________________________________________________________________/\\\____
  \ \\\                                                                                                         \ \\\
   \ \\\                        End of file. Thanks for using QuickFF! Come back soon!!                          \ \\\
____\///__________________________________________________________________________________________________________\///__
'''

def sysinfo():
    info  = '\nUser:           ' + getpass.getuser() + '\n'
    info += 'Machine info:   ' + ' '.join(os.uname()) + '\n'
    info += 'Time:           ' + datetime.datetime.now().isoformat().replace('T', ' ') + '\n'
    info += 'Python version: ' + sys.version.replace('\n', '') + '\n'
    info += 'Current Dir:    ' + os.getcwd() + '\n'
    info += 'Command line:   ' + ' '.join(sys.argv) + '\n'
    return info

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
            Second Step of force field development: calculate harmonic force field 
            parameters for every internal coordinate separately from perturbation
            trajectories.

            **Optional Arguments**

            skip_dihedrals
                If set to True, the dihedral ff parameters will not
                be calculated.
        '''
        ff = FFTable()
        maxlength = max([len(icname) for icname in self.model.val.pot.terms.keys()]) + 2
        for icname in sorted(self.model.val.pot.terms.keys()):
            ics = self.system.ics[icname]
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

    def refine_cost(self):
        '''
            Second step of force field development: refine the force constants 
            using a Hessian least squares cost function.
        '''
        fcs = self.cost.estimate()
        self.model.val.update_fcs(fcs)
        self.model.val.get_fftable().print_screen()

    def run(self):
        print header
        print sysinfo()
        print '~'*120+'\n'
        print 'System information:\n'
        self.system.print_atom_info()
        print '\nModel information:\n'
        self.model.print_info()
        print '\nDetermine dihedral potentials\n'
        self.model.val.determine_dihedral_potentials(self.system)
        print '\nEstimating all pars for bonds, bends and opdists\n'
        self.estimate_from_pt()
        print '\nRefining force constants using a Hessian LSQ cost\n'
        self.refine_cost()
        print '\n'+'~'*120+'\n'
        print 'Time:           ' + datetime.datetime.now().isoformat().replace('T', ' ') + '\n'
        print footer
        return self.model.val.get_fftable()

    def plot_pt(self, icname):
        '''
            Generate and plot the perturbation trajectories for all ics with a
            name compatible with icname.
        '''
        for ic in self.system.ics[icname]:
            name = ic.name.replace('/', '-')
            trajectory = self.pert_theory.generate(ic, steps=51)
            self.pert_theory.plot(ic, trajectory, 'energies-'+name+'.pdf')
            self.pert_theory.write(trajectory, 'trajectory-'+name+'.xyz')
