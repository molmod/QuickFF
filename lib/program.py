from quickff.fftable import DataArray, FFTable
from quickff.perturbation import RelaxedGeoPertTheory
from quickff.cost import HessianFCCost

import cPickle, os, sys, getpass, datetime

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
    def __init__(self, system, model, fns_traj=None):
        '''
            The central class to manage the entire program.

            **Arguments**

            system
                An instance of the System class and contains all the
                system information

            model
                An instance of the Model class and contains all the info
                to define the total PES and its electrostatic contribution.
            
            **Optional Arguments**
            
            fns_traj
                A file name to store the perturbation trajectories to. The
                trajectories are stored after Pickling.
        '''
        self.system = system
        self.model = model
        self.pert_theory = RelaxedGeoPertTheory(system, model)
        self.cost = HessianFCCost(system, model)
        self.fns_traj = fns_traj
        
    def generate_trajectories(self, skip_dihedrals=True):
        '''
            Generate a perturbation trajectory for all ics (dihedrals can be
            excluded and store the coordinates in a dictionary.

            **Optional Arguments**

            skip_dihedrals
                If set to True, the dihedral ff parameters will not
                be calculated.
        '''
        maxlength = max([len(icname) for icname in self.model.val.pot.terms.keys()]) + 2
        #Check if a filename with trajectories is given. If the file exists,
        #read it and return the trajectories
        if self.fns_traj is not None:
            if os.path.isfile(self.fns_traj):
                with open(self.fns_traj,'r') as f:
                    trajectories = cPickle.load(f)
                return trajectories 
        #Generate trajectories from scratch
        trajectories = {}
        for icname in sorted(self.model.val.pot.terms.keys()):
            ics = self.system.ics[icname]
            if skip_dihedrals and icname.startswith('dihed'):
                continue
            for i_ics, ic in enumerate(ics):
                sys.stdout.write('\r    %s Processing %2i/%i' %(
                    icname+' '*(maxlength-len(icname)), i_ics+1, len(ics)
                ))
                sys.stdout.flush()
                trajectories[ic.name] = self.pert_theory.generate(ic)
            print ''
        #Check if we need to write the generated trajectories to a file
        if self.fns_traj is not None:
            with open(self.fns_traj,'w') as f:
                cPickle.dump(trajectories,f)
        return trajectories

    def estimate_from_pt(self, trajectories, skip_dihedrals=True):
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
                k, q0 = self.pert_theory.estimate(ic,trajectories[ic.name])
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
        print '\nDetermine the coordinates of the perturbation trajectories\n'
        self.trajectories = self.generate_trajectories()
        print '\nEstimating all pars for bonds, bends and opdists\n'
        self.estimate_from_pt(self.trajectories)
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
