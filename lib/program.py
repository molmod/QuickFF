# -*- coding: utf-8 -*-
#QuickFF is a code to quickly derive accurate force fields from ab initio input.
#Copyright (C) 2012 - 2014 Louis Vanduyfhuys <Louis.Vanduyfhuys@UGent.be>
#Steven Vandenbrande <Steven.Vandenbrande@UGent.be>, 
#Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center for Molecular Modeling
#(CMM), Ghent University, Ghent, Belgium; all rights reserved unless otherwise
#stated.
#
#This file is part of QuickFF.
#
#QuickFF is free software; you can redistribute it and/or
#modify it under the terms of the GNU General Public License
#as published by the Free Software Foundation; either version 3
#of the License, or (at your option) any later version.
#
#QuickFF is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, see <http://www.gnu.org/licenses/>
#
#--

from quickff.fftable import DataArray, FFTable
from quickff.perturbation import RelaxedGeoPertTheory
from quickff.cost import HessianFCCost
from quickff.paracontext import *

import cPickle, os, sys, getpass, datetime, numpy, scipy

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
    info += 'Numpy version:  ' + numpy.__version__ + '\n'
    info += 'Scipy version:  ' + scipy.__version__ + '\n'
    info += 'Current Dir:    ' + os.getcwd() + '\n'
    info += 'Command line:   ' + ' '.join(sys.argv) + '\n'
    return info

class Program(object):
    '''
        The central class to manage the entire program.
    '''
    def __init__(self, system, model, fns_traj=None):
        '''
            **Arguments**

            system
                An instance of the System class and contains all the
                system information

            model
                An instance of the Model class and contains all the info
                to define the total PES and its electrostatic contribution.

            **Optional Arguments**

            fns_traj
                A file name to store the perturbation trajectories to or to
                read the trajectories from if the file exists. The trajectories
                are stored/read after Pickling.
        '''
        self.system = system
        self.model = model
        self.pert_theory = RelaxedGeoPertTheory(system, model)
        self.cost = HessianFCCost(system, model)
        self.fns_traj = fns_traj

    def generate_trajectories(self, skip_dihedrals=True, verbose=True):
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
        all_ics = []
        for icname in sorted(self.model.val.pot.terms.keys()):
            if skip_dihedrals and icname.startswith('dihed'):
                continue
            ics = self.system.ics[icname]
            for ic in ics:
                all_ics.append(ic)
            if verbose:
                print '    %s will generate %i trajectories' %(
                    icname+' '*(maxlength-len(icname)), len(ics)
                )
        results = paracontext.map(self.pert_theory.generate, all_ics)
        for i in xrange(len(all_ics)):
            trajectories[all_ics[i].name] = results[i]
        #Check if we need to write the generated trajectories to a file
        if self.fns_traj is not None:
            with open(self.fns_traj,'w') as f:
                cPickle.dump(trajectories,f)
        return trajectories

    def estimate_from_pt(self, trajectories, skip_dihedrals=True, verbose=True):
        '''
            Second Step of force field development: calculate harmonic force field
            parameters for every internal coordinate separately from perturbation
            trajectories.

            **Arguments**

            trajectories
                A dictionairy containing numpy arrays representing perturbation
                trajectories for each icname.

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
            if verbose:
                print '    %s   K = %s    q0 = %s' % (
                    descr, ks.string(), q0s.string()
                )
        self.model.val.update_fftable(ff)
        return ff

    def refine_cost(self, verbose=True):
        '''
            Second step of force field development: refine the force constants
            using a Hessian least squares cost function.
        '''
        fcs = self.cost.estimate()
        self.model.val.update_fcs(fcs)
        fftab = self.model.val.get_fftable()
        if verbose:
            fftab.print_screen()
        return fftab

    def run(self):
        '''
            Run all steps of the QuickFF methodology to derive a covalent
            force field. This method returns an instance of the class
            :class:`quickff.fftable.FFTable`, which contains all force field
            parameters.
        '''
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
        fftab = self.estimate_from_pt(self.trajectories)
        print '\nRefining force constants using a Hessian LSQ cost\n'
        fftab = self.refine_cost()
        print '\n'+'~'*120+'\n'
        print 'Time:           ' + datetime.datetime.now().isoformat().replace('T', ' ') + '\n'
        print footer
        return fftab

    def plot_pt(self, icname, start=None, end=None, steps=51, verbose=True):
        '''
            Generate the perturbation trajectories and plot the energy
            contributions along the trajectory for all ics with a name
            compatible with icname.

            **Arguments**

            icname
                A string describing for which ics the perturbation trajectories
                should be generated.
        '''
        #Logging
        if verbose:
            print header
            print sysinfo()
            print '~'*120+'\n'
            print 'System information:\n'
            self.system.print_atom_info()
            print '\nModel information:\n'
            self.model.print_info()
            print '\nDetermine the coordinates of the perturbation trajectories\n'
        #Reading/generating trajectories
        trajectories = {}
        if self.fns_traj is not None:
            if os.path.isfile(self.fns_traj):
                with open(self.fns_traj,'r') as f:
                    trajectories = cPickle.load(f)
        for i, ic in enumerate(self.system.ics[icname]):
            if ic.name in trajectories.keys():
                #already read
                if verbose:
                    print '    %s Read %2i/%i from %s' %(
                        icname, i+1, len(self.system.ics[icname]), self.fns_traj
                    )
            else:
                #generating
                if verbose:
                    sys.stdout.write('    %s Generating %2i/%i' %(
                        icname, i+1, len(self.system.ics[icname])
                    ))
                    sys.stdout.flush()
                try:
                    trajectories[ic.name] = self.pert_theory.generate(ic, start=start, end=end, steps=51)
                    print ''
                except KeyboardInterrupt:
                    if verbose:
                        sys.stdout.write(' INTERRUPTED\n')
                        sys.stdout.flush()
        #Writing trajectories
        if self.fns_traj is not None:
            with open(self.fns_traj,'w') as f:
                cPickle.dump(trajectories, f)
        #Plotting/writing output
        for ic in self.system.ics[icname]:
            if ic.name in trajectories.keys():
                name = ic.name.replace('/', '-')
                self.pert_theory.plot(ic, trajectories[ic.name], 'energies-'+name+'.pdf')
                self.pert_theory.write(trajectories[ic.name], 'trajectory-'+name+'.xyz')
        #Logging
        if verbose:
            print ''
            print '\n'+'~'*120+'\n'
            print 'Time:           ' + datetime.datetime.now().isoformat().replace('T', ' ') + '\n'
            print footer
