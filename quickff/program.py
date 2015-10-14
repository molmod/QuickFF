# -*- coding: utf-8 -*-
# QuickFF is a code to quickly derive accurate force fields from ab initio input.
# Copyright (C) 2012 - 2015 Louis Vanduyfhuys <Louis.Vanduyfhuys@UGent.be>
# Steven Vandenbrande <Steven.Vandenbrande@UGent.be>,
# Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center for Molecular Modeling
# (CMM), Ghent University, Ghent, Belgium; all rights reserved unless otherwise
# stated.
#
# This file is part of QuickFF.
#
# QuickFF is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# QuickFF is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
#--

'''Class to combine all steps of the methodology in convenient functions.
'''


import numpy as np
import cPickle, os, sys, getpass, datetime, numpy, scipy

from molmod.units import deg
from molmod import parse_unit
from yaff import DeltaList, Harmonic, Cosine, ForceField, ForcePartValence, Chebychev1

from quickff.fftable import DataArray, FFTable
from quickff.perturbation import RelaxedGeoPertTheory
from quickff.cost import HessianFCCost
from quickff.paracontext import *
from quickff.iclist import ICList
from quickff.tools import get_vterm

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

    Welcome to QuickFF 1.0 - a Python package to quickly derive force fields from ab initio input data

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
    def __init__(self, system, refdata, fn_traj=None, skip_ics=[], refineq=False):
        '''
            **Arguments**

            system
                An instance of the Yaff System class containing all the
                system information.

            refdata
                An instance of the ReferenceData class containing all
                information about the reference PES.

            **Optional Arguments**

            fn_traj
                A file name to store the perturbation trajectories to or to
                read the trajectories from if the file exists. The trajectories
                are stored/read after Pickling.
        '''
        self.system = system
        self.refdata = refdata
        self.skip_ics = skip_ics
        self.make_iclist(skip_ics)
        self.refineq = refineq
        self.pert_theory = RelaxedGeoPertTheory(self.system, self.refdata, self.iclist)
        self.cost = HessianFCCost(self.system, self.refdata, self.iclist)
        self.fn_traj = fn_traj
        self._check()

    def _check(self):
        assert len(self.iclist.ics)>1, "Strain function not defined for systems with only one IC"
        assert np.sum(np.array(self.refdata.pbc))==self.system.cell.nvec
        if self.system.bonds is None:
            self.system.detect_bonds()
        assert self.system.ffatypes is not None

    def make_iclist(self, skip_ics):
        self.dlist = DeltaList(self.system)
        self.iclist = ICList(self.dlist)
        # Add all internal coordinates
        for bond in self.system.iter_bonds():
            self.iclist.add_ic(bond, 'bond', skip_ics)
        for angle in self.system.iter_angles():
            self.iclist.add_ic(angle, 'angle', skip_ics)
        for dihed in self.system.iter_dihedrals():
            self.iclist.add_ic(dihed, 'dihed', skip_ics)
        for oop in self.system.iter_oops():
            self.iclist.add_ic(oop, 'opdist', skip_ics)
        self.iclist.finalize_icnames()
        self.system.pos[:] = self.refdata.coords
        self.iclist.dlist.forward()
        self.iclist.forward()

    def generate_trajectories(self, skip_dihedrals=True, verbose=True):
        '''
            Generate a perturbation trajectory for all ics (dihedrals can be
            excluded and store the coordinates in a dictionary.

            **Optional Arguments**

            skip_dihedrals
                If set to True, the dihedral ff parameters will not
                be calculated.
        '''
        maxlength = max([len(icname) for icname in self.iclist.icnames]) + 2
        #Check if a filename with trajectories is given. If the file exists,
        #read it and return the trajectories
        if self.fn_traj is not None:
            if os.path.isfile(self.fn_traj):
                with open(self.fn_traj,'r') as f:
                    trajectories = cPickle.load(f)
                return trajectories
        #Generate trajectories from scratch
        trajectories = {}
        all_ics = []
        for iic in xrange(self.iclist.nic):
            if skip_dihedrals and self.iclist.icnames[self.iclist.icname_ids[iic]].startswith('dihed'):
                continue
            all_ics.append(iic)
        print "Generating %d trajectories" % (len(all_ics))
        results = paracontext.map(self.pert_theory.generate, all_ics)
        for i, iic in enumerate(all_ics):
            trajectories[iic] = results[i]
        #Check if we need to write the generated trajectories to a file
        if self.fn_traj is not None:
            with open(self.fn_traj,'w') as f:
                cPickle.dump(trajectories,f)
        return trajectories

    def estimate_from_pt(self, trajectories, skip_dihedrals=True, verbose=True, ff=None):
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

            ff
                FFTable instance to which the parameters found here will be added
                If not given, a new FFTable will be created.
        '''
        if ff is None: ff = FFTable()
        maxlength = max([len(icname) for icname in self.iclist.icnames]) + 2
        for iname, icname in enumerate(self.iclist.icnames):
            iics = np.where(self.iclist.icname_ids == iname)[0]
            if skip_dihedrals and icname.startswith('dihed'):
                continue
            ks  = DataArray(unit=self.iclist.kunits[iname])
            q0s = DataArray(unit=self.iclist.qunits[iname])
            for iic in iics:
                k, q0 = self.pert_theory.estimate(iic,trajectories[iic])
                ks.append(k)
                q0s.append(q0)
            ff.add(icname, ks, q0s)
            descr = icname + ' '*(maxlength-len(icname))
        if verbose:
            ff.print_screen()
        return ff

    def estimate_linear_model(self, skip_dihedrals=True, verbose=True, ff=None):
        '''
            Second step of force-field development.

            **Optional Arguments**

            skip_dihedrals
                If set to True, the dihedral ff parameters will not
                be calculated.

            ff
                FFTable instance to which the parameters found here will be added
                If not given, a new FFTable will be created.
        '''
        if ff is None: ff = FFTable()
        maxlength = max([len(icname) for icname in self.iclist.icnames]) + 2
        ndof = 3*self.system.natom
        # Construct matrix containing gradients of ics wrt cartesian
        # coordinates as columns
        grads = np.zeros((self.iclist.nic,self.system.natom,3), float)
        qstars = np.zeros((self.iclist.nic,), float)
        for iic in xrange(self.iclist.nic):
            q = ForcePartValence(self.system)
            q.add_term(Chebychev1(-2.0,self.iclist.ics[iic]))
            q_ff = ForceField(self.system,[q])
            q_ff.update_pos(self.refdata.coords)
            qstars[iic] =  q_ff.compute(grads[iic,:]) + 1.0
            if np.amax(np.abs(grads[iic]))>1e5:
                #TODO Usually happens for dihed with phi=180,
                #check for other cases where this might happen
                if self.iclist.icnames[self.iclist.icname_ids[iic]].startswith('dihed') and skip_dihedrals:
                    grads[iic] = 0.0
                else:
                    raise UserWarning, "Singularity in gradient for ic %d: %30s, q0 = %f, gmax = %f" % \
                        (iic,self.iclist.icnames[self.iclist.icname_ids[iic]],qstars[iic],np.amax(np.abs(grads[iic])))
            if np.all(np.abs(grads[iic])<1e-10):
                if self.iclist.icnames[self.iclist.icname_ids[iic]].startswith('dihed') and\
                     skip_dihedrals: pass
                else: raise UserWarning, "Gradient for ic %d: %30s is zero, the linear model will not work..." \
                     % (iic,self.iclist.icnames[self.iclist.icname_ids[iic]])
        grads = grads.reshape((self.iclist.nic, self.system.natom*3))
        #TODO Figure out difference between phess and hess on final force field
        hess = self.refdata.phess.reshape((ndof, ndof))
        if self.refdata.ncff is not None: hess -= self.refdata.ff_hessian
        grad = self.refdata.grad
        if self.refdata.ncff is not None: grad -= self.refdata.ff_grad
        grad = grad.reshape((ndof,))
        # Estimate force constant and rest value for all ics
        for iname, icname in enumerate(self.iclist.icnames):
            iics = np.where(self.iclist.icname_ids == iname)[0]
            if skip_dihedrals and icname.startswith('dihed'):
                continue
            ks  = DataArray(unit=self.iclist.kunits[iname])
            q0s = DataArray(unit=self.iclist.qunits[iname])
            rhs = np.zeros(self.iclist.nic)
            for iic in iics:
                # Determine k and q0 by solving an overdetermined system of
                # equations. This will go haywire for ics with grad zero...
                rhs[:] = 0.0
                rhs[iic] = 1.0
                sol = np.linalg.lstsq(grads, rhs, rcond=1e-9)
                x = sol[0]
                k = np.dot(x.transpose(),np.dot(hess,x))
                q0 = qstars[iic] - np.dot(grad,x)/k
                ks.append(k)
                q0s.append(q0)
            ff.add(icname, ks, q0s)
            descr = icname + ' '*(maxlength-len(icname))
        if verbose:
            ff.print_screen()
        return ff

    def refine_cost(self, fftab, verbose=True):
        '''
            Second step of force field development: refine the force constants
            using a Hessian least squares cost function.

            ff:
                FFTable object containing current force-field parameters
        '''
        self.cost.estimate(fftab)
        if verbose:
            fftab.print_screen()
        return fftab

    def refine_rvs(self, fftab, verbose=True):
        '''
            Third step of force field development: refine the rest values by
            revisiting the perturbation trajectories of the first step with the
            force constants of the second step.
        '''
        # Construct a valence force field based on the current parameters
        val = ForcePartValence(self.system)
        for iic in xrange(self.iclist.nic):
            icname = self.iclist.icnames[self.iclist.icname_ids[iic]]
            if not icname in fftab.pars.keys(): continue
            if icname.startswith('dihed'):
                pars = [fftab.pars[icname]['m'].mean,fftab.pars[icname]['k'].mean,fftab.pars[icname]['q0'].mean]
                indexes = self.iclist.ics[iic].index_pairs
                term = get_vterm(pars, [indexes[0][1],indexes[0][0],indexes[2][0],indexes[2][1]])
                if term is None:
                    term = Cosine(pars[0], pars[1], pars[2], self.iclist.ics[iic])
            else:
                term = Harmonic(fftab.pars[icname]['k'].mean, fftab.pars[icname]['q0'].mean, self.iclist.ics[iic])
            val.add_term(term)
        ff = ForceField(self.system,[val])
        for iname, icname in enumerate(self.iclist.icnames):
            iics = np.where(self.iclist.icname_ids == iname)[0]
            if icname.startswith('dihed'): continue
            ks  = DataArray(unit=self.iclist.kunits[iname])
            q0s = DataArray(unit=self.iclist.qunits[iname])
            for iic in iics:
                k, q0 = self.pert_theory.refineq(iic,self.trajectories[iic],(fftab.pars[icname]['k'].mean,fftab.pars[icname]['q0'].mean), ff)
                ks.append(k)
                q0s.append(q0)
            fftab.pars[icname]['q0'].data[:] = q0s.data[:]
            fftab.pars[icname]['q0'].update_statistics()
        if verbose:
            fftab.print_screen()
        return fftab

    def run(self,linear_model=False):
        '''
            Run all steps of the QuickFF methodology to derive a covalent
            force field. This method returns an instance of the class
            :class:`quickff.fftable.FFTable`, which contains all force field
            parameters.
        '''
        print header
        print sysinfo()
        print '~'*120+'\n'
        print '\nDetermine dihedral potentials\n'
        fftab, deleted_diheds = self.determine_dihedral_potentials()
        if not linear_model:
            print '\nDetermine the coordinates of the perturbation trajectories\n'
            self.trajectories = self.generate_trajectories()
            print '\nEstimating all pars for bonds, bends and opdists\n'
            fftab = self.estimate_from_pt(self.trajectories, ff=fftab)
        else:
            print '\nEstimating all pars for bonds, bends and opdists\n'
            fftab = self.estimate_linear_model(ff=fftab)
        print '\nRefining force constants using a Hessian LSQ cost\n'
        fftab = self.refine_cost(fftab)
        if self.refineq:
            print '\nRefining rest values\n'
            fftab = self.refine_rvs(fftab)
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
            print '\nDetermine the coordinates of the perturbation trajectories\n'
        #Reading/generating trajectories
        trajectories = {}
        if self.fn_traj is not None:
            if os.path.isfile(self.fn_traj):
                with open(self.fn_traj,'r') as f:
                    trajectories = cPickle.load(f)
        iname = np.where(self.iclist.icnames==icname)[0][0]
        iics = np.where(self.iclist.icname_ids == iname)[0]
        print iics
        for i,iic in enumerate(iics):
            if iic in trajectories.keys():
                #already read
                if verbose:
                    print '    %s Read %2i/%i from %s' %(
                        icname, i+1, iics.shape[0], self.fn_traj
                    )
            else:
                #generating
                if verbose:
                    sys.stdout.write('    %s Generating %2i/%i' %(
                        icname, i+1, iics.shape[0]
                    ))
                    sys.stdout.flush()
                try:
                    trajectories[iic] = self.pert_theory.generate(iic, start=start, end=end, steps=51)
                    print ''
                except KeyboardInterrupt:
                    if verbose:
                        sys.stdout.write(' INTERRUPTED\n')
                        sys.stdout.flush()
        #Writing trajectories
        if self.fn_traj is not None:
            with open(self.fn_traj,'w') as f:
                cPickle.dump(trajectories, f)
        #Plotting/writing output
        for i,iic in enumerate(iics):
            if iic in trajectories.keys():
                name = "%s_%05d" % (icname.replace('/', '-'),i)
                self.pert_theory.plot(iic, trajectories[iic], 'energies-'+name+'.pdf')
                self.pert_theory.write(trajectories[iic], 'trajectory-'+name+'.xyz')
        #Logging
        if verbose:
            print ''
            print '\n'+'~'*120+'\n'
            print 'Time:           ' + datetime.datetime.now().isoformat().replace('T', ' ') + '\n'
            print footer

    def determine_dihedral_potentials(self, marge2=15*deg, marge3=15*deg, verbose=True):
        '''
            Determine the potential of every dihedral based on the values of
            the dihedral angles in the geometry. First try if a cosine potential
            of the form 0.5*K*[1 - cos(m(psi-psi0))] works well with m=2,3 and
            psi0 = 0,pi/m. If this doesn't work, raise a warning and ignore
            dihedral.
        '''
        ff = FFTable()
        maxlength = max([len(icname) for icname in self.iclist.icnames]) + 2
        deleted_diheds = []
        def determine_mrv(m, psi, verbose=False):
            per = 360*deg/m
            val = psi%per
            if verbose: print per, val
            if (val>=0 and val<=per/6.0) or (val>=5*per/6.0 and val<=per):
                return m, 0.0
            elif val>=2*per/6.0 and val<=4*per/6.0:
                return m, per/2.0
            else:
                return -1, 0.0
        self.system.pos[:] = self.refdata.coords
        self.iclist.dlist.forward()
        self.iclist.forward()
        for icname in self.iclist.icnames:
            if not icname.startswith('dihed'):
                continue
            ms = []
            rvs = []
            ks = []
            descr = icname + ' '*(maxlength-len(icname))
            iics = np.where(self.iclist.icname_ids == np.where(self.iclist.icnames==icname)[0][0])[0]
            for iic in iics:
                psi0 = self.iclist.ictab[iic]['value']
                i = self.iclist.ics[iic].index_pairs[0][0]
                j = self.iclist.ics[iic].index_pairs[1][1]
                n1 = len(self.system.neighs1[i])
                n2 = len(self.system.neighs1[j])

                if set([n1,n2])==set([4,4]):
                    m, rv = determine_mrv(3, psi0)
                elif set([n1,n2])==set([3,4]):
                    m, rv = determine_mrv(6, psi0)
                elif set([n1,n2])==set([2,4]):
                    m, rv = determine_mrv(3, psi0)
                elif set([n1,n2])==set([3,3]):
                    m, rv = determine_mrv(2, psi0)
                elif set([n1,n2])==set([2,3]):
                    m, rv = determine_mrv(2, psi0)#, verbose=(icname == "dihed/Al.O_CA.C_CA.C_PC"))
                elif set([n1,n2])==set([2,2]):
                    m, rv = determine_mrv(1, psi0)
                else:
                    m, rv = -1, 0.0
                ms.append(m)
                rvs.append(rv)
                ks.append(0.0)
            m = DataArray(ms, unit='au')
            rv = DataArray(rvs, unit='deg')
            k = DataArray(ks, unit='kjmol')

            if m.mean == -1 or m.std > 0.0 or rv.std > 0.0:
                if verbose:
                    print '    %s   WARNING: ' % descr +\
                          'could not determine clear trent in dihedral patterns, ' +\
                          'dihedral is ignored in force field !!!'
                    #print m.data, rv.data
                deleted_diheds.append(icname)
            else:
                if verbose:
                    print '    %s   0.5*K*[1 - cos(%i(psi - %5.1f))]' % (
                        descr, m.mean, rv.mean/deg
                    )
                for iic in enumerate(iics):
                    pass
                ff.add( icname, k, rv,m=m)
        return ff, deleted_diheds
