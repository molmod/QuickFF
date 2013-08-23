#! /usr/bin/env python

from molmod.units import deg
from molmod.ic import bond_length, dihed_angle

import numpy as np

from quickff.ic import IC
from quickff.terms import HarmonicTerm, CosineTerm
from quickff.fftable import DataArray, FFTable

__all__ = [
    'Model', 'ZeroPart', 'HarmonicPart', 'CoulombPart', 'ValencePart'
]


class Model(object):
    def __init__(self, total, val, ei):
        '''
           A class defining the ab initio total energy of the system,
           the force field electrostatic contribution and the
           force field valence terms.

           **Arguments**

           total
                A model for the ab initio total energy, should be an instance
                of HarmonicPart

           ei
                A model for the force field electrostatic energy, should be
                an instance of HarmonicPart or CoulombPart

           val
                an instance of the ValencePart class containing all details
                of the force field valence contributions.
        '''
        self.total = total
        self.ei = ei
        self.val = val

    @classmethod
    def from_system(cls, system, project=True, eirule=0, eikind='Harmonic'):
        '''
           **Arguments**

           system
                An instance of the System class containing all the
                information of the system.

           **Optional Arguments**

           project
                If True, project the translational and rotational
                degrees of freedom out of the hessian.

           eirule
                an integer defining the exculsion rule of the
                electrostatic interactions. Can be -1,0,1,2,3

           eikind
                a string defining the model kind of the electrostatic
                interactions. Can be 'Harmonic' or 'Coulomb'
        '''
        #Total model
        if project:
            total = HarmonicPart('Total Harmonic', system.ref.coords, 0.0, system.ref.grad, system.ref.phess)
        else:
            total = HarmonicPart('Total Harmonic', system.ref.coords, 0.0, system.ref.grad, system.ref.hess)
        #EI model
        epairs = _get_excluded_pairs(system, eirule)
        if eirule == -1:
            ei = ZeroPart('EI Zero')
        elif eikind == 'Harmonic':
            grad_ei, hess_ei = _calc_ei_grad_hess(system, epairs)
            ei = HarmonicPart('EI Harmonic', system.ref.coords, 0.0, grad_ei, hess_ei)
        elif eikind == 'Coulomb':
            ei = CoulombPart('EI Coulomb', system.ref.coords, system.charges, epairs)
        #Valence terms
        val = ValencePart(system)
        return cls(total, val, ei)


class ZeroPart(object):
    def __init__(self, name):
        self.name = name

    def calc_energy(self, coords):
        return 0.0

    def calc_gradient(self, coords):
        return np.zeros([len(coords), 3], float)

    def calc_hessian(self, coords):
        return np.zeros([len(coords), 3, len(coords), 3], float)


class HarmonicPart(object):
    def __init__(self, name, coords0, energy0, gradient, hessian):
        self.name = name
        self.coords0 = coords0
        self.energy0 = energy0
        self.gradient = gradient
        self.hessian = hessian
        self.natoms = len(coords0)
        self.evals, self.evecs = np.linalg.eigh(self.hessian.reshape([3*self.natoms, 3*self.natoms]))
        self.ievals = np.zeros(len(self.evals), float)
        for i, eigval in enumerate(self.evals):
            if abs(eigval)>1e-10:
                self.ievals[i] = 1.0/eigval
        self.ihessian = np.dot(self.evecs, np.dot(np.diag(self.ievals), self.evecs.T)).reshape([self.natoms, 3, self.natoms, 3])

    def calc_energy(self, coords):
        energy = self.energy0
        dx = (coords - self.coords0).reshape([3*self.natoms])
        energy += np.dot(self.gradient.reshape([3*self.natoms]), dx)
        energy += 0.5*np.dot(dx, np.dot(self.hessian.reshape([3*self.natoms, 3*self.natoms]), dx))
        return energy

    def calc_gradient(self, coords):
        dx = (coords - self.coords0).reshape([3*self.natoms])
        return self.gradient + np.dot(self.hessian.reshape([3*self.natoms, 3*self.natoms]), dx).reshape([self.natoms, 3])

    def calc_hessian(self, coords):
        return self.hessian


class CoulombPart(object):
    def __init__(self, name, coords0, charges, epairs, shift=True):
        self.name = name
        self.charges = charges
        self.epairs = epairs
        self.shift = 0.0
        if shift:
            self.shift = -self.calc_energy(coords0)

    def calc_energy(self, coords):
        energy = self.shift
        for i, qi in enumerate(self.charges):
            for j, qj in enumerate(self.charges):
                if j >= i: break
                if [i, j] in self.epairs or [j, i] in self.epairs: continue
                bond = IC('_internal_ei_bond', [i, j], bond_length)
                energy += qi*qj/bond.value(coords)
        return energy

    def calc_gradient(self, coords):
        grad = np.zeros(3*len(self.charges), float)
        for i, qi in enumerate(self.charges):
            for j, qj in enumerate(self.charges):
                if j >= i: break
                if [i, j] in self.epairs or [j, i] in self.epairs: continue
                bond = IC('_internal_ei_bond', [i, j], bond_length)
                r = bond.value(coords)
                grad += -qi*qj/(r**2)*bond.grad(coords)
        return grad

    def calc_hessian(self, coords):
        hess = np.zeros([3*len(self.charges), 3*len(self.charges)], float)
        for i, qi in enumerate(self.charges):
            for j, qj in enumerate(self.charges):
                if j >= i: break
                if [i, j] in self.epairs or [j, i] in self.epairs: continue
                bond = IC('_internal_ei_bond', [i, j], bond_length)
                r = bond.value(coords)
                qgrad = bond.grad(coords)
                hess += qi*qj/(r**2)*(2.0/r*np.outer(qgrad, qgrad) - bond.hess(coords))
        return hess


class ValencePart(object):
    '''
        A class managing all valence force field terms. This class will mainly
        be used in the second step of the fitting procedure, when the force
        constants are refined at fixed values for the rest values.
    '''
    def __init__(self, system):
        self.vterms = {}
        for icname, ics in sorted(system.ics.iteritems()):
            terms = []
            for ic in ics:
                if icname.startswith('dihed'):
                    #Dihedral potential is determined later based on the geometry
                    terms.append(None)
                else:
                    terms.append(HarmonicTerm(ic, system.ref.coords, None, None))
            self.vterms[icname] = terms

    def determine_dihedral_potentials(self, system, marge2=15*deg, marge3=15*deg):
        '''
            Determine the potential of every dihedral based on the values of
            the dihedral angles in the geometry. First try if a cosine potential
            of the form 0.5*K*[1 - cos(m(psi-psi0))] works well with m=2,3 and
            psi0 = 0,pi/m. If this doesn't work, raise a warning and ignore
            dihedral.
        '''
        maxlength = max([len(icname) for icname in system.ics.keys()]) + 2
        deleted_diheds = []
        for icname in sorted(self.vterms.keys()):
            if not icname.startswith('dihed'):
                continue
            ms = []
            rvs = []
            descr = icname + ' '*(maxlength-len(icname))
            ics = system.ics[icname]
            for ic in ics:
                psi0 = abs(ic.value(system.ref.coords))
                n1 = len(system.nlist[ic.indexes[1]])
                n2 = len(system.nlist[ic.indexes[2]])
                if psi0 >= 0 and psi0 <= max([marge2, marge3]):
                    #use m=3 if at least one of the central atoms
                    #has 4 neighbors
                    if 4 in [n1, n2]:
                        ms.append(3)
                        rvs.append(0.0)
                    #use m=2 if at least one of the central atoms
                    #has 3 neighbors
                    elif 3 in [n1, n2]:
                        ms.append(2)
                        rvs.append(0.0)
                    #use m=1 in all other cases
                    else:
                        ms.append(1)
                        rvs.append(0.0)
                elif psi0 >= 60*deg-marge3 and psi0 <= 60*deg+marge3:
                    ms.append(3)
                    rvs.append(60.0*deg)
                elif psi0 >= 90*deg-marge2 and psi0 <= 90*deg+marge2:
                    ms.append(2)
                    rvs.append(90.0*deg)
                elif psi0 >= 120*deg-marge3 and psi0 <= 120*deg+marge3:
                    ms.append(3)
                    rvs.append(0.0*deg)
                elif psi0 >= 180*deg-marge2 and psi0 <= 180*deg+marge2:
                    #use m=3 if at least one of the central atoms
                    #has 4 neighbors
                    if 4 in [n1, n2]:
                        ms.append(3)
                        rvs.append(60.0*deg)
                    #use m=2 if at least one of the central atoms
                    #has 3 neighbors
                    elif 3 in [n1, n2]:
                        ms.append(2)
                        rvs.append(0.0)
                    #use m=1 in all other cases
                    else:
                        ms.append(1)
                        rvs.append(180.0*deg)
                else:
                    ms.append(-1)
                    rvs.append(np.cos(psi0))
            m = DataArray(ms, unit='au')
            rv = DataArray(rvs, unit='deg')
            if m.mean == -1 or m.std > 0.0 or rv.std > 0.0:
                print '    %s   WARNING: ' % descr +\
                      'could not determine clear trent in dihedral angles, ' +\
                      'dihedral is ignored in force field !!!'
                deleted_diheds.append(icname)
            else:
                print '    %s   0.5*K*[1 - cos(%i(psi - psi0))]    psi0 = %s' % (
                    descr, m.mean, rv.string()
                )
                for i, ic in enumerate(ics):
                    ic.icf = dihed_angle
                    ic.qunit = 'deg'
                    self.vterms[icname][i] = CosineTerm(
                        ic, system.ref.coords, 0.0, rv.mean, m.mean
                    )
        for icname in deleted_diheds:
            del system.ics[icname]
            del self.vterms[icname]

    def _get_nterms(self):
        'Method that returns the number of valence terms in the force field.'
        return len(self.vterms)

    nterms = property(_get_nterms)

    def update_fftable(self, fftab):
        '''
            A method to update all force field parameters (force constants and
            rest values) with the values from the given FFTable.

            **Arguments**

            fftab
                An instance of the FFTable class containing force field
                parameters.
        '''
        for icname in sorted(self.vterms.keys()):
            if not icname in fftab.pars.keys():
                continue
            for term in self.vterms[icname]:
                k, q0 = fftab[icname]
                term.k = k
                term.q0 = q0

    def get_fftable(self):
        '''
            A method to return a FFTable instance containing all force field
            parameters (force constants and rest values).
        '''
        fftab = FFTable()
        for icname in sorted(self.vterms.keys()):
            ks = []
            q0s = []
            ms = []
            for term in self.vterms[icname]:
                ks.append(term.k)
                q0s.append(term.q0)
                if isinstance(term, CosineTerm):
                    ms.append(term.A)
            if isinstance(term, CosineTerm):
                fftab.add(icname,
                    DataArray(data=ks, unit=term.ic.kunit),
                    DataArray(data=q0s, unit=term.ic.qunit),
                    m=DataArray(data=ms, unit='au')
                )
            else:
                fftab.add(icname,
                    DataArray(data=ks, unit=term.ic.kunit),
                    DataArray(data=q0s, unit=term.ic.qunit)
                )
        return fftab

    def update_fcs(self, fcs):
        '''
            A method to update the force constants of the valence terms. The
            ordering of fcs in the input argument should be the same as the
            ordering of sorted(system.ics.keys()).
        '''
        for i, icname in enumerate(sorted(self.vterms.keys())):
            for term in self.vterms[icname]:
                term.k = fcs[i]

    def get_fcs(self):
        '''
            A method to return the force constants of the valence terms. The
            ordering of fcs in the input argument should be the same as the
            ordering of sorted(system.ics.keys()).
        '''
        fcs = np.zeros(self.nterms, float)
        for i, icname in enumerate(sorted(self.vterms.keys())):
            for term in self.vterms[icname]:
                fcs[i] = term.k
        return fcs

    def calc_energy(self, coords):
        energy = 0.0
        for icname, vterms in sorted(self.vterms.iteritems()):
            for vterm in vterms:
                energy += vterm.calc_energy(coords=coords)
        return energy

    def calc_gradient(self, coords):
        natoms = len(coords)
        gradient = np.zeros([natoms, 3], float)
        for icname, vterms in sorted(self.vterms.iteritems()):
            for vterm in vterms:
                gradient += vterm.calc_gradient(coords=coords)
        return gradient

    def calc_hessian(self, coords):
        natoms = len(coords)
        hessian = np.zeros([natoms, 3, natoms, 3], float)
        for icname, vterms in sorted(self.vterms.iteritems()):
            for vterm in vterms:
                hessian += vterm.calc_hessian(coords=coords)
        return hessian



def _calc_ei_grad_hess(system, epairs):
    grad = np.zeros(3*system.natoms, float)
    hess = np.zeros([3*system.natoms, 3*system.natoms], float)
    for i in xrange(system.natoms):
        for j in xrange(i):
            qi = system.charges[i]
            qj = system.charges[j]
            if [i, j] in epairs or [j, i] in epairs: continue
            bond = IC('_inter_ei_bond', [i, j], bond_length)
            r = bond.value(system.ref.coords)
            qgrad = bond.grad(system.ref.coords)
            grad += -qi*qj/(r**2)*qgrad
            hess += qi*qj/(r**2)*(2.0/r*np.outer(qgrad, qgrad) - bond.hess(system.ref.coords))
    return grad, hess


def _get_excluded_pairs(system, exclude_rule):
    epairs = []
    if exclude_rule > 0:
        for bond in system.bonds:
            epairs.append([bond[0], bond[1]])
    if exclude_rule > 1:
        for bend in system.bends:
            epairs.append([bend[0], bend[2]])
    if exclude_rule > 2:
        for dihed in system.diheds:
            epairs.append([dihed[0], dihed[3]])
    return epairs
