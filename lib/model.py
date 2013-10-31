from molmod.units import deg
from molmod.ic import bond_length, dihed_angle

import numpy as np

from quickff.ic import IC
from quickff.terms import HarmonicTerm, CosineTerm
from quickff.fftable import DataArray, FFTable

__all__ = [
    'Model', 'AIPart', 'EIPart', 'ValencePart',
    'ZeroPot', 'HarmonicPot', 'CoulombPot', 'TermListPot',
]


class Model(object):
    def __init__(self, ai, val, ei, vdw):
        '''
           A class defining the ab initio total energy of the system,
           the force field electrostatic contribution and the
           force field valence contribution.

           **Arguments**

           ai
                A model for the ab initio total energy, should be an instance
                of HarmonicPart

           val
                an instance of the ValencePart class containing all details
                of the force field valence terms.

           ei
                A model for the force field electrostatic energy, should be
                an instance of HarmonicPart or CoulombPart

           vdw
                A model for the force field van der Waals energy, should be
                an instance of HarmonicPart or LennartJonesPart
        '''
        self.ai = ai
        self.val = val
        self.ei = ei
        self.vdw = vdw

    @classmethod
    def from_system(cls, system, ei_pot_kind='Harmonic', ei_scales=[1.0,1.0,1.0],
        vdw_scales=[0.0,0.0,1.0], vdw_pot_kind='Harmonic', ai_project=True):
        '''
            **Arguments**

            system
                An instance of the System class containing all the
                information of the system.

            **Optional Arguments**

            ei_pot_kind
                a string defining the potential kind of the electrostatic
                interactions. Can be 'Zero', 'Harmonic' or 'Coulomb'. If Coulomb
                is chosen, the exact Coulombic potential will be used to
                evaluate EI interactions. If Harmonic is chosen, a second
                order Taylor expansion is used. Harmonic is a lot faster and
                should already give accurate results.

            ei_scales
                a list containing the scales for the 1-2, 1-3 and 1-4
                contribution to the electrostatic interactions

            vdw_model
                a string defining the model kind of the van der Waals
                interactions. Can be Harmonic, Exact or Zero. If Exact
                is chose, the exact  potential will be used to evaluate
                van der Waals interactions. If Harmonic is chosen, a second
                order Taylor expansion is used. Harmonic is a lot faster and
                should already give accurate results.

            vdw_scales
                a list containing the scales for the 1-2, 1-3 and 1-4
                contribution to the van der Waals interactions

            ai_project
                If True, project the translational and rotational
                degrees of freedom out of the hessian.
        '''
        ai  = AIPart.from_system(system, ai_project)
        ei  = EIPart.from_system(system, ei_scales, ei_pot_kind)
        vdw = VDWPart.from_system(system, vdw_scales, vdw_pot_kind)
        val = ValencePart.from_system(system)
        return cls(ai, val, ei, vdw)

    def print_info(self):
        self.ai.print_info()
        self.ei.print_info()
        self.vdw.print_info()
        self.val.print_info()


#####  PES Parts


class BasePart(object):
    '''
        A base class for several parts to the AI and/or FF PES
    '''
    def __init__(self, name, pot):
        self.name = name
        self.pot = pot

    def calc_energy(self, coords):
        return self.pot.calc_energy(coords)

    def calc_gradient(self, coords):
        return self.pot.calc_gradient(coords)

    def calc_hessian(self, coords):
        return self.pot.calc_hessian(coords)

    def print_info(self):
        print '    %s kind = %s' %(self.name, self.pot.kind)


class AIPart(BasePart):
    '''
        A class for describing the Ab initio PES.
    '''
    def __init__(self, pot, project=None):
        BasePart.__init__(self, 'AI Total', pot)
        self.project = project

    @classmethod
    def from_system(cls, system, project=True):
        if project:
            hess = system.ref.phess.copy()
        else:
            hess = system.ref.hess.copy()
        pot = HarmonicPot(
            system.ref.coords.copy(), 0.0,
            system.ref.grad.copy(), hess
        )
        return cls(pot, project)

    def print_info(self):
        print '    %s project Rot/Trans = ' %self.name, self.project
        BasePart.print_info(self)


class EIPart(BasePart):
    '''
        A class for describing the electrostatic part to the FF PES.
    '''
    def __init__(self, pot, scales):
         BasePart.__init__(self, 'FF Electrostatic', pot)
         self.scales = scales

    @classmethod
    def from_system(cls, system, scales, pot_kind):
        if pot_kind.lower() == 'zero':
            pot = ZeroPot()
        else:
            #Generate list of atom pairs subject ot EI-scaling according to scales
            scaled_pairs = [[],[],[]]
            for bond in system.bonds:
                scaled_pairs[0].append([bond[0], bond[1]])
            for bend in system.bends:
                scaled_pairs[1].append([bend[0], bend[2]])
            for dihed in system.diheds:
                scaled_pairs[2].append([dihed[0], dihed[3]])
            #Construct the potential
            exact = CoulombPot(system.charges.copy(), scales, scaled_pairs, coords0=system.ref.coords.copy())
            if pot_kind.lower() in ['harm', 'harmonic']:
                grad = exact.calc_gradient(system.ref.coords.copy())
                hess = exact.calc_hessian(system.ref.coords.copy())
                pot = HarmonicPot(system.ref.coords.copy(), 0.0, grad, hess)
            elif pot_kind.lower() in ['coul', 'coulomb', 'exact']:
                pot = exact
        return cls(pot, scales)

    def print_info(self):
        print '    %s scales = %.2f %.2f %.2f' %(self.name, self.scales[0], self.scales[1], self.scales[2])
        BasePart.print_info(self)


class VDWPart(BasePart):
    '''
        A class for describing the van der Waals part to the FF PES.
    '''
    def __init__(self, pot, scales):
         BasePart.__init__(self, 'FF van der Waals', pot)
         self.scales = scales

    @classmethod
    def from_system(cls, system, scales, pot_kind):
        if pot_kind.lower() == 'zero':
            pot = ZeroPot()
        elif pot_kind.lower() in ['exact', 'harmonic', 'harm']:
            #Generate list of atom pairs subject ot VDW-scaling according to scales
            scaled_pairs = [[],[],[]]
            for bond in system.bonds:
                scaled_pairs[0].append([bond[0], bond[1]])
            for bend in system.bends:
                scaled_pairs[1].append([bend[0], bend[2]])
            for dihed in system.diheds:
                scaled_pairs[2].append([dihed[0], dihed[3]])
            #Construct the potential
            exact = LennartJonesPot(system.sigmas.copy(), system.epsilons.copy(), scales, scaled_pairs, coords0=system.ref.coords.copy())
            if pot_kind.lower() in ['harm', 'harmonic']:
                grad = exact.calc_gradient(system.ref.coords.copy())
                hess = exact.calc_hessian(system.ref.coords.copy())
                pot = HarmonicPot(system.ref.coords.copy(), 0.0, grad, hess)
            elif pot_kind.lower() in ['lj', 'lennartjones', 'exact']:
                pot = exact
        return cls(pot, scales)

    def print_info(self):
        print '    %s scales = %.2f %.2f %.2f' %(self.name, self.scales[0], self.scales[1], self.scales[2])
        BasePart.print_info(self)


class ValencePart(BasePart):
    '''
        A class managing all valence force field terms. This class will mainly
        be used in the second step of the fitting procedure, when the force
        constants are refined at fixed values for the rest values.
    '''
    def __init__(self, pot):
        BasePart.__init__(self, 'FF Covalent', pot)

    @classmethod
    def from_system(cls, system):
        vterms = {}
        for icname, ics in sorted(system.ics.iteritems()):
            terms = []
            for ic in ics:
                if icname.startswith('dihed'):
                    #Dihedral potential is determined later based on the geometry
                    terms.append(None)
                else:
                    terms.append(HarmonicTerm(ic, system.ref.coords, None, None))
            vterms[icname] = terms
        pot = TermListPot(vterms)
        return cls(pot)

    def determine_dihedral_potentials(self, system, marge2=15*deg, marge3=15*deg, verbose=True):
        '''
            Determine the potential of every dihedral based on the values of
            the dihedral angles in the geometry. First try if a cosine potential
            of the form 0.5*K*[1 - cos(m(psi-psi0))] works well with m=2,3 and
            psi0 = 0,pi/m. If this doesn't work, raise a warning and ignore
            dihedral.
        '''
        maxlength = max([len(icname) for icname in system.ics.keys()]) + 2
        deleted_diheds = []
        for icname in sorted(self.pot.terms.keys()):
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
                if verbose:
                    print '    %s   WARNING: ' % descr +\
                          'could not determine clear trent in dihedral angles, ' +\
                          'dihedral is ignored in force field !!!'
                deleted_diheds.append(icname)
            else:
                if verbose:
                    print '    %s   0.5*K*[1 - cos(%i(psi - %5.1f))]' % (
                        descr, m.mean, rv.mean
                    )
                for i, ic in enumerate(ics):
                    ic.icf = dihed_angle
                    ic.qunit = 'deg'
                    self.pot.terms[icname][i] = CosineTerm(
                        ic, system.ref.coords, 0.0, rv.mean, m.mean
                    )
        for icname in deleted_diheds:
            del system.ics[icname]
            del self.pot.terms[icname]

    def _get_nterms(self):
        'Method that returns the number of valence terms in the force field.'
        return len(self.pot.terms)

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
        for icname in sorted(self.pot.terms.keys()):
            if not icname in fftab.pars.keys():
                continue
            for term in self.pot.terms[icname]:
                k, q0 = fftab[icname]
                term.k = k
                term.q0 = q0

    def get_fftable(self):
        '''
            A method to return a FFTable instance containing all force field
            parameters (force constants and rest values).
        '''
        fftab = FFTable()
        for icname in sorted(self.pot.terms.keys()):
            ks = []
            q0s = []
            ms = []
            for term in self.pot.terms[icname]:
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
        for i, icname in enumerate(sorted(self.pot.terms.keys())):
            for term in self.pot.terms[icname]:
                term.k = fcs[i]

    def get_fcs(self):
        '''
            A method to return the force constants of the valence terms. The
            ordering of fcs in the input argument should be the same as the
            ordering of sorted(system.ics.keys()).
        '''
        fcs = np.zeros(self.nterms, float)
        for i, icname in enumerate(sorted(self.pot.terms.keys())):
            for term in self.pot.terms[icname]:
                fcs[i] = term.k
        return fcs


#####  Potentials


class BasePot(object):
    def __init__(self, kind):
        self.kind = kind

    def  calc_energy(self, coords):
        raise NotImplementedError

    def  calc_gradient(self, coords):
        raise NotImplementedError

    def  calc_hessian(self, coords):
        raise NotImplementedError


class ZeroPot(BasePot):
    def __init__(self):
        BasePot.__init__(self, 'Zero')

    def calc_energy(self, coords):
        return 0.0

    def calc_gradient(self, coords):
        return np.zeros([len(coords), 3], float)

    def calc_hessian(self, coords):
        return np.zeros([len(coords), 3, len(coords), 3], float)


class HarmonicPot(BasePot):
    def __init__(self, coords0, energy0, grad0, hess0):
        BasePot.__init__(self, 'Harmonic')
        self.coords0 = coords0
        self.energy0 = energy0
        self.grad0 = grad0
        self.hess0 = hess0
        self.natoms = len(coords0)

    def calc_energy(self, coords):
        energy = self.energy0
        dx = (coords - self.coords0).reshape([3*self.natoms])
        energy += np.dot(self.grad0.reshape([3*self.natoms]), dx)
        energy += 0.5*np.dot(dx, np.dot(self.hess0.reshape([3*self.natoms, 3*self.natoms]), dx))
        return energy

    def calc_gradient(self, coords):
        dx = (coords - self.coords0).reshape([3*self.natoms])
        return self.grad0 + np.dot(self.hess0.reshape([3*self.natoms, 3*self.natoms]), dx)

    def calc_hessian(self, coords):
        return self.hess0


class CoulombPot(BasePot):
    def __init__(self, charges, scales, scaled_pairs, coords0=None):
        BasePot.__init__(self, 'Coulomb')
        self.charges = charges
        self.scales = scales
        self.scaled_pairs = scaled_pairs
        self.shift = 0.0
        if coords0 is not None:
            self.shift = -self.calc_energy(coords0)

    def _get_scale(self, i, j):
        if [i, j] in self.scaled_pairs[0] or [j, i] in self.scaled_pairs[0]:
            return self.scales[0]
        elif [i, j] in self.scaled_pairs[1] or [j, i] in self.scaled_pairs[1]:
            return self.scales[1]
        elif [i, j] in self.scaled_pairs[2] or [j, i] in self.scaled_pairs[2]:
            return self.scales[2]
        else:
            return 1.0

    def calc_energy(self, coords):
        energy = self.shift
        for i, qi in enumerate(self.charges):
            for j, qj in enumerate(self.charges):
                if j >= i: break
                scale = self._get_scale(i, j)
                if scale==0.0: continue
                bond = IC('_internal_ei_bond', [i, j], bond_length)
                energy += qi*qj/bond.value(coords)*scale
        return energy

    def calc_gradient(self, coords):
        grad = np.zeros(3*len(self.charges), float)
        for i, qi in enumerate(self.charges):
            for j, qj in enumerate(self.charges):
                if j >= i: break
                scale = self._get_scale(i, j)
                if scale==0.0: continue
                bond = IC('_internal_ei_bond', [i, j], bond_length)
                r = bond.value(coords)
                grad += -qi*qj/(r**2)*bond.grad(coords)*scale
        return grad

    def calc_hessian(self, coords):
        hess = np.zeros([3*len(self.charges), 3*len(self.charges)], float)
        for i, qi in enumerate(self.charges):
            for j, qj in enumerate(self.charges):
                if j >= i: break
                scale = self._get_scale(i, j)
                if scale==0.0: continue
                bond = IC('_internal_ei_bond', [i, j], bond_length)
                r = bond.value(coords)
                qgrad = bond.grad(coords)
                hess += qi*qj/(r**2)*(2.0/r*np.outer(qgrad, qgrad) - bond.hess(coords))*scale
        return hess



class LennartJonesPot(BasePot):
    def __init__(self, sigmas, epsilons, scales, scaled_pairs, coords0=None):
        BasePot.__init__(self, 'LennartJones')
        self.sigmas = sigmas
        self.epsilons = epsilons
        self.scales = scales
        self.scaled_pairs = scaled_pairs
        self.shift = 0.0
        if coords0 is not None:
            self.shift = -self.calc_energy(coords0)

    def _get_scale(self, i, j):
        if [i, j] in self.scaled_pairs[0] or [j, i] in self.scaled_pairs[0]:
            return self.scales[0]
        elif [i, j] in self.scaled_pairs[1] or [j, i] in self.scaled_pairs[1]:
            return self.scales[1]
        elif [i, j] in self.scaled_pairs[2] or [j, i] in self.scaled_pairs[2]:
            return self.scales[2]
        else:
            return 1.0

    def calc_energy(self, coords):
        energy = self.shift
        for i, (si, ei) in enumerate(zip(self.sigmas, self.epsilons)):
            for j, (sj, ej) in enumerate(zip(self.sigmas, self.epsilons)):
                if j >= i: break
                scale = self._get_scale(i, j)
                if scale==0.0: continue
                sigma = 0.5*(si+sj)
                epsilon = np.sqrt(ei*ej)
                bond = IC('_internal_ei_bond', [i, j], bond_length)
                x = (sigma/bond.value(coords))
                energy += 4.0*epsilon*(x**12-x**6)*scale
        return energy

    def calc_gradient(self, coords):
        grad = np.zeros(3*len(coords), float)
        for i, (si, ei) in enumerate(zip(self.sigmas, self.epsilons)):
            for j, (sj, ej) in enumerate(zip(self.sigmas, self.epsilons)):
                if j >= i: break
                scale = self._get_scale(i, j)
                if scale==0.0: continue
                sigma = 0.5*(si+sj)
                epsilon = np.sqrt(ei*ej)
                bond = IC('_internal_ei_bond', [i, j], bond_length)
                x = (sigma/bond.value(coords))
                grad -= 24.0*epsilon/sigma*(2.0*x**13-x**7)*bond.grad(coords)*scale
        return grad

    def calc_hessian(self, coords):
        hess = np.zeros([3*len(coords), 3*len(coords)], float)
        for i, (si, ei) in enumerate(zip(self.sigmas, self.epsilons)):
            for j, (sj, ej) in enumerate(zip(self.sigmas, self.epsilons)):
                if j >= i: break
                scale = self._get_scale(i, j)
                if scale==0.0: continue
                sigma = 0.5*(si+sj)
                epsilon = np.sqrt(ei*ej)
                bond = IC('_internal_ei_bond', [i, j], bond_length)
                x = (sigma/bond.value(coords))
                qgrad = bond.grad(coords)
                hess += 24.0*epsilon/sigma**2*(26*x**14-7*x**8)*np.outer(qgrad, qgrad)*scale
                hess -= 24.0*epsilon/sigma*(2*x**13-x**7)*bond.hess(coords)*scale
        return hess


class TermListPot(BasePot):
    '''
        A potential specified by a list of terms to describe the FF valence
        energy.
    '''
    def __init__(self, terms):
        BasePot.__init__(self, 'TermList')
        self.terms = terms

    def calc_energy(self, coords):
        energy = 0.0
        for icname, terms in sorted(self.terms.iteritems()):
            for term in terms:
                energy += term.calc_energy(coords=coords)
        return energy

    def calc_gradient(self, coords):
        natoms = len(coords)
        gradient = np.zeros([natoms, 3], float)
        for icname, terms in sorted(self.terms.iteritems()):
            for term in terms:
                gradient += term.calc_gradient(coords=coords)
        return gradient

    def calc_hessian(self, coords):
        natoms = len(coords)
        hessian = np.zeros([natoms, 3, natoms, 3], float)
        for icname, terms in sorted(self.terms.iteritems()):
            for term in terms:
                hessian += term.calc_hessian(coords=coords)
        return hessian
