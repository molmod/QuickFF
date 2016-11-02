# -*- coding: utf-8 -*-
# QuickFF is a code to quickly derive accurate force fields from ab initio input.
# Copyright (C) 2012 - 2016 Louis Vanduyfhuys <Louis.Vanduyfhuys@UGent.be>
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


from quickff.log import log

from molmod.units import kcalmol, angstrom, deg, amu


__all__ = ['dump_charmm22_prm', 'dump_charmm22_psf']



def _check_charm22(valence):
    """Print warnings for all kinds of energy terms not supported by CHARMM22.

       **Arguments**

       valence
            Instance of ValenceFF, which defines the force field.
    """
    # First report all types of terms not supported by CHARMM22.
    skipped = set([])
    supported_kinds = set(['BONDHARM', 'BENDAHARM', 'TORSION'])
    for term in valence.iter_terms():
        kind = term.basename[:term.basename.find('/')].upper()
        if kind not in supported_kinds and kind not in skipped:
            log.warning('Not writing {} term to CHARMM22 parameter file.'.format(kind))
            skipped.add(kind)



# Template based on the PROTEINS-ZINC-NUCLEID-ACIDS force field file by MacKerell.
# (par_all27_prot_na.prm)
charmm22_template = '''\

BONDS
!
!V(bond) = Kb(b - b0)**2
!
!Kb: kcal/mole/A**2
!b0: A
!
!atom type Kb          b0
!
{}


ANGLES
!
!V(angle) = Ktheta(Theta - Theta0)**2
!
!V(Urey-Bradley) = Kub(S - S0)**2
!
!Ktheta: kcal/mole/rad**2
!Theta0: degrees
!Kub: kcal/mole/A**2 (Urey-Bradley)
!S0: A
!
!atom types     Ktheta    Theta0   Kub     S0
!
{}


DIHEDRALS
!
!V(dihedral) = Kchi(1 + cos(n(chi) - delta))
!
!Kchi: kcal/mole
!n: multiplicity
!delta: degrees
!
!atom types             Kchi    n   delta
{}


IMPROPER
!
!V(improper) = Kpsi(psi - psi0)**2
!
!Kpsi: kcal/mole/rad**2
!psi0: degrees
!note that the second column of numbers (0) is ignored
!
!atom types           Kpsi                   psi0
!

'''


def _bonds_to_charmm22(valence):
    """Construct a BONDS section of a CHARMM22 parameter file.

       **Arguments**

       valence
            Instance of ValenceFF, which defines the force field.
    """
    result = []
    for master in valence.iter_masters(label='BONDHARM'):
        ffatypes = master.basename.split('/')[1].split('.')
        K, q0 = valence.get_params(master.index)
        if K < 1e-6*kcalmol/angstrom**2:
            continue
        result.append('{:15s}  {:15s}  {:9.3f}  {:10.4f}'.format(
            ffatypes[0], ffatypes[1], 0.5*K/(kcalmol/angstrom**2), q0/angstrom
        ))
    return '\n'.join(result)


def _angles_to_charmm22(valence):
    """Construct a ANGLES section of a CHARMM22 parameter file.

       **Arguments**

       valence
            Instance of ValenceFF, which defines the force field.
    """
    result = []
    for master in valence.iter_masters(label='BENDAHARM'):
        ffatypes = master.basename.split('/')[1].split('.')
        K, q0 = valence.get_params(master.index)
        if K < 1e-6*kcalmol:
            continue
        result.append('{:15s}  {:15s}  {:15s}  {:9.3f}  {:10.4f}'.format(
            ffatypes[0], ffatypes[1], ffatypes[2], 0.5*K/kcalmol, q0/deg
        ))
    return '\n'.join(result)


def _dihedrals_to_charmm22(valence):
    """Construct a DIHEDRALS section of a CHARMM22 parameter file.

       **Arguments**

       valence
            Instance of ValenceFF, which defines the force field.
    """
    result = []
    for master in valence.iter_masters(label='TORSION'):
        ffatypes = master.basename.split('/')[1].split('.')
        m, K, q0 = valence.get_params(master.index)
        if K < 1e-6*kcalmol: continue
        result.append('{:15s}  {:15s}  {:15s}  {:15s}  {:10.4f}  {:2.0f}  {:8.2f}'.format(
            ffatypes[0], ffatypes[1], ffatypes[2], ffatypes[3], 0.5*K/kcalmol, m, q0/deg+180
        ))
    return '\n'.join(result)


def dump_charmm22_prm(valence, fn):
    """Dump supported parameters in a CHARMM22 parameter file.

       **Arguments**

       valence
            Instance of ValenceFF, which defines the force field.

       fn
            The filename to write to.
    """
    # First report all types of terms not supported by CHARMM22.
    _check_charm22(valence)

    # Dump supported parameters in prm file.
    with open(fn, 'w') as f:
        f.write(charmm22_template.format(
            _bonds_to_charmm22(valence),
            _angles_to_charmm22(valence),
            _dihedrals_to_charmm22(valence)))


psf_template = '''\
PSF

       1 !NTITLE
 REMARKS Generated with QuickFF

{}

{}

{}

{}

       0 !NIMPHI: impropers


       0 !NDON: donors


       0 !NACC: acceptors


       0 !NNB


'''

def _atoms_to_psf(system):
    result = ['{:8d} !NATOM'.format(system.natom)]
    for iatom in xrange(system.natom):
        ffatype = system.get_ffatype(iatom)
        if len(ffatype) > 4:
            log.warning('Atom type too long for CHARMM PSF file: {}'.format(ffatype))
        result.append('{:8d} A    1    MOL  {:4} {:4} {:10.6f} {:13.4f}           0'.format(
            iatom+1, ffatype, ffatype, system.charges[iatom], system.masses[iatom]/amu))
    return '\n'.join(result)


def _ics_to_charmm22(valence, label, maxwidth, header):
    result = []
    width = maxwidth
    nic = 0
    for term in valence.iter_terms(label=label):
        # Extract iatom0, iatom1
        iatoms = term.get_atoms()
        # Predict what the width will be and fix if needed
        width += len(iatoms)
        if width >= maxwidth:
            width = 0
        nic += 1
        # Convert to string and add to result
        s = ''.join(['{:8d}'.format(iatom) for iatom in iatoms])
        if width == 0:
            result.append(s)
        else:
            result[-1] += s
    result.insert(0, header.format(nic))
    return '\n'.join(result)


def _bonds_to_charmm22(valence):
    return _ics_to_charmm22(valence, 'BONDHARM', 8, '{:8d} !NBOND: bonds')


def _angles_to_charmm22(valence):
    return _ics_to_charmm22(valence, 'BENDAHARM', 9, '{:8d} !NTHETA: angles')


def _dihedrals_to_charmm22(valence):
    return _ics_to_charmm22(valence, 'TORSION', 8, '{:8d} !NPHI: dihedrals')


def dump_charmm22_psf(system, valence, fn):
    """Dump supported internal coordinates in a CHARMM psf file.

       **Arguments**

       system
            Instance of yaff.System class

       valence
            Instance of ValenceFF, which defines the force field.

       fn
            The filename to write to.
    """
    # First report all types of terms not supported by CHARMM22.
    _check_charm22(valence)

    # Dump supported internal coordinates into PSF file.
    with open(fn, 'w') as f:
        f.write(psf_template.format(
            _atoms_to_psf(system),
            _bonds_to_charmm22(valence),
            _angles_to_charmm22(valence),
            _dihedrals_to_charmm22(valence)))
