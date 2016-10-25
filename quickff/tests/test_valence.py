from molmod.units import angstrom, kjmol, rad, deg

from quickff.valence import ValenceFF
from quickff.tools import guess_ffatypes

from itertools import permutations

from common import log, read_system

import numpy as np

def check_terms(name):
    'Check whether all ICs are present in ValenceFF instance'
    #TODO: CROSS terms
    with log.section('NOSETST', 2):
        system, ref = read_system(name)
        guess_ffatypes(system, 'high')
        valence = ValenceFF(system)
    #check if every bond is present and harmonic
    for bond in system.iter_bonds():
        found = False
        for term in valence.iter_terms('BONDHARM'):
            at0, at1 = term.get_atoms()
            if bond[0]==at0 and bond[1]==at1 \
            or bond[0]==at1 and bond[1]==at0:
                assert not found, 'BondHarm term %s was already found!' %str(bond)
                found = True
        assert found, 'No BondHarm term found for bond %s' %str(bond)
    #check if every bend is present
    for angle in system.iter_angles():
        found = False
        for term in valence.iter_terms('BENDAHARM'):
            at0, at1, at2 = term.get_atoms()
            if angle[0]==at0 and angle[1]==at1 and angle[2]==at2 \
            or angle[0]==at2 and angle[1]==at1 and angle[2]==at0:
                assert not found, 'BendAHarm term %s was already found!' %str(angle)
                found = True
        assert found, 'No BendAHarm term found for bond %s' %str(angle)
    #check if every dihedral is present
    for dihed in system.iter_dihedrals():
        found = False
        for term in valence.iter_terms('TORSION'):
            at0, at1, at2, at3 = term.get_atoms()
            if dihed[0]==at0 and dihed[1]==at1 and dihed[2]==at2 and dihed[3]==at3\
            or dihed[0]==at3 and dihed[1]==at2 and dihed[2]==at1 and dihed[3]==at0:
                assert not found, 'Torsion term %s was already found!' %str(dihed)
                found = True
        assert found, 'No Torsion term found for bond %s' %str(dihedral)
    #check if every oop distance is present and Harm for rv of 0 and SQHARM else
    for oop in system.iter_oops():
        found = False
        for term in valence.iter_terms('/OOPDIST'):
            at0, at1, at2, at3 = term.get_atoms()
            for p0, p1, p2 in permutations([at0, at1, at2]):
                if oop[0]==p0 and oop[1]==p1 and oop[2]==p2 and oop[3]==at3:
                    assert not found, 'OopDist term %s was already found!' %str(oop)
                    found = True
        for term in valence.iter_terms('SQOOPDIST'):
            at0, at1, at2, at3 = term.get_atoms()
            for p0, p1, p2 in permutations([at0, at1, at2]):
                if oop[0]==p0 and oop[1]==p1 and oop[2]==p2 and oop[3]==at3:
                    assert not found, 'SqOopDist term %s was already found!' %str(oop)
                    found = True
        assert found, 'No (Sq)OopDist term found for bond %s' %str(oop)

def get_analytic_numeric_hessian(valence, term, **ffpars):
    #setup ff
    valence.set_params(term.index, **ffpars)
    for term2 in valence.iter_terms():
        assert len(term2.slaves)==0
        if term2.index!=term.index:
            valence.set_params(term2.index, rv0=0.0, fc=0.0)
    #print [(i,valence.system.ffatypes[valence.system.ffatype_ids[i]]) for i in term.get_atoms()]
    #compute hcov using built-in function (which uses yaff routine)
    ref = valence.get_hessian_contrib(term.index)
    #compute hcov numerically using valence.calc_energy
    eps = 1e-4
    natoms = len(valence.system.pos)
    num = np.zeros([3*natoms, 3*natoms], float)
    for i in xrange(3*natoms):
        Di = np.zeros(3*natoms, float)
        Di[i] = eps
        Di = Di.reshape([natoms, 3])
        for j in xrange(3*natoms):
            Dj = np.zeros(3*natoms, float)
            Dj[j] = eps
            Dj = Dj.reshape([natoms, 3])
            tmp  = valence.calc_energy(valence.system.pos + Di + Dj)
            tmp -= valence.calc_energy(valence.system.pos + Di - Dj)
            tmp -= valence.calc_energy(valence.system.pos - Di + Dj)
            tmp += valence.calc_energy(valence.system.pos - Di - Dj)
            num[i,j] = tmp/(2.0*eps)**2
    num = 0.5*(num+num.T)
    return ref, num

def check_hessian_bonds(name, tol=1e-3*kjmol/angstrom**2):
    with log.section('NOSETST', 2):
        system, ref = read_system(name)
        guess_ffatypes(system, 'highest')
        valence = ValenceFF(system)
    for term in valence.iter_terms('BONDHARM'):
        rv = np.random.uniform(low=1.00, high=2.00)*angstrom
        fc = np.random.uniform(low=1000, high=3000)*kjmol/angstrom**2
        ref, num = get_analytic_numeric_hessian(valence, term, fc=fc, rv0=rv)
        nonzeros = np.where(ref!=0.0)
        M = (abs(ref-num)).max()
        iM, jM = np.where(abs(ref-num)==M)[0][0], np.where(abs(ref-num)==M)[1][0]
        print 'Bond    %2i (random FC=%8.3f kjmol/A^2    RV=%7.3f A  ):  MaxDev(%2i,%2i)=%.3e kjmol/A^2' %(term.index, fc/(kjmol/angstrom**2), rv/angstrom, iM, jM, M/(kjmol/angstrom**2))
        assert M<tol
    del system, valence, ref, num

def check_hessian_bends(name, tol=1e-3*kjmol/angstrom**2):
    with log.section('NOSETST', 2):
        system, ref = read_system(name)
        guess_ffatypes(system, 'highest')
        valence = ValenceFF(system)
    for term in valence.iter_terms('BENDAHARM'):
        rv = np.random.uniform(low=10, high=170)*deg
        fc = np.random.uniform(low=100, high=1000)*kjmol/rad**2
        ref, num = get_analytic_numeric_hessian(valence, term, fc=fc, rv0=rv)
        nonzeros = np.where(ref!=0.0)
        M = (abs(ref-num)).max()
        iM, jM = np.where(abs(ref-num)==M)[0][0], np.where(abs(ref-num)==M)[1][0]
        print 'Bend    %2i (random FC=%8.3f kjmol/rad^2  RV=%7.3f deg):  MaxDev(%2i,%2i)=%.3e kjmol/A^2' %(term.index, fc/(kjmol/rad**2), rv/deg, iM, jM, M/(kjmol/angstrom**2))
        assert M<tol
    del system, valence, ref, num

def check_hessian_dihedrals(name, tol=1e-3*kjmol/angstrom**2):
    with log.section('NOSETST', 2):
        system, ref = read_system(name)
        guess_ffatypes(system, 'highest')
        valence = ValenceFF(system)
    for term in valence.iter_terms('TORSION'):
        #if not term.index==21: continue
        #vterm = valence.vlist.vtab[term.index]
        #q0 = valence.iclist.ictab[vterm['ic0']]['value']
        rv = np.random.uniform(low=0, high=180)*deg #q0
        fc = np.random.uniform(low=10, high=50)*kjmol
        ref, num = get_analytic_numeric_hessian(valence, term, fc=fc, rv0=rv)
        nonzeros = np.where(ref!=0.0)
        #print ref[nonzeros]/(kjmol/angstrom**2)
        #print num[nonzeros]/(kjmol/angstrom**2)
        M = (abs(ref-num)).max()
        iM, jM = np.where(abs(ref-num)==M)[0][0], np.where(abs(ref-num)==M)[1][0]
        print 'Torsion %2i (random FC=%8.3f kjmol        RV=%7.3f deg):  MaxDev(%2i,%2i)=%.3e kjmol/A^2' %(term.index, fc/kjmol, rv/deg, iM, jM, M/(kjmol/angstrom**2))
        assert M<tol
    del system, valence, ref, num

def check_hessian_oops(name, tol=1e-3*kjmol/angstrom**2):
    with log.section('PROGRAM', 2):
        system, ref = read_system(name)
        guess_ffatypes(system, 'highest')
        valence = ValenceFF(system)
    for term in valence.iter_terms('/OOPDIST'):
        rv = 0.0
        fc = np.random.uniform(low=500, high=5000)*kjmol/angstrom**2
        ref, num = get_analytic_numeric_hessian(valence, term, fc=fc, rv0=rv)
        nonzeros = np.where(ref!=0.0)
        M = (abs(ref-num)).max()
        iM, jM = np.where(abs(ref-num)==M)[0][0], np.where(abs(ref-num)==M)[1][0]
        print 'OopDist %2i (random FC=%8.3f kjmol/A^2    RV=%7.3f A  ):  MaxDev(%2i,%2i)=%.3e kjmol/A^2' %(term.index, fc/(kjmol/angstrom**2), rv/angstrom, iM, jM, M/(kjmol/angstrom**2))
        assert M<tol
    for term in valence.iter_terms('SQOOPDIST'):
        rv = np.random.uniform(low=0.01, high=0.1)*angstrom**2
        fc = np.random.uniform(low=500, high=5000)*kjmol/angstrom**4
        ref, num = get_analytic_numeric_hessian(valence, term, fc=fc, rv0=rv)
        nonzeros = np.where(ref!=0.0)
        M = (abs(ref-num)).max()
        iM, jM = np.where(abs(ref-num)==M)[0][0], np.where(abs(ref-num)==M)[1][0]
        print 'SQOpDst %2i (random FC=%8.3f kjmol/A^4    RV=%7.3f A^2):   MaxDev(%2i,%2i)=%.3e kjmol/A^2' %(term.index, fc/(kjmol/angstrom**4), rv/angstrom**2, iM, jM, M/(kjmol/angstrom**2))
        assert M<tol
    del system, valence, ref, num



def test_terms_water():
    check_terms('water/gaussian.fchk')

def test_terms_methane():
    check_terms('methane/gaussian.fchk')

def test_terms_ethene():
    check_terms('ethene/gaussian.fchk')

def test_terms_ethanol():
    check_terms('ethanol/gaussian.fchk')

def test_terms_amoniak():
    check_terms('amoniak/gaussian.fchk')

def test_terms_benzene():
    check_terms('benzene/gaussian.fchk')




def test_hessian_bonds_water():
    check_hessian_bonds('water/gaussian.fchk')

def test_hessian_bends_water():
    check_hessian_bends('water/gaussian.fchk')


def test_hessian_bonds_methane():
    check_hessian_bonds('methane/gaussian.fchk')

def test_hessian_bends_methane():
    check_hessian_bends('methane/gaussian.fchk')


def test_hessian_bonds_ethane():
    check_hessian_bonds('ethane/gaussian.fchk')

def test_hessian_bends_ethane():
    check_hessian_bends('ethane/gaussian.fchk')

def test_hessian_dihedrals_ethane():
    check_hessian_dihedrals('ethane/gaussian.fchk')


def test_hessian_bonds_ethene():
    check_hessian_bonds('ethene/gaussian.fchk')

def test_hessian_bends_ethene():
    check_hessian_bends('ethene/gaussian.fchk')

def test_hessian_dihedrals_ethene():
    check_hessian_dihedrals('ethene/gaussian.fchk')

def test_hessian_oops_ethene():
    check_hessian_oops('ethene/gaussian.fchk')


def test_hessian_bonds_ethanol():
    check_hessian_bonds('ethanol/gaussian.fchk')

def test_hessian_bends_ethanol():
    check_hessian_bends('ethanol/gaussian.fchk')

def test_hessian_dihedrals_ethanol():
    check_hessian_dihedrals('ethanol/gaussian.fchk')


def test_hessian_bonds_amoniak():
    check_hessian_bonds('amoniak/gaussian.fchk')

def test_hessian_bends_amoniak():
    check_hessian_bends('amoniak/gaussian.fchk')

def test_hessian_oops_amoniak():
    check_hessian_oops('amoniak/gaussian.fchk')


def test_hessian_bonds_benzene():
    check_hessian_bonds('benzene/gaussian.fchk')

def test_hessian_bends_benzene():
    check_hessian_bends('benzene/gaussian.fchk')

def test_hessian_dihedrals_benzene():
    check_hessian_dihedrals('benzene/gaussian.fchk')

def test_hessian_oops_benzene():
    check_hessian_oops('benzene/gaussian.fchk')

