from quickff.valence import ValenceFF
from quickff.tools import guess_ffatypes

from itertools import permutations

from common import log, read_system

def check_terms(name):
    'Check whether all ICs are present in ValenceFF instance'
    #TODO: CROSS terms
    with log.section('PROGRAM', 2):
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
