from molmod.units import kjmol, angstrom, deg, rad

from common import get_program, get_ref_ffs

def run_example(molecule, ei_rule, ei_scheme, atypes_level):
    #get system and reference force field
    program = get_program(molecule, atypes_level, ei_scheme, ei_rule)
    refs = get_ref_ffs(molecule, atypes_level, ei_scheme, ei_rule)
    print 'Comparing pre ...'
    trajectories = program.generate_trajectories(verbose=False)
    fftab1 = program.estimate_from_pt(trajectories, verbose=False)
    fftab1.print_screen()
    compare(fftab1, refs['pre'])
    print 'Comparing post ...'
    fftab2 = program.refine_cost(verbose=False)
    fftab2.print_screen()
    compare(fftab2, refs['post'])

def compare(fftable, ref, tol_bond_k=1e-1*kjmol/(angstrom**2), tol_bond_q=1e-6*angstrom,
        tol_bend_k=1e-1*kjmol/(rad**2), tol_bend_q=1e-3*deg, tol_dihed_k=1e-1*kjmol,
        tol_oopdist_k=1e-1*kjmol/(angstrom**2), tol_oopdist_q=1e-6*angstrom):
    for icname in sorted(ref.keys()):
        k, q0 = fftable[icname]
        kref, q0ref = ref[icname]
        if icname.startswith('bond'):
            print '    %s k  = %12.6f    ref = %12.6f  [kjmol/A^2]' %(icname, k/(kjmol/angstrom**2), kref/(kjmol/angstrom**2))
            print '    %s r0 = %12.6f    ref = %12.6f  [A]' %(icname, q0/angstrom, q0ref/angstrom)
            print ''
            assert abs(k-kref)<tol_bond_k
            assert abs(q0-q0ref)<tol_bond_q
        elif icname.startswith('angle'):
            print '    %s k      = %12.6f    ref = %12.6f  [kjmol/rad^2]' %(icname, k/(kjmol/rad**2), kref/(kjmol/rad**2))
            print '    %s theta0 = %12.6f    ref = %12.6f  [deg]' %(icname, q0/deg, q0ref/deg)
            print ''
            assert abs(k-kref)<tol_bend_k
            assert abs(q0-q0ref)<tol_bend_q
        elif icname.startswith('dihed'):
            print '    %s k = %12.6f    ref = %12.6f  [kjmol]' %(icname, k/kjmol, kref/kjmol)
            print ''
            assert abs(k-kref)<tol_dihed_k
        elif icname.startswith('oopdist'):
            print '    %s k  = %12.6f    ref = %12.6f  [kjmol/A^2]' %(icname, k/(kjmol/angstrom**2), kref/(kjmol/angstrom**2))
            print '    %s d0 = %12.6f    ref = %12.6f  [A]' %(icname, q0/angstrom, q0ref/angstrom)
            print ''
            assert abs(k-kref)<tol_oopdist_k
            assert abs(q0-q0ref)<tol_oopdist_q
        else:
            raise ValueError('Recieved invalid ic %s' %icname)

def test_water_high_he_noei():
    run_example('water', -1, 'he', 'high')

def test_water_high_he_0():
    run_example('water', 0, 'he', 'high')

def test_water_high_he_1():
    run_example('water', 1, 'he', 'high')

def test_methane_high_he_noei():
    run_example('methane', -1, 'he', 'high')

def test_methane_high_he_0():
    run_example('methane', 0, 'he', 'high')

def test_methane_high_he_1():
    run_example('methane', 1, 'he', 'high')

def test_ethanol_high_he_noei():
    run_example('ethanol', -1, 'he', 'high')

def test_ethanol_high_he_0():
    run_example('ethanol', 0, 'he', 'high')

def test_ethanol_high_he_1():
    run_example('ethanol', 1, 'he', 'high')

def test_ethanol_high_he_2():
    run_example('ethanol', 2, 'he', 'high')

def test_ethanol_high_he_3():
    run_example('ethanol', 3, 'he', 'high')
