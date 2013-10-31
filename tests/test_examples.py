from molmod.units import kjmol, angstrom, deg, rad

from quickff.model import Model
from quickff.program import Program

from common import get_system, get_reference_ff

def run_example(molecule, ei_rule, ei_scheme, atypes_level):
    #get system and reference force field
    system = get_system(molecule, atypes_level, ei_scheme)
    reference = get_reference_ff(molecule, atypes_level, ei_scheme, ei_rule)   
    #Construct model
    if ei_rule==-1:
        ei_pot_kind='Zero'
        ei_scales = [0.0,0.0,0.0]
    elif ei_rule==0:
        ei_pot_kind='Harm'
        ei_scales = [1.0,1.0,1.0]
    elif ei_rule==1:
        ei_pot_kind='Harm'
        ei_scales = [0.0,1.0,1.0]
    elif ei_rule==2:
        ei_pot_kind='Harm'
        ei_scales = [0.0,0.0,1.0]
    elif ei_rule==3:
        ei_pot_kind='Harm'
        ei_scales = [0.0,0.0,0.0]
    else:
        raise ValueError('Invalid ei-rule %i' %ei_rule)
    model = Model.from_system(system, ei_pot_kind=ei_pot_kind, ei_scales=ei_scales)
    model.val.determine_dihedral_potentials(system, verbose=False)
    program = Program(system, model)
    print 'Comparing ff from perturbation theory ...'
    trajectories = program.generate_trajectories(verbose=False)
    fftab1 = program.estimate_from_pt(trajectories, verbose=False)
    compare(fftab1, reference['pt'])
    print 'Comparing ff after refinement ...'
    fftab2 = program.refine_cost(verbose=False)
    compare(fftab2, reference['refined'])

def compare(fftable, ref):
    assert sorted(ref.keys())==sorted(fftable.pars.keys())
    tolerance = {
        'bond'  : [1e-1*kjmol/angstrom**2, 1e-4*angstrom],
        'bend'  : [1e-1*kjmol/rad**2     , 1e-1*deg     ],
        'dihed' : [1e-1*kjmol            , None         ],
        'opdist': [1e-1*kjmol/angstrom**2, 1e-4*angstrom],
    }
    for icname in sorted(ref.keys()):
        k, q0 = fftable[icname]
        kref, q0ref = ref[icname]
        if icname.startswith('bond'):
            print '    %s k  = %12.6f    ref = %12.6f  [kjmol/A^2]' %(icname, k/(kjmol/angstrom**2), kref/(kjmol/angstrom**2))
            print '    %s r0 = %12.6f    ref = %12.6f  [A]' %(icname, q0/angstrom, q0ref/angstrom)
            print ''
            assert abs(k-kref)<tolerance['bond'][0]
            assert abs(q0-q0ref)<tolerance['bond'][1]
        elif icname.startswith('angle'):
            print '    %s k      = %12.6f    ref = %12.6f  [kjmol/rad^2]' %(icname, k/(kjmol/rad**2), kref/(kjmol/rad**2))
            print '    %s theta0 = %12.6f    ref = %12.6f  [deg]' %(icname, q0/deg, q0ref/deg)
            print ''
            assert abs(k-kref)<tolerance['bend'][0]
            assert abs(q0-q0ref)<tolerance['bend'][1]
        elif icname.startswith('dihed'):
            print '    %s k = %12.6f    ref = %12.6f  [kjmol]' %(icname, k/kjmol, kref/kjmol)
            print ''
            assert abs(k-kref)<tolerance['dihed'][0]
        elif icname.startswith('opdist'):
            print '    %s k  = %12.6f    ref = %12.6f  [kjmol/A^2]' %(icname, k/(kjmol/angstrom**2), kref/(kjmol/angstrom**2))
            print '    %s d0 = %12.6f    ref = %12.6f  [A]' %(icname, q0/angstrom, q0ref/angstrom)
            print ''
            assert abs(k-kref)<tolerance['opdist'][0]
            assert abs(q0-q0ref)<tolerance['opdist'][1]
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

def test_ethene_high_he_noei():
    run_example('ethene', -1, 'he', 'high')

def test_ethene_high_he_0():
    run_example('ethene', 0, 'he', 'high')

def test_ethene_high_he_1():
    run_example('ethene', 1, 'he', 'high')

def test_ethene_high_he_2():
    run_example('ethene', 2, 'he', 'high')

def test_ethene_high_he_3():
    run_example('ethene', 3, 'he', 'high')

def test_benzene_high_he_noei():
    run_example('benzene', -1, 'he', 'high')

def test_benzene_high_he_0():
    run_example('benzene', 0, 'he', 'high')

def test_benzene_high_he_1():
    run_example('benzene', 1, 'he', 'high')

def test_benzene_high_he_2():
    run_example('benzene', 2, 'he', 'high')

def test_benzene_high_he_3():
    run_example('benzene', 3, 'he', 'high')

def test_amoniak_high_he_noei():
    run_example('amoniak', -1, 'he', 'high')

def test_amoniak_high_he_0():
    run_example('amoniak', 0, 'he', 'high')

def test_amoniak_high_he_1():
    run_example('amoniak', 1, 'he', 'high')

def test_amoniak_high_he_2():
    run_example('amoniak', 2, 'he', 'high')

def test_amoniak_high_he_3():
    run_example('amoniak', 3, 'he', 'high')
