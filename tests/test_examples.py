from molmod.units import kjmol, angstrom, deg, rad
import os

from quickff.model import Model
from quickff.program import Program
from quickff.context import context

from common import get_system, read_qff_out, translate_rule

def run_example(molecule, ei_rule=2, vdw_rule=2, vdw_from='uff', ei_scheme='he', atypes_level='high', verbose=True):
    #Preparing calc
    moldir = context.get_fn('examples/%s' %molecule)
    ei_scales, ei_pot_kind = translate_rule(ei_rule)
    vdw_scales, vdw_pot_kind = translate_rule(vdw_rule)
    command  = 'qff-est.py --atypes-level=%s ' % atypes_level
    command += '--ei-scheme=%s --ei-scales=%s --ei-model=%s ' %(ei_scheme, ','.join([str(x) for x in ei_scales]), ei_pot_kind)
    command += '--vdw-from=%s --vdw-scales=%s --vdw-model=%s ' %(vdw_from, ','.join([str(x) for x in vdw_scales]), vdw_pot_kind)
    command += '%s/gaussian.fchk %s/gaussian.fchk.h5' %(moldir, moldir)
    #Running quickff calc
    cwd = os.getcwd()
    os.system('mkdir -p %s/testwork/%s' %(cwd,molecule))
    os.system('cd %s/testwork/%s; %s > %s_%s%s_%s%s.qff; rm pars_*.txt system.chk; cd %s' %(
        cwd, molecule, command, atypes_level, ei_scheme, ei_rule, vdw_from, vdw_rule, cwd
    ))
    #Reading new quickff calc and comparing with reference
    fn_ref = os.path.join(moldir, '%s_%s%s_%s%s.qff' %(atypes_level, ei_scheme, ei_rule, vdw_from, vdw_rule) )
    fn_new = '%s/testwork/%s/%s_%s%s_%s%s.qff' %(cwd, molecule, atypes_level, ei_scheme, ei_rule, vdw_from, vdw_rule)
    ref = read_qff_out(fn_ref)
    new = read_qff_out(fn_new)
    print 'Comparing ff from perturbation theory ...'
    compare(new['pt'], ref['pt'])
    print 'Comparing ff after refinement ...'
    compare(new['refined'], ref['refined'])
    #Remove directory with temporary files
    os.system('rm -r %s/testwork' %cwd)

def compare(new, ref):
    assert sorted(ref.keys())==sorted(new.keys())
    tolerance = {
        'bond'  : [1e-6*kjmol/angstrom**2, 1e-6*angstrom],
        'bend'  : [1e-6*kjmol/rad**2     , 1e-6*deg     ],
        'dihed' : [1e-6*kjmol            , None         ],
        'opdist': [1e-6*kjmol/angstrom**2, 1e-6*angstrom],
    }
    for icname in sorted(ref.keys()):
        k, q0 = new[icname]
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


#water
def test_water_noei_novdw():
    run_example('water', ei_rule=-1, vdw_rule=-1)

def test_water_ei0_novdw():
    run_example('water', ei_rule=0, vdw_rule=-1)


#methane
def test_methane_noei_novdw():
    run_example('methane', ei_rule=-1, vdw_rule=-1)

def test_methane_ei0_novdw():
    run_example('methane', ei_rule=0, vdw_rule=-1)


#ethanol
def test_ethanol_noei_novdw():
    run_example('ethanol', ei_rule=-1, vdw_rule=-1)

def test_ethanol_ei0_novdw():
    run_example('ethanol', ei_rule=0, vdw_rule=-1)

def test_ethanol_ei2_uff2():
    run_example('ethanol', ei_rule=2, vdw_rule=2)


#ethane
def test_ethene_noei_novdw():
    run_example('ethene', ei_rule=-1, vdw_rule=-1)

def test_ethene_ei0_novdw():
    run_example('ethene', ei_rule=0, vdw_rule=-1)

def test_ethene_ei2_uff2():
    run_example('ethene', ei_rule=2, vdw_rule=2)


#benzene
def test_benzene_noei_novdw():
    run_example('benzene', ei_rule=-1, vdw_rule=-1)

def test_benzene_ei0_novdw():
    run_example('benzene', ei_rule=0, vdw_rule=-1)


def test_benzene_ei2_uff2():
    run_example('benzene', ei_rule=2, vdw_rule=2)


#amoniak
def test_amoniak_noei_novdw():
    run_example('amoniak', ei_rule=-1, vdw_rule=-1)

def test_amoniak_ei0_novdw():
    run_example('amoniak', ei_rule=0, vdw_rule=-1)

def test_amoniak_ei2_uff2():
    run_example('amoniak', ei_rule=2, vdw_rule=2)
