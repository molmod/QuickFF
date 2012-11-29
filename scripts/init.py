#! /usr/bin/env python

from lib.system import System
import numpy as np, cPickle, os

verbose = False
do_ic_tests = False
root_model = '/home/louis/mil53/model/metals/covalent/est'
root_data = '/home/louis/mil53/data/metals'
sysnames = ['Al/linker', 'Al/oxide', 'Cr/linker', 'Cr/oxide', 'Ga/linker', 'Ga/oxide', 'In/linker', 'In/oxide']
icnames = {
    'bond/%s.O_CA/harm/dist/K'                  : ('kjmol/A**2'  , 'A'   ),
    'bond/%s.O_HY/harm/dist/K'                  : ('kjmol/A**2'  , 'A'   ),
    'bond/O_CA.C_CA/harm/dist/K'                : ('kjmol/A**2'  , 'A'   ),
    'bond/O_HY.H_HY/harm/dist/K'                : ('kjmol/A**2'  , 'A'   ),
    'bond/C_PH.C_PC/harm/dist/K'                : ('kjmol/A**2'  , 'A'   ),
    'bond/C_PH.C_PH/harm/dist/K'                : ('kjmol/A**2'  , 'A'   ),
    'bond/C_PH.H_PH/harm/dist/K'                : ('kjmol/A**2'  , 'A'   ),
    'bond/C_PC.C_CA/harm/dist/K'                : ('kjmol/A**2'  , 'A'   ),
    'angle/%s.O_HY.%s/harm/angle/K'             : ('kjmol/rad**2', 'deg' ),
    'angle/O_HY.%s.O_CA/harm/angle/K'           : ('kjmol/rad**2', 'deg' ),
    'angle/O_CA.%s.O_CA/sqbend/angle/K'         : ('kjmol/rad**2', 'deg' ),
    'angle/O_CA.C_CA.O_CA/harm/angle/K'         : ('kjmol/rad**2', 'deg' ),
    'angle/C_PC.C_CA.O_CA/harm/angle/K'         : ('kjmol/rad**2', 'deg' ),
    'angle/C_PH.C_PC.C_PH/harm/angle/K'         : ('kjmol/rad**2', 'deg' ),
    'angle/C_PH.C_PH.C_PC/harm/angle/K'         : ('kjmol/rad**2', 'deg' ),
    'angle/C_CA.O_CA.%s/harm/angle/K'           : ('kjmol/rad**2', 'deg' ),
    'angle/C_CA.C_PC.C_PH/harm/angle/K'         : ('kjmol/rad**2', 'deg' ),
    'angle/H_PH.C_PH.C_PH/harm/angle/K'         : ('kjmol/rad**2', 'deg' ),
    'angle/H_PH.C_PH.C_PC/harm/angle/K'         : ('kjmol/rad**2', 'deg' ),
    'angle/H_HY.O_HY.%s/harm/angle/K'           : ('kjmol/rad**2', 'deg' ),
    'dihed/O_CA.%s.O_HY.%s/cos-m2-0/dihed/K'    : ('kjmol/rad**2', 'deg' ),
    'dihed/O_CA.C_CA.O_CA.%s/cos-m2-0/dihed/K'  : ('kjmol/rad**2', 'deg' ),
    'dihed/C_PH.C_PC.C_CA.O_CA/cos-m2-0/dihed/K': ('kjmol/rad**2', 'deg' ),
    'dihed/C_PC.C_CA.O_CA.%s/cos-m2-0/dihed/K'  : ('kjmol/rad**2', 'deg' ),
    'dihed/C_PH.C_PH.C_PC.C_PH/cos-m2-0/dihed/K': ('kjmol/rad**2', 'deg' ),
    'dihed/C_PC.C_PH.C_PH.C_PC/cos-m2-0/dihed/K': ('kjmol/rad**2', 'deg' ),
    'dihed/C_CA.O_CA.%s.O_HY/cos-m2-0/dihed/K'  : ('kjmol/rad**2', 'deg' ),
    'dihed/C_CA.O_CA.%s.O_CA/cos-m2-0/dihed/K'  : ('kjmol/rad**2', 'deg' ),
    'dihed/C_CA.C_PC.C_PH.C_PH/cos-m2-0/dihed/K': ('kjmol/rad**2', 'deg' ),
    'dihed/H_PH.C_PH.C_PH.C_PC/cos-m2-0/dihed/K': ('kjmol/rad**2', 'deg' ),
    'dihed/H_PH.C_PH.C_PC.C_PH/cos-m2-0/dihed/K': ('kjmol/rad**2', 'deg' ),
    'dihed/H_PH.C_PH.C_PC.C_CA/cos-m2-0/dihed/K': ('kjmol/rad**2', 'deg' ),
    'dihed/H_PH.C_PH.C_PH.H_PH/cos-m2-0/dihed/K': ('kjmol/rad**2', 'deg' ),
    'dihed/H_HY.O_HY.%s.O_CA/cos-m2-0/dihed/K'  : ('kjmol/rad**2', 'deg' ),
}

if __name__=='__main__':
    for name in sysnames:
        os.system('rm -r %s/%s' %(root_model,name))
        os.system('mkdir -p %s/data/%s' %(root_model,name))
        system = System(name, '%s/%s/freq.mfs' %(root_data,name), fn_psf='%s/%s/opt.psf' %(root_data,name), ei_exclude=0)
        if verbose:
            system.totharm.print_hessian()
            system.eiharm.print_hessian()
        system.find_ic_patterns(icnames)
        if do_ic_tests:
            system.test_ics()
        with open('%s/data/%s/system.pp' %(root_model,name), 'w') as f: cPickle.dump(system, f)
