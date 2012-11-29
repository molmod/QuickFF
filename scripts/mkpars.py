#! /usr/bin/env python

from molmod.units import *
from lib import *

import numpy as np, sys, cPickle
root_model = '/home/louis/mil53/model/metals/covalent/est'

icnames = [
    'bond/%s.O_CA/harm/dist/K'                  ,
    'bond/%s.O_HY/harm/dist/K'                  ,
    'bond/O_CA.C_CA/harm/dist/K'                ,
    'bond/O_HY.H_HY/harm/dist/K'                ,
    'bond/C_PH.C_PC/harm/dist/K'                ,
    'bond/C_PH.C_PH/harm/dist/K'                ,
    'bond/C_PH.H_PH/harm/dist/K'                ,
    'bond/C_PC.C_CA/harm/dist/K'                ,
    'angle/%s.O_HY.%s/harm/angle/K'             ,
    'angle/O_HY.%s.O_CA/harm/angle/K'           ,
    'angle/O_CA.%s.O_CA/sqbend/angle/K'         ,
    'angle/O_CA.C_CA.O_CA/harm/angle/K'         ,
    'angle/C_PC.C_CA.O_CA/harm/angle/K'         ,
    'angle/C_PH.C_PC.C_PH/harm/angle/K'         ,
    'angle/C_PH.C_PH.C_PC/harm/angle/K'         ,
    'angle/C_CA.O_CA.%s/harm/angle/K'           ,
    'angle/C_CA.C_PC.C_PH/harm/angle/K'         ,
    'angle/H_PH.C_PH.C_PH/harm/angle/K'         ,
    'angle/H_PH.C_PH.C_PC/harm/angle/K'         ,
    'angle/H_HY.O_HY.%s/harm/angle/K'           ,
    'dihed/O_CA.%s.O_HY.%s/cos-m2-0/dihed/K'    ,
    'dihed/O_CA.C_CA.O_CA.%s/cos-m2-0/dihed/K'  ,
    'dihed/C_PH.C_PC.C_CA.O_CA/cos-m2-0/dihed/K',
    'dihed/C_PC.C_CA.O_CA.%s/cos-m2-0/dihed/K'  ,
    'dihed/C_PH.C_PH.C_PC.C_PH/cos-m2-0/dihed/K',
    'dihed/C_PC.C_PH.C_PH.C_PC/cos-m2-0/dihed/K',
    'dihed/C_CA.O_CA.%s.O_HY/cos-m2-0/dihed/K'  ,
    'dihed/C_CA.O_CA.%s.O_CA/cos-m2-0/dihed/K'  ,
    'dihed/C_CA.C_PC.C_PH.C_PH/cos-m2-0/dihed/K',
    'dihed/H_PH.C_PH.C_PH.C_PC/cos-m2-0/dihed/K',
    'dihed/H_PH.C_PH.C_PC.C_PH/cos-m2-0/dihed/K',
    'dihed/H_PH.C_PH.C_PC.C_CA/cos-m2-0/dihed/K',
    'dihed/H_PH.C_PH.C_PH.H_PH/cos-m2-0/dihed/K',
    'dihed/H_HY.O_HY.%s.O_CA/cos-m2-0/dihed/K'  ,
]

metal = sys.argv[1]
ffkind = sys.argv[2]
qunits = {'bond': 'A'         , 'angle': 'deg'         , 'dihed': 'deg'  }
kunits = {'bond': 'kjmol/A**2', 'angle': 'kjmol/rad**2', 'dihed': 'kjmol'}
with open('%s/data/%s/linker/ff_%s.pp' %(root_model,metal,ffkind), 'r') as f: linker = cPickle.load(f)
with open('%s/data/%s/oxide/ff_%s.pp' %(root_model,metal,ffkind), 'r') as f: oxide = cPickle.load(f)


f = open('%s/data/%s/pars_%s.txt' %(root_model,metal,ffkind), 'w')
print >> f, """# longname                                     unit                      value # fx/fr
# -----------------------------------------------------------------------------#------"""
for name in icnames:
    icname = name
    if '%s' in name: icname = name %( (metal.upper(),)*icname.count('%') )
    if icname in linker.k_cov and len(linker.k_cov[icname])>0:
        k = linker.k_cov[icname].mean
        q = linker.q_cov[icname].mean
    elif icname in oxide.k_cov and len(oxide.k_cov[icname])>0:
        k = oxide.k_cov[icname].mean
        q = oxide.q_cov[icname].mean
    else:
        print 'Warning: skipping %s, no compatible system found.' %icname
        continue
    
    kunit = kunits[icname.split('/')[0]]
    description  = icname 
    description += ' '*(45-len(icname))
    description += kunit.replace('**', '^')
    description += ' '*(15-len(kunit))
    print >> f, '  %s  % 15.6f' %(description, k/parse_unit(kunit))
    
    qunit = qunits[icname.split('/')[0]]
    description  = icname 
    description += ' '*(45-len(icname))
    description += qunit.replace('**', '^')
    description += ' '*(15-len(qunit))
    print >> f, '  %s  % 15.6f' %(description, q/parse_unit(qunit))
f.close()

