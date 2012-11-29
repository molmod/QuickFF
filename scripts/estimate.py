#! /usr/bin/env python

from lib.system import System, FFTable
import numpy as np, cPickle, os

from init import root_model, sysnames

for name in sysnames:
    with open('%s/data/%s/system.pp' %(root_model,name), 'r') as f: system = cPickle.load(f)
    uf = FFTable(system, coupling=None, free_depth=None)
    ux = FFTable(system, coupling=None, free_depth=3)
    sf = FFTable(system, coupling='symmetric', free_depth=None)
    sx = FFTable(system, coupling='symmetric', free_depth=3)
    with open('%s/data/%s/ff_uf.pp' %(root_model,name), 'w') as f: cPickle.dump(uf, f)
    with open('%s/data/%s/ff_ux.pp' %(root_model,name), 'w') as f: cPickle.dump(ux, f)
    with open('%s/data/%s/ff_sf.pp' %(root_model,name), 'w') as f: cPickle.dump(sf, f)
    with open('%s/data/%s/ff_sx.pp' %(root_model,name), 'w') as f: cPickle.dump(sx, f)
