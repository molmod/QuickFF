#!/usr/bin/env python
'''Remove largest pieces of unused data from gaussian_mbis.h5 files, and repack them.

All gaussian_mbis.h5 files were first generated with the following bash loop:

    for d in */; do (cd $d; horton-wpart.py gaussian.fchk gaussian_mbis.h5 mbis
                     --grid=ultrafine --lmax=1 --overwrite); done

using HORTON 2.0.1. Then this clean_h5.py script is executed to get rid of excess data.
'''

from glob import glob
import os

import h5py as h5

for fn_h5 in glob('*/*.h5'):
    with h5.File(fn_h5, 'r') as f, h5.File(fn_h5 + '.tmp', 'w') as fnew:
        for key in f:
            if isinstance(f[key], h5.Dataset):
                f.copy(key, fnew)
    os.remove(fn_h5)
    os.rename(fn_h5 + '.tmp', fn_h5)
