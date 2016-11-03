from contextlib import contextmanager
from glob import glob
import os
import shutil
import tempfile

import h5py as h5

from quickff.io import read_abinitio

from quickff.context import context
from quickff.io import read_abinitio
from quickff.reference import SecondOrderTaylor
from quickff.valence import ValenceFF

from yaff import System

from quickff.log import log
log.set_level('silent')

__all__ = ['log', 'read_system', 'tmpdir']

def read_system(name):
    # Load system data.
    fn = context.get_fn(os.path.join('systems', name))
    numbers, coords, energy, grad, hess, masses, rvecs, pbc = read_abinitio(fn)
    # Try to load charges.
    charges = None
    fns_wpart = glob(os.path.join(os.path.dirname(fn), '*wpart.h5'))
    print os.path.join(os.path.dirname(name), '*wpart.h5')
    if len(fns_wpart) > 0:
        with h5.File(fns_wpart[0], 'r') as f:
            if 'wpart/hi/charges' in f:
                charges = f['wpart/hi/charges'][:]
    # Create system object.
    system = System(numbers, coords, charges=charges)
    system.detect_bonds()
    system.set_standard_masses()
    # Load ab initio data.
    ai = SecondOrderTaylor('ai', coords=system.pos.copy(), energy=energy, grad=grad, hess=hess, pbc=pbc)
    return system, ai

@contextmanager
def tmpdir(name):
    """Create a temporary directory where output files of a test can be written to.

    The name argument can be used to obtained a recognizable test directory name. The
    temporary directory and its contents are automatically removed.
    """
    dn = tempfile.mkdtemp(name)
    try:
        yield dn
    finally:
        shutil.rmtree(dn)
