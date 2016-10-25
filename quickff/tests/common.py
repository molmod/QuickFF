from quickff.io import read_abinitio

from quickff.context import context
from quickff.io import read_abinitio
from quickff.reference import SecondOrderTaylor
from quickff.valence import ValenceFF

from yaff import System

from quickff.log import log
log.set_level('silent')

__all__ = ['log', 'read_system']

def read_system(name):
    fn = context.get_fn('systems/%s' %name)
    numbers, coords, energy, grad, hess, masses, rvecs, pbc = read_abinitio(fn)
    system = System(numbers, coords)
    system.detect_bonds()
    system.set_standard_masses()
    ai = SecondOrderTaylor('ai', coords=system.pos.copy(), energy=energy, grad=grad, hess=hess, pbc=pbc)
    return system, ai
