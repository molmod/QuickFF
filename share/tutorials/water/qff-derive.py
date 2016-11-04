from yaff import System
from molmod.io.fchk import FCHKFile

from quickff.reference import SecondOrderTaylor, get_ei_ff
from quickff.program import BaseProgram
from quickff.log import log

import h5py as h5

#define class for deriving the force field
class Program(BaseProgram):
    def run(self):
        with log.section('PROGRAM', 2):
            #deriving diagonal force field
            self.do_pt_generate()
            self.do_pt_estimate()
            self.average_pars()
            self.do_hc_estimatefc(['HC_FC_DIAG'])
            #adding and fitting cross terms
            self.do_cross_init()
            self.do_hc_estimatefc(['HC_FC_CROSS'], logger_level=1)
            #write output
            self.make_output()

#load Gaussian Formatted Checkpoint file
fchk = FCHKFile('/home/louis/build/quickff/share/systems/water/gaussian.fchk')
numbers = fchk.fields.get('Atomic numbers')
energy = fchk.fields.get('Total Energy')
coords = fchk.fields.get('Current cartesian coordinates').reshape([len(numbers), 3])
grad = fchk.fields.get('Cartesian Gradient').reshape([len(numbers), 3])
hess = fchk.get_hessian().reshape([len(numbers), 3, len(numbers), 3])

#Construct Yaff System file
system = System(numbers, coords)
system.detect_bonds()
system.set_standard_masses()

#Construct a QuickFF SecondOrderTaylor object containing the AI reference
ai = SecondOrderTaylor('ai', coords=coords, energy=energy, grad=grad, hess=hess)

#define atom types
rules = [
    ('H', '1 & =1%8'), #hydrogen atom with one oxygen neighbor
    ('O', '8 & =2%1'), #oxygen atom with two hydrogen neighbors
]
system.detect_ffatypes(rules)

#construct electrostatic force field from HE charges in gaussian_wpart.h5
with hf.File('/home/louis/build/quickff/share/systems/water/gaussian_mbis.h5') as f:
    charges = f['charges'][:]
scales = [1.0, 1.0, 1.0, 1.0]
ff_ei = get_ei_ff('EI', system, charges, scales, radii=None, average=True, pbc=[0,0,0])

#initialize and run program
program = Program(system, ai, ffrefs=[ff_ei], fn_yaff='pars_cov.txt', plot_traj=True, xyz_traj=True)
program.run()
