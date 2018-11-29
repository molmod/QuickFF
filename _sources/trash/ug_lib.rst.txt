.. _seclab_ug_lib:

Importing QuickFF as a library
##############################

QuickFF can also be treated as library of classes and methods that is imported 
in a script written by the user. This procedure allows more control over the
core features of QuickFF and also allows more complex force fields. In this
User Guide, we will show how to write a custom script for deriving a force field
using QuickFF. The :ref:`tutorials <seclab_tutorials>` will include examples of
specific systems, such as :ref:`water <seclab_tu_water>` and 
:ref:`biphenyl <seclab_tu_biphenyl>`.

Define the system
=================

First we need to define an instance of the 
`Yaff System class <http://molmod.github.io/yaff/ug_system.html>`_ to define the
molecular system::

    from molmod.units import angstrom
    from yaff import System
    import numpy as np
    #initialize system
    numbers = np.array([1,1])
    coords = np.array([[-0.37*angstrom,0.0,0.0],[0.37*angstrom,0.0,0.0]])
    system = System(numbers, coords, rvecs=None)

Here, we considered a (non-periodic) hydrogen molecule oriented along the 
x-axis and a H-H bond length of 0.74 angstrom.

Define the ab initio reference
==============================

Second, we define the ab initio geometry, gradient and Hessian in equilibrium::

    from quickff.reference import SecondOrderTaylor
    #initialize numpy arrays
    energy = 0.0
    grad = np.zeros([2,3], float)
    hess = np.zeros([...], float)
    #initialize the ab initio reference
    ai = SecondOrderTaylor('ai', coords=coords, energy=energy, grad=grad, hess=hess)

Define force field reference objects
====================================

Next, we define an extra reference object for each force field reference we 
want to include. Suppose we want to take an a priori derived electrostatic
contribution into account. If a Yaff parameter file `pars_ei.txt` is at hand for
the EI contribution, we can add it as follows::

    from yaff.pes.ff import ForceField
    from quickff.reference import YaffForceField
    ff = ForceField.generate(system, 'pars_ei.txt', rcut=10*angstrom)
    ei = YaffForceField('ei', ff)

The keyword argument `rcut` represents the cutoff for the electrostatic
interactions, which was given a high enough value to include all interactions
of our gas-phase system. Instead of reading the electrostatic force field (`ff`
in the code block above) from a given parameter file, one can also define it 
using Yaff as a library. For more information see page in the Yaff manual on 
`Beyond force field parameter files <http://molmod.github.io/yaff/ug_forcefield.html#beyond-force-field-parameter-files>`_.

Define the program and run it
=============================

Finally, one has to define a program, which is a sequence of routines of the 
BaseProgram class (see :ref:`Reference Guide <seclab_rg_modules_program>` for 
more information). For our example of a single hydrogen molecule, one could get
an accurate estimate of both the bond rest value and its force constant from the
perturbation trajectories without force constant refinement::

    from quickff.program import BaseProgram
    class CustomProgram(BaseProgram):
        def run(self):
            self.do_pt_generate()
            self.do_pt_estimate()
            self.make_output()
    
    program = CustomProgram(system, ai, ffrefs=[ei])
    program.run()

The arguments of the initializer of CustomProgram, both mandatory and keyword
arguments, are the same as the BaseProgram class (see
:ref:`Reference Guide  <seclab_rg_modules_program>` for more information).

Logging
=======

Controlling the amount of logging can be done at the beginning of your custom
script as follows::

    from quickff.log import log
    log.set_level(level)

In the second line of code, `level` should be either the string `silent`, `low`,
`medium`, `high` or `highest`, or an integer withing the range 0-4 (0 
corresponds to `silent`, 4 corresponds to `highest`). One can introduce custom
logging sections in the script as follows::

    with log.section(key, print_level, timer=timer_description):
        #insert code here
        ...
        #dump string to logger
        log.dump(string)

The meaning of the various variables is:

* **key**: A short description of the section. This will be repeated in the begining of each line dumped to the logger under that section.

* **print_level**: the minimum level that should be assigned to the logger (through the use of `log.set_level`) to actually print the strings passed through `log.dump`.

* **timer_description**: The description of this section that will be used in the timings of each section. If `None` is given, no timing will be included for this section.

* **string**: a string that will be dumped to the logger

For example, the code below::

    with log.section('TEST', 2, timer='Testing the logger'):
        log.dump('This is a test')

will generate the following line in the log output if `log.level` is set to `medium` or higher::

     TEST   This is a test

and the following line at the end of the log output (actual timing below is not
representative)::

     TIMING Testing the logger              0:00:00.000001

