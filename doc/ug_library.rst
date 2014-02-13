Importing QuickFF as a library
##############################

QuickFF can also be used as a library that can be imported in external scripts.
This approach is particulary usefull to process large amounts of molecules,
derive a force field and perform force field simulations for each of them. The 
basics if such an approach will be discussed in this section. 

Defining the system
===================

One should always start by defining the system of interest, which is done by
creating an instance of the System class. This can be done by calling the
:meth:`quickff.system.System` constructor::

    from quickff import *
    import numpy as np

    #define potential energy surface
    coords = np.array([])
    grad = np.array([])
    hess = np.array([])
    ref = ReferenceData(coords=coords, grad=grad, hess=hess)

    #define molecular system
    numbers = np.array([])
    bonds = np.array([])
    system = System(numbers, ref)

    #guess atom types automatically
    system.guess_ffatypes('high')
    
    #determine internal coordinates from topology
    system.determine_ics_from_topology()


The system can also be defined by reading input files (similarly as with the
`qff-est.py` script)::

    from quickff import *
    
    #read system from input files
    system = System.from_files(['gaussian.fchk', 'gaussian_wpart.h5'], ei_scheme='he')
    
    #guess atom types
    system.guess_ffatypes('high')

    #determine internal coordinates from topology
    system.determine_ics_from_topology()

Defining the model and program
==============================

Next, a model can be constucted using the method 
:meth:`quickff.model.Model.from_system`::

    model = Model.from_system(system)

Define an instance of the class :class:`quickff.program.Program` to manage the
force field derivation::

    program = Program(system, model)

Constructing the force field
============================
    
Now we can apply the QuickFF methodology and derive a force field for the 
system. This can be done all at once using the method 
:meth:`quickff.program.Program.run`::

    fftab = program.run()

Or the user can perform every step manually::

    #estimate rest angle and multiplicity of dihedral potentials from geometry
    model.val.determine_dihedral_potentials(system)
    
    #Determine the coordinates of the perturbation trajectories
    trajectories = program.generate_trajectories()
    
    #Estimating all pars for bonds, bends and opdists
    fftab = program.estimate_from_pt(trajectories)
    
    #Refining force constants using a Hessian LSQ cost
    fftab = program.refine_cost()

Generating output
=================

The force field parameters can be dumped to the screen::

    fftab.print_screen()

The parameters can also be writen to files in various formats::

    fftab.dump_ffit2('pars_ffit2.txt')
    fftab.dump_yaff('pars_yaff.txt')
