Tutorial 1 - Benzene
####################

This tutorial will describe in detail the use of QuickFF to generate a force
field for benzene. This tutorial assumes that the following to input files are
available (examples of such file are provided with the source code in the 
directory `share/systems/benzene`):

* gaussian.fchk
    A Gaussian Formatted Checkpoint file containing the ab initio equilibrium
    geometry and ab initio Hessian in equilibrium. This file can be generated
    by performing a freq job in Gaussian using your level of theory and basis
    set of choice. The following line should be added prior to the route
    section (i.e. the lines starting with #)::
    
        %chk=gaussian.chk
    
    When the calculation terminated succesfully, make the fchk file using the
    following command::
    
        formchk gaussian.chk gaussian.fchk

* gaussian_wpart.h5
    A HDF5 file containing the atomic charges. Such a file can, for example, be 
    generated using `Horton <http://molmod.github.com/horton/>`_. To determine 
    the path in the HDF5 file to these charges, one can use the `h5dump` 
    command (part of the Linux package `hdf5-tools`)::
    
        h5dump -n gaussian_wpart.h5
    
    The output of this command for the gaussian_wpart.h5 file provided in the 
    directory `share/systems/benzene` is:

    .. program-output:: h5dump -n /home/louis/build/quickff/share/systems/benzene/gaussian_wpart.h5

    From this output, we can conclude that, for example, Hirshfeld-I charges can
    be found in the path `/wpart/hi/charges`.

Covalent force field without non-bonding terms
==============================================

As a first example, we will derive a force field containing only covalent terms,
i.e. harmonic bonds, harmonic bends, single-cosine dihedrals and harmonic 
out-of-plane distances. No electrostatic nor van der Waals interactions will be 
included in the force field. Furthermore, we use the :ref:`built-in feature 
<seclab_ug_atype_estimator>` to automatically determine atom types according to 
the level `low` (for benzene it does not make any difference which level is 
chosen). This can very easily be done using the :ref:`qff-est.py 
<seclab_rg_qffest>` script (documentation on the available options can be 
accessed by `qff-est.py` :option:`--help`)::

    qff-est.py --ei-model=Zero --vdw-model=Zero --atypes-level=low --ic-ids=all gaussian.fchk

Note that the option :option:`--vdw-model=Zero` is not necessary here, because 
`Zero` is the default value. The same applies for the option 
:option:`--ic-ids=all`. The script will dump all relevant information to the 
screen, for this tutorial, the output is as follows:

.. program-output:: qff-est.py --ei-model=Zero --vdw-model=Zero --atypes-level=low --ic-ids=all /home/louis/build/quickff/share/systems/benzene/gaussian.fchk

This output contains all relevant information concerning: the machine
environment, molecular system, model and resulting force field. Let's break down
the output in these 4 components (as the first 18 lines are a simple header, we
will ignore them in this tutorial):

* Machine environment:
    This section contains some general information about the machine environment
    on which you ran the quickff job. It contains the user name, machine info,
    starting time of job, Python, Numpy and Scipy version, current directory and
    the command used to envoke quickff.
   
* System information:
    This section contains relevant information about the molecular system, this
    includes atom number, atom type, charge, charge radius, van der Waals
    epsilon and van der Waals sigma parameters for every atom in the system.

* Model information:
    This section contains relevant information about ab initio model used as
    training data and the resulting force field model. The following entries 
    can be recognized:
    
    * AI Total project Rot/Trans =  True
        The rotational and translational degrees of freedom are projected out
        of the ab initio Hessian before it is introduced in the Taylor expansion
    
    * AI Total kind = AbInitio (Harmonic)
        The ab initio energy is approximated as a second order Taylor expansion,
        i.e. a harmonic energy expression
        
    * FF Electrostatic scales = 0.00 0.00 1.00
        1-2 and 1-3 electrostatic interactions are rescaled with a factor 0.00,
        i.e. they are switched of, and 1-4 electrostatic interactions are 
        rescaled with a factor 1.00, i.e. they are all fully accounted for.
    
    * FF Electrostatic kind = Zero
        The electrostatic interactions in the force field are modeled by means
        of a Zero model. This means that they are completely neglected
        (disregarding the van der Waals scales values).
    
    * FF van der Waals scales = 0.00 0.00 1.00
        1-2 and 1-3 van der Waals interactions are rescaled by a factor of 0.00,
        hence, they are switched off. 1-4 interactions are fully accounted for.
    
    * FF van der Waals kind = Zero
        The van der Waals interactions in the force field are modeled by means
        of a Zero model. This means that they are completely neglected
        (disregarding the van der Waals scales values).
    
    * FF Covalent kind = TermList
        The covalent part of the force field is modeled by means of a list of
        covalent force field terms. One for each kind of internal coordinates.
        See `FF Covalent term icnames` for a complete list of these terms.
    
    * FF Covalent term icnames:
        A list of all covalent terms included in the force field.

Force field with electrostatics
===============================

As a second example, we will derive a force field containing covalent terms and
electrostatic interactions. We wish to include the electrostatic interactions
between all atom pairs. This means we will use an ei-scale of 1.0 for 1-2, 1-3 
and 1-4 atom pairs (i.e. :option:`--ei-scales=1.0,1.0,1.0`). The charges will be
taken from the Hirshfeld-I path (`/wpart/hi`, see beginning of this section) in 
the gaussian_wpart.h5 file. No van der Waals interactions will be included in 
the force field. Furthermore, we again use the :ref:`built-in feature 
<seclab_ug_atype_estimator>` to automatically determine atom types according to 
the level `low`. This time, the required command is::

    qff-est.py --ei-model=HarmPoint --ei-path=/wpart/hi --ei-scales=1.0,1.0,1.0 --vdw-model=Zero --atypes-level=low --ic-ids=all gaussian.fchk gaussian_wpart.h5

Again, some options are not necessary (:option:`--ei-scales=1.0,1.0,1.0`, 
:option:`--vdw-model=Zero` and :option:`--ic-ids=all`), because their value 
matches the default value. The output of the script this time is as follows:

.. program-output:: qff-est.py --ei-model=HarmPoint --ei-path=/wpart/hi --ei-scales=1.0,1.0,1.0 --vdw-model=Zero --atypes-level=low --ic-ids=all /home/louis/build/quickff/share/systems/benzene/gaussian.fchk /home/louis/build/quickff/share/systems/benzene/gaussian_wpart.h5

In the model information section, we now see an entry::

    FF Electrostatic kind = CoulombPoint (Harmonic)

This means that the electrostatic interactions in the force field are modeled by
means of Coulomb interactions between point charges. However, the Coulombic 
energy expression is Taylor-expanded up to second order. This drastically 
speeds up de processing of the perturbation trajectories without losing
accuracy.
