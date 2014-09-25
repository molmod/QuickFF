Tutorial 1 - Benzene
####################

This tutorial will describe in detail the use of QuickFF to generate a force
field for benzene. This tutorial assumes that the following to input files are
available (examples of such file are provided with the source code in the 
directory `share/systems/benzene`):

* gaussian.fchk
    A Gaussian Formatted Checkpoint file containing the ab initio equilibrium
    geometry and ab initio Hessian in equilibrium. This file can be generated
    by performing a freq job in Gaussian using you level of theory and basis
    set of choice. The following line should be added prior to the route
    section (i.e. the lines starting with #)::
    
        %chk=gaussian.chk
    
    When the calculation terminated succesfully, make the fchk file using the
    following command::
    
        formchk gaussian.chk gaussian.fchk

* gaussian_wpart.h5
    A HDF5 file containing the atomic charges. To determine the path in the HDF5
    file to these charges, one can use the `h5dump` command (part of the Linux
    package `hdf5-tools`)::
    
        h5dump -n gaussian_wpart.h5
    
    The output of this command for the gaussian_wpart.h5 file provided in the 
    directory `share/systems/benzene` is:

    .. program-output:: h5dump -n /home/louis/build/quickff/share/systems/benzene/gaussian_wpart.h5

    From this output, we can conclude that Hirshfeld-I charges can be found in
    the path `/wpart/hi/charges`.

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
