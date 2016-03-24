Tutorial 1 - Benzene
####################

This tutorial will describe in detail the use of QuickFF to generate a force
field for benzene. This tutorial assumes that the following input files are
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

As a first example, we will derive a force field containing only covalent terms.
No electrostatic nor van der Waals interactions will be included in the force 
field. Furthermore, we use the :ref:`built-in feature 
<seclab_ug_atype_estimator>` to automatically determine atom types according to 
the level `low` (for benzene it does not make any difference which level is 
chosen). This can very easily be done using the :ref:`qff.py 
<seclab_rg_scripts_qff>` script (documentation on the available options can be 
accessed by `qff.py` :option:`--help`)::

    qff.py --ffatypes=low gaussian.fchk

The script will dump all relevant information to the screen, for this tutorial,
the output is as follows:

.. program-output:: qff.py --ffatypes=low /home/louis/build/quickff/share/systems/benzene/gaussian.fchk

The logger will dump to following information to the screen (or a file if the
:option:`--logfile` option was used):

* Machine environment:
    This section contains some general information about the machine environment
    on which you ran the quickff job. It contains the user name, machine info,
    starting time of job, Python version, current directory and the command used
    to envoke quickff.


Force field with electrostatics
===============================

TODO
