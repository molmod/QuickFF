Tutorial 1 - Benzene
####################

This tutorial will describe in detail the use of QuickFF by means of the command
line scripts, to generate a force field for benzene. This tutorial assumes that
the following input files are available (examples of such file are provided with
the source code in the directory `share/systems/benzene`):

* ``gaussian.fchk``
    A Gaussian Formatted Checkpoint file containing the ab initio equilibrium
    geometry and ab initio Hessian in equilibrium. This file can be generated
    by performing a freq job in Gaussian using your level of theory and basis
    set of choice. The following line should be added prior to the route
    section (i.e. the lines starting with #)::

        %chk=gaussian.chk

    When the calculation terminated succesfully, make the fchk file using the
    following command::

        formchk gaussian.chk gaussian.fchk

* ``gaussian_mbis.h5``
    A HDF5 file containing the atomic charges. Such a file can, for example, be
    generated using `HORTON <http://molmod.github.com/horton/>`_. If you have HORTON 2.x
    installed, the following command will derive the atomic charges in this example from
    the wavefunction in the ``gaussian.fchk`` file::

        horton-wpart.py gaussian.fchk gaussian_mbis.h5 mbis --grid=ultrafine

    To determine
    the path in the HDF5 file to these charges, one can use the `h5dump`
    command (part of the Linux package `hdf5-tools`)::

        h5dump -n gaussian_mbis.h5

    The output of this command for the gaussian_mbis.h5 file provided in the
    directory ``share/systems/benzene`` is:

    .. program-output:: h5dump -n ../share/systems/benzene/gaussian_mbis.h5

    From this output, we can conclude that, for this example, MBIS charges can
    be found in the path ``/charges``.


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

.. program-output:: qff.py --ffatypes=low --suffix=_noei ../share/systems/benzene/gaussian.fchk

The logger will dump to following information to the screen (or a file if the
:option:`--logfile` option was used):

* Machine environment:
    This section contains some general information about the machine environment
    on which you ran the quickff job. It contains the user name, machine info,
    starting time of job, Python version, Numpy version, Scipy versoin,
    Matplotlib version, current directory and the command used
    to envoke quickff.

* Routine sequence:
    Every task that is performed will be shown in the log output. For example,
    the line ``PTGEN  Constructing trajectories`` indicates that QuickFF is
    busy constructing the perturbation trajectories. As a second example, the
    line ``HCEST  Estimating force constants from Hessian cost for tasks
    HC_FC_DIAG HC_FC_CROSS`` indicates QuickFF is estimating force constants by
    means of a least-square fit of the FF Hessian to the AI Hessian for each
    term for which the task `HC_FC_DIAG` and `HC_FC_CROSS` was selected.

* Final force field parameters:
    The final FF parameters are shown at the end of the log output.

* Timings:
    Finally, a summary of the timings for several steps during run of the
    program is provided.

The force field parameters were also written to the Yaff parameters file
`pars_cov.txt`:

.. program-output:: cat pars_cov_noei.txt


Force field with electrostatics
===============================

* Generating Yaff input
    First we need to generate the Yaff input file for the electrostatic
    contribution to the force field. This is done using the script
    :ref:`qff_input-ei.py <seclab_rg_scripts_inputei>`. For this tutorial,
    we will convert the charges given in the dataset ``/charges`` of the
    file ``gaussian_mbis.h5`` for the atoms in gaussian.fchk with atom types
    according to the level `medium` and use Gaussian charge distributions::

        qff-input-ei.py --ffatypes=low --gaussian gaussian.fchk gaussian_mbis.h5:charges

    This command dumped the following output to the screen, indicating wheter or
    not the atom types are well chosen from the point of view of electrostatics
    (see second remark in :ref:`qff-input-ei.py <seclab_ug_tools_inputei>`):

    .. program-output:: qff-input-ei.py --ffatypes=low --gaussian ../share/systems/benzene/gaussian.fchk ../share/systems/benzene/gaussian_mbis.h5:charges

    Furthermore, the following Yaff parameter (`pars_ei.txt`) file was written:

    .. program-output:: cat pars_ei.txt

* Constructing the covalent contribution
    Now, we generate a covalent force field on top of the previously derived
    electrostatic contribution using the qff.py script::

        qff.py --ffatype=low --ei=pars_ei.txt gaussian.fchk

    The logging output for this job is:

    .. program-output:: qff.py --ffatypes=low --suffix=_ei --ei=pars_ei.txt ../share/systems/benzene/gaussian.fchk

    An extra line appeared in the beginning of the log output, i.e.
    ``QFF    Initializing Yaff force field reference for EI``. This indicates
    that an extra reference instance was created to represent the EI
    contribution to the force field. Furthermore, the covalent parameters are
    almost identical compared to the FF without electrostatics. This is indeed
    what we expect due to the charges being so small.

    The force field parameters were also written to the Yaff parameters file
    `pars_cov.txt`:

    .. program-output:: cat pars_cov_ei.txt
