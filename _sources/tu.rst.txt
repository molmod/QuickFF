Tutorials
#########

Tutorial 1 - Benzene
********************

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
    directory ``share/tuturials/benzene`` is:

    .. program-output:: h5dump -n gaussian_mbis.h5
       :cwd: ../share/tutorials/benzene

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

    qff.py --ffatypes=low --suffix=_noei gaussian.fchk

The script will dump all relevant information to the screen, for this tutorial,
the output is as follows:

.. program-output:: qff.py --ffatypes=low --suffix=_noei gaussian.fchk
   :cwd: ../share/tutorials/benzene

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
`pars_yaff_noei.txt`:

.. program-output:: cat pars_yaff_noei.txt
   :cwd: ../share/tutorials/benzene


Force field with electrostatics
===============================

* Generating Yaff input
    First we need to generate the Yaff input file for the electrostatic
    contribution to the force field. This is done using the script
    :ref:`qff_input-ei.py <seclab_rg_scripts_inputei>`. For this tutorial,
    we will convert the charges given in the dataset ``/charges`` of the
    file ``gaussian_mbis.h5`` for the atoms in gaussian.fchk with atom types
    according to the level `medium` and use Gaussian charge distributions::

        qff-input-ei.py  -v--ffatypes=low --gaussian gaussian.fchk gaussian_mbis.h5:charges pars_ei_mbisgauss.txt

    This command dumped the following output to the screen, indicating wheter or
    not the atom types are well chosen from the point of view of electrostatics
    (see second remark in :ref:`qff-input-ei.py <seclab_ug_tools_inputei>`):

    .. program-output:: qff-input-ei.py -v --ffatypes=low --gaussian gaussian.fchk gaussian_mbis.h5:charges pars_ei_mbisgauss.txt
       :cwd: ../share/tutorials/benzene

    Furthermore, the following Yaff parameter (`pars_ei.txt`) file was written:

    .. program-output:: cat pars_ei_mbisgauss.txt
       :cwd: ../share/tutorials/benzene

* Constructing the covalent contribution
    Now, we generate a covalent force field on top of the previously derived
    electrostatic contribution using the qff.py script::

        qff.py --ffatype=low --ei=pars_ei_mbisgauss.txt --suffix=_mbisgauss gaussian.fchk

    The logging output for this job is:

    .. program-output:: qff.py --ffatypes=low --ei=pars_ei_mbisgauss.txt --suffix=_mbisgauss gaussian.fchk
      :cwd: ../share/tutorials/benzene

    An extra line appeared in the beginning of the log output, i.e.
    ``QFF    Initializing Yaff force field reference for EI``. This indicates
    that an extra reference instance was created to represent the EI
    contribution to the force field. Furthermore, the covalent parameters are
    almost identical compared to the FF without electrostatics. This is indeed
    what we expect due to the charges being so small.

    The force field parameters were also written to the Yaff parameters file
    `pars_yaff_mbisgauss.txt`:

    .. program-output:: cat pars_yaff_mbisgauss.txt
       :cwd: ../share/tutorials/benzene

|

.. _seclab_tu_water:

Tutorial 2 - Water
********************

In the current tutorial, we will illustrate the usage of QuickFF as a library
by constructing a basic force field for water. In
:ref:`Tutorial 3 <seclab_tu_biphenyl>`, we will cover the more advanced features
of QuickFF.

For this tutorial we again assume to have the input files ``gaussian.fchk`` and
``gaussian_mbis.h5``. We will construct a script (which we will give the name
``qff-derive.py``) that will read the input files, convert them to data objects
QuickFF can read and implement a small program using QuickFF modules to derive
the force field. These files, together with the expected output, can be found
in ``share/tutorials/water``.

First, we define the required :ref:`settings <_seclab_ug_settings>` using the 
:ref:`Settings <_seclab_rg_modules_settings>` class:

.. literalinclude:: ../share/tutorials/water/qff-derive.py
   :language: python
   :lines: 28-30

Second, we load the input files and convert them to a Yaff
`system <http://molmod.github.io/yaff/rg_yaff.html#module-yaff.system>`_
instance and a QuickFF
`SecondOrderTaylor <http://molmod.github.io/QuickFF/rg_module.html#module-quickff.reference>`_
instance. Therefore, we use routines implemented in the
`fchk <http://molmod.github.io/molmod/reference/io.html#module-molmod.io.fchk>`_
module of MolMod.

.. literalinclude:: ../share/tutorials/water/qff-derive.py
   :language: python
   :lines: 32-47

At this instance, the system does not yet contain force field atom types. To
define these, we could use the `QuickFF feature guess_ffatypes <http://molmod.github.io/QuickFF/rg_module.html#quickff.tools.guess_ffatypes>`_.
However, we can also use the `ATSELECT <http://molmod.github.io/yaff/ug_atselect.html>`_
language of Yaff:

.. literalinclude:: ../share/tutorials/water/qff-derive.py
   :language: python
   :lines: 49-54

The final input we need to parse, is the set of atomic charges. These charges
are read from the dataset ``charges`` in the HDF5 file ``guassian_mbis.h5``.

.. literalinclude:: ../share/tutorials/water/qff-derive.py
   :language: python
   :lines: 56-61

``scales`` is a list of 4 floats, representing the electrostatic scalings of 1-2,
1-3, 1-4 and 1-5 pairs. In this case, no scaling is applied to the electrostatic
interactions. ``radii=None`` indicates that point charges are used.
``Average=True`` indicates that the set of charges are averaged over atoms of the same type.
Finally, ``pbc=[0,0,0]`` indicates that the system is not periodic. The returned
object ``ff_ei`` is an instance of the
`YaffForceField class <http://molmod.github.io/QuickFF/rg_module.html#quickff.reference.YaffForceField>`_.

Now we have all required input, and we will write a new Program class, which
will inherits from the
`BaseProgram class <http://molmod.github.io/QuickFF/rg_module.html#quickff.program.BaseProgram>`_,
and contain a sequence of instructions that will be executed consequently. In
this case, first a diagonal force field will be constructed (i.e. without any
cross terms), then we will add cross terms and fit the cross terms using a
Hessian least squares fit while keeping the diagonal terms fixed:

.. literalinclude:: ../share/tutorials/water/qff-derive.py
   :language: python
   :pyobject: Program

Finally, we can now setup the program and run it:

.. literalinclude:: ../share/tutorials/water/qff-derive.py
   :language: python
   :lines: 63-65

To summarize, the entire script as well as the logger output is given below.

.. container:: toggle

   .. container:: header

      **Click to show/hide the entire script**

   .. literalinclude:: ../share/tutorials/water/qff-derive.py
      :language: python

.. container:: toggle

   .. container:: header

      **Click to show/hide the logger output**

   .. program-output:: python qff-derive.py
      :cwd: ../share/tutorials/water/

|

.. _seclab_tu_biphenyl:

Tutorial 3 - Biphenyl
********************

In the current tutorial, we will illustrate advanced features of QuickFF by 
constructing a force field for biphenyl.

COMING SOON
