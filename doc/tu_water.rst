.. _seclab_tu_water:

Tutorial 2 - Water
###################

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

First we load the input files and convert them to a Yaff
`system <http://molmod.github.io/yaff/rg_yaff.html#module-yaff.system>`_
instance and a QuickFF
`SecondOrderTaylor <http://molmod.github.io/QuickFF/rg_module.html#module-quickff.reference>`_
instance. Therefore, we use routines implemented in the
`fchk <http://molmod.github.io/molmod/reference/io.html#module-molmod.io.fchk>`_
module of MolMod.

.. literalinclude:: ../share/tutorials/water/qff-derive.py
   :language: python
   :lines: 25-39

At this instance, the system does not yet contain force field atom types. To
define these, we could use the `QuickFF feature guess_ffatypes <http://molmod.github.io/QuickFF/rg_module.html#quickff.tools.guess_ffatypes>`_.
However, we can also use the `ATSELECT <http://molmod.github.io/yaff/ug_atselect.html>`_
language of Yaff:

.. literalinclude:: ../share/tutorials/water/qff-derive.py
   :language: python
   :lines: 41-46

The final input we need to parse, is the set of atomic charges. These charges
are read from the dataset ``charges`` in the HDF5 file ``guassian_mbis.h5``.

.. literalinclude:: ../share/tutorials/water/qff-derive.py
   :language: python
   :lines: 48-52

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
   :lines: 54-56

To summarize, the entire script as well as the logger output is given below.

.. container:: toggle

   .. container:: header

      **Click to show/hide the entire script**

   .. literalinclude:: ../share/tutorials/water/qff-derive.py
      :language: python

.. container:: toggle

   .. container:: header

      **Click to show/hide the logger output**

   .. program-output:: python /home/louis/build/quickff/share/tutorials/water/qff-derive.py
