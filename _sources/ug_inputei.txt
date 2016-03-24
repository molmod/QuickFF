.. _seclab_ug_tools:

.. _seclab_ug_tools_inputei:

Generating Yaff parameter file for the electrostatic contribution
##################################################################

A Yaff force field file for the electrostatic contribution can be constructed
from the output of `Horton <https://theochem.github.io/horton/>`_ using the 
:ref:`qff-input-ei.py <seclab_rg_scripts_inputei>` script::

    qff-input-ei.py [options] fn_sys fn_wpart path

**Arguments**

* ``fn_sys``
    a file containing the system information. Any file type accepted by the 
    main script :ref:`qff.py <seclab_rg_scripts_qff>` is also accepted here
    (see :ref:`input files <seclab_inputfiles>`). 

* ``fn_wpart``
    is a `Horton <https://theochem.github.io/horton/>`_: HDF5 file containing
    the charges of each atom. The charges in this file should have the same 
    ordering as the atoms in the system file `fn_sys`.

* ``path``
    is the path in the HDF5 file which contains a dataset `charges`. In other
    words, the file `fn_wpart` should contain a dataset `/path/charges`
    containing the charges for each atom in the system.

**Options**

* Force field atom types (:option:`--ffatypes=LEVEL`):
    Determine the force field atom types according to the given level (see
    :ref:`automatic estimation of atom types <seclab_ug_atype_estimator>`). By
    default, the atom types are assumed to be defined in `fn_sys`.

* Gaussian distributed charges (:option:`--gaussian`):
    Treat the atomic charges as Gaussian charge distributions with atomic radii
    determined according to the precedure of 
    `Chen et al. <http://www.sciencedirect.com/science/article/pii/S0009261407002618>`_,
    implemented in the routine `get_ei_radii` from the
    :ref:`tools <seclab_rg_modules_tools>` module.

**Remarks**

All electrostatic interactions will be included in the force field, i.e. no
scaling is applied. This means that atoms that are connected to each other 1, 2 
or 3 bonds will also interact electrostatically. To change this behaviour, open
the resulting Yaff parameter file after applying this script, and change the
corresponding scales in the file. For example: `FIXQ:SCALE 1 0.0` will turn of 
electrostatic interactions for bonded atoms. See also the 
`Yaff page <http://molmod.github.io/yaff/ug_forcefield.html#prefix-fixq>`_ on 
the parameter file for electrostatic interactions.

The script might print a WARNING to the screen indicating the for a certain atom
type the standard deviation of the charges for this atom type is rather high.
This could indicate that the current atom types are to general, i.e. not
specific enough for the current chemical configuration. A higher level of atom 
type specification might then be adequate.
