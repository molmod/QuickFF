.. _seclab_vl_main:

Validation methodology
######################

In this section, a force field is constructed using the current version
of QuickFF (i.e. |release|), after which it is tested and evaluated for its
ability to reproduce the *ab initio* input it was fitted to. The force field was
constructed using the QuickFF script :ref:`qff.py <seclab_ug_script>`. The
required *ab initio* input is contained in the file ``gaussian.fchk`` and the
electrostatic partitioning in ``gaussian_mbis.h5``. The Yaff parameter file for
the electrostatic contribution is generated using

    qff-input-ei.py -f --ffatypes=low --gaussian gaussian.fchk gaussian_mbis.h5:charges pars_ei_mbisgauss.txt

The config file containing the settings that are used for the force field is 
given by ``config.txt`` and can be viewed below. Furthermore, all required input
files as well as the generated output can be found in the directory 
``share/validation``.

.. container:: toggle

   .. container:: header

      **Click to show/hide the QuickFF config file used for generating the force field**
        
   .. include:: ../share/validation/config.txt
      :literal:

To evaluate the performance of the force field, it is applied to perform a basic
geometry optimization followed by a normal mode analysis using Yaff and TAMkin.
Finally, two validations are performed for each system:

* the force field energy is compared to the *ab initio* energy along the 
  perturbation trajectories that were constructed by QuickFF.

* the geometry (i.e. bond lengths, bending angles, dihedral angles and 
  ut-of-plane distances) and normal mode frequencies as predicted by the force 
  field are compared with the given *ab initio* input. To this end, the
  following statistics are used:
  
  - Root Mean Square Deviation (RMSD)

    .. math:: 
       RMSD = \sqrt{\frac{1}{N}\sum_{i=1}^N \left(A_i - F_i\right)^2}

  - Mean Deviation (RMSD)

    .. math:: 
       MD = \frac{1}{N}\sum_{i=1}^N \left(A_i-F_i\right)

  - Root Mean Square Error (RMSE)

    .. math:: 
       RMSE = \sqrt{\frac{1}{N}\sum_{i=1}^N \left(A_i-F_i-MD\right)^2}

  These statistics satisfy the condition :math:`RMDS^2=RMSE^2+MD^2`. 

Validated systems
#################

The following systems are considered for these validation tests:

.. toctree::
   :maxdepth: 1

   vl_water.rst
   vl_benzene.rst
