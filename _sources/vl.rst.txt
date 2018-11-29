.. _seclab_vl_main:

Validation
##########

Validation methodology
**********************

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

|
|

.. _seclab_vl_water:

Validation 1 - water
**********************

Force field construction
========================

.. container:: toggle

   .. container:: header

      **Click to show/hide the output of qff-input-ei.py**
      
   .. program-output:: qff-input-ei.py -v --ffatypes=low --gaussian gaussian.fchk gaussian_mbis.h5:charges pars_ei_mbisgauss.txt
      :cwd: ../share/validation/water


.. container:: toggle

   .. container:: header

      **Click to show/hide the Yaff parameter file of the electrostatic contribution**
        
   .. program-output:: cat pars_ei_mbisgauss.txt
      :cwd: ../share/validation/water


.. container:: toggle

   .. container:: header

      **Click to show/hide the QuickFF log file**
        
   .. program-output:: qff.py -c ../config.txt gaussian.fchk
      :cwd: ../share/validation/water


.. container:: toggle

   .. container:: header

      **Click to show/hide the Yaff parameter file of the covalent contribution**
        
   .. program-output:: cat pars_cov.txt
      :cwd: ../share/validation/water

|

Force field evaluation
======================

Perturbation trajectories
-------------------------

The figures below illustrate the reproduction of the *ab initio* energy along 
the constructed perturbation trajectories. For water, there are 3 such
trajectories, two bonds and one bend.


* **First O-H bond**:

.. image:: /_static/water/trajectory-BondHarm-H.O-0_Ehc3.png
    :width: 400
    :target: /_static/water/trajectory-BondHarm-H.O-0_Ehc3.png    

* **Second O-H bond**:

.. image:: /_static/water/trajectory-BondHarm-H.O-1_Ehc3.png
    :width: 400
    :target: /_static/water/trajectory-BondHarm-H.O-1_Ehc3.png
    
* **H-O-H bend**:

.. image:: /_static/water/trajectory-BendAHarm-H.O.H-2_Ehc3.png
    :width: 400
    :target: /_static/water/trajectory-BendAHarm-H.O.H-2_Ehc3.png


Geometry and frequencies
------------------------

.. program-output:: python simulate.py
  :cwd: ../share/validation/water

|

.. _seclab_vl__static/benzene:

Validation 2 - benzene
**********************

Force field construction
========================

.. container:: toggle

   .. container:: header

      **Click to show/hide the output of qff-input-ei.py**
      
   .. program-output:: qff-input-ei.py -v --ffatypes=low --gaussian gaussian.fchk gaussian_mbis.h5:charges pars_ei_mbisgauss.txt
      :cwd: ../share/validation/benzene


.. container:: toggle

   .. container:: header

      **Click to show/hide the Yaff parameter file of the electrostatic contribution**
        
   .. program-output:: cat pars_ei_mbisgauss.txt
      :cwd: ../share/validation/benzene


.. container:: toggle

   .. container:: header

      **Click to show/hide the QuickFF log file**
        
   .. program-output:: qff.py -c ../config.txt gaussian.fchk
      :cwd: ../share/validation/benzene


.. container:: toggle

   .. container:: header

      **Click to show/hide the Yaff parameter file of the covalent contribution**
        
   .. program-output:: cat pars_cov.txt
      :cwd: ../share/validation/benzene

|

Force field evaluation
======================

Perturbation trajectories
-------------------------

The figures below illustrate the reproduction of the *ab initio* energy along 
the constructed perturbation trajectories. For water, there are 3 such
trajectories, two bonds and one bend.


* **C-H bonds**:

.. image:: /_static/benzene/trajectory-BondHarm-C.H-0_Ehc3.png
    :width: 400
    :target: /_static/benzene/trajectory-BondHarm-C.H-0_Ehc3.png    

.. image:: /_static/benzene/trajectory-BondHarm-C.H-7_Ehc3.png
    :width: 400
    :target: /_static/benzene/trajectory-BondHarm-C.H-7_Ehc3.png

.. image:: ./_static/benzene/trajectory-BondHarm-C.H-8_Ehc3.png
    :width: 400
    :target: /_static/benzene/trajectory-BondHarm-C.H-8_Ehc3.png

.. image:: /_static/benzene/trajectory-BondHarm-C.H-9_Ehc3.png
    :width: 400
    :target: /_static/benzene/trajectory-BondHarm-C.H-9_Ehc3.png

.. image:: /_static/benzene/trajectory-BondHarm-C.H-10_Ehc3.png
    :width: 400
    :target: /_static/benzene/trajectory-BondHarm-C.H-10_Ehc3.png

.. image:: /_static/benzene/trajectory-BondHarm-C.H-11_Ehc3.png
    :width: 400
    :target: /_static/benzene/trajectory-BondHarm-C.H-11_Ehc3.png


* **C-C bonds**:

.. image:: /_static/benzene/trajectory-BondHarm-C.C-1_Ehc3.png
    :width: 400
    :target: /_static/benzene/trajectory-BondHarm-C.C-1_Ehc3.png

.. image:: /_static/benzene/trajectory-BondHarm-C.C-2_Ehc3.png
    :width: 400
    :target: /_static/benzene/trajectory-BondHarm-C.C-2_Ehc3.png

.. image:: /_static/benzene/trajectory-BondHarm-C.C-3_Ehc3.png
    :width: 400
    :target: /_static/benzene/trajectory-BondHarm-C.C-3_Ehc3.png

.. image:: /_static/benzene/trajectory-BondHarm-C.C-4_Ehc3.png
    :width: 400
    :target: /_static/benzene/trajectory-BondHarm-C.C-4_Ehc3.png

.. image:: /_static/benzene/trajectory-BondHarm-C.C-5_Ehc3.png
    :width: 400
    :target: /_static/benzene/trajectory-BondHarm-C.C-5_Ehc3.png

.. image:: /_static/benzene/trajectory-BondHarm-C.C-6_Ehc3.png
    :width: 400
    :target: /_static/benzene/trajectory-BondHarm-C.C-6_Ehc3.png   


* **C-C-C bends**:

.. image:: /_static/benzene/trajectory-BendAHarm-C.C.C-24_Ehc3.png
    :width: 400
    :target: /_static/benzene/trajectory-BendAHarm-C.C.C-24_Ehc3.png

.. image:: /_static/benzene/trajectory-BendAHarm-C.C.C-25_Ehc3.png
    :width: 400
    :target: /_static/benzene/trajectory-BendAHarm-C.C.C-25_Ehc3.png

.. image:: /_static/benzene/trajectory-BendAHarm-C.C.C-26_Ehc3.png
    :width: 400
    :target: /_static/benzene/trajectory-BendAHarm-C.C.C-26_Ehc3.png

.. image:: /_static/benzene/trajectory-BendAHarm-C.C.C-27_Ehc3.png
    :width: 400
    :target: /_static/benzene/trajectory-BendAHarm-C.C.C-27_Ehc3.png

.. image:: /_static/benzene/trajectory-BendAHarm-C.C.C-28_Ehc3.png
    :width: 400
    :target: /_static/benzene/trajectory-BendAHarm-C.C.C-28_Ehc3.png

.. image:: /_static/benzene/trajectory-BendAHarm-C.C.C-29_Ehc3.png
    :width: 400
    :target: /_static/benzene/trajectory-BendAHarm-C.C.C-29_Ehc3.png


* **C-C-H bends**:

.. image:: /_static/benzene/trajectory-BendAHarm-C.C.H-12_Ehc3.png
    :width: 400
    :target: /_static/benzene/trajectory-BendAHarm-C.C.H-12_Ehc3.png

.. image:: /_static/benzene/trajectory-BendAHarm-C.C.H-13_Ehc3.png
    :width: 400
    :target: /_static/benzene/trajectory-BendAHarm-C.C.H-13_Ehc3.png

.. image:: /_static/benzene/trajectory-BendAHarm-C.C.H-14_Ehc3.png
    :width: 400
    :target: /_static/benzene/trajectory-BendAHarm-C.C.H-14_Ehc3.png

.. image:: /_static/benzene/trajectory-BendAHarm-C.C.H-15_Ehc3.png
    :width: 400
    :target: /_static/benzene/trajectory-BendAHarm-C.C.H-15_Ehc3.png

.. image:: /_static/benzene/trajectory-BendAHarm-C.C.H-16_Ehc3.png
    :width: 400
    :target: /_static/benzene/trajectory-BendAHarm-C.C.H-16_Ehc3.png

.. image:: /_static/benzene/trajectory-BendAHarm-C.C.H-17_Ehc3.png
    :width: 400
    :target: /_static/benzene/trajectory-BendAHarm-C.C.H-17_Ehc3.png


* **C-C-H-C out-of-plane distances**:

.. image:: /_static/benzene/trajectory-Oopdist-C.C.H.C-54_Ehc3.png
    :width: 400
    :target: /_static/benzene/trajectory-Oopdist-C.C.H.C-54_Ehc3.png

.. image:: /_static/benzene/trajectory-Oopdist-C.C.H.C-55_Ehc3.png
    :width: 400
    :target: /_static/benzene/trajectory-Oopdist-C.C.H.C-55_Ehc3.png

.. image:: /_static/benzene/trajectory-Oopdist-C.C.H.C-56_Ehc3.png
    :width: 400
    :target: /_static/benzene/trajectory-Oopdist-C.C.H.C-56_Ehc3.png

.. image:: /_static/benzene/trajectory-Oopdist-C.C.H.C-57_Ehc3.png
    :width: 400
    :target: /_static/benzene/trajectory-Oopdist-C.C.H.C-57_Ehc3.png

.. image:: /_static/benzene/trajectory-Oopdist-C.C.H.C-58_Ehc3.png
    :width: 400
    :target: /_static/benzene/trajectory-Oopdist-C.C.H.C-58_Ehc3.png

.. image:: /_static/benzene/trajectory-Oopdist-C.C.H.C-59_Ehc3.png
    :width: 400
    :target: /_static/benzene/trajectory-Oopdist-C.C.H.C-59_Ehc3.png


Geometry and frequencies
------------------------

.. program-output:: python simulate.py
  :cwd: ../share/validation/benzene
