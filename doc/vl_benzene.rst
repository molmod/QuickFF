.. _seclab_vl_benzene:

Validation 2 - Benzene
######################

Force field construction
************************

.. container:: toggle

   .. container:: header

      **Click to show/hide the output of qff-input-ei.py**
      
   .. program-output:: qff-input-ei.py -v --ffatypes=low --gaussian gaussian.fchk gaussian_mbis.h5:charges pars_ei_mbisgauss.txt
      :cwd: ../share/validation/benzene

|

.. container:: toggle

   .. container:: header

      **Click to show/hide the Yaff parameter file of the electrostatic contribution**
        
   .. program-output:: cat pars_ei_mbisgauss.txt
      :cwd: ../share/validation/benzene

|

.. container:: toggle

   .. container:: header

      **Click to show/hide the QuickFF log file**
        
   .. program-output:: qff.py -c ../config.txt gaussian.fchk
      :cwd: ../share/validation/benzene

|

.. container:: toggle

   .. container:: header

      **Click to show/hide the Yaff parameter file of the covalent contribution**
        
   .. program-output:: cat pars_cov.txt
      :cwd: ../share/validation/benzene

Force field evaluation
**********************

Perturbation trajectories
=========================

The figures below illustrate the reproduction of the *ab initio* energy along 
the constructed perturbation trajectories. For water, there are 3 such
trajectories, two bonds and one bend.


* **C-H bonds**:

.. image:: ../../../share/validation/benzene/trajectory-BondHarm-C.H-0_Ehc3.png
    :width: 400
    :target: ../../../share/validation/benzene/trajectory-BondHarm-C.H-0_Ehc3.png    

.. image:: ../../../share/validation/benzene/trajectory-BondHarm-C.H-7_Ehc3.png
    :width: 400
    :target: ../../../share/validation/benzene/trajectory-BondHarm-C.H-7_Ehc3.png

.. image:: ../../../share/validation/benzene/trajectory-BondHarm-C.H-8_Ehc3.png
    :width: 400
    :target: ../../../share/validation/benzene/trajectory-BondHarm-C.H-8_Ehc3.png

.. image:: ../../../share/validation/benzene/trajectory-BondHarm-C.H-9_Ehc3.png
    :width: 400
    :target: ../../../share/validation/benzene/trajectory-BondHarm-C.H-9_Ehc3.png

.. image:: ../../../share/validation/benzene/trajectory-BondHarm-C.H-10_Ehc3.png
    :width: 400
    :target: ../../../share/validation/benzene/trajectory-BondHarm-C.H-10_Ehc3.png

.. image:: ../../../share/validation/benzene/trajectory-BondHarm-C.H-11_Ehc3.png
    :width: 400
    :target: ../../../share/validation/benzene/trajectory-BondHarm-C.H-11_Ehc3.png


* **C-C bonds**:

.. image:: ../../../share/validation/benzene/trajectory-BondHarm-C.C-1_Ehc3.png
    :width: 400
    :target: ../../../share/validation/benzene/trajectory-BondHarm-C.C-1_Ehc3.png

.. image:: ../../../share/validation/benzene/trajectory-BondHarm-C.C-2_Ehc3.png
    :width: 400
    :target: ../../../share/validation/benzene/trajectory-BondHarm-C.C-2_Ehc3.png

.. image:: ../../../share/validation/benzene/trajectory-BondHarm-C.C-3_Ehc3.png
    :width: 400
    :target: ../../../share/validation/benzene/trajectory-BondHarm-C.C-3_Ehc3.png

.. image:: ../../../share/validation/benzene/trajectory-BondHarm-C.C-4_Ehc3.png
    :width: 400
    :target: ../../../share/validation/benzene/trajectory-BondHarm-C.C-4_Ehc3.png

.. image:: ../../../share/validation/benzene/trajectory-BondHarm-C.C-5_Ehc3.png
    :width: 400
    :target: ../../../share/validation/benzene/trajectory-BondHarm-C.C-5_Ehc3.png

.. image:: ../../../share/validation/benzene/trajectory-BondHarm-C.C-6_Ehc3.png
    :width: 400
    :target: ../../../share/validation/benzene/trajectory-BondHarm-C.C-6_Ehc3.png   


* **C-C-C bends**:

.. image:: ../../../share/validation/benzene/trajectory-BendAHarm-C.C.C-24_Ehc3.png
    :width: 400
    :target: ../../../share/validation/benzene/trajectory-BendAHarm-C.C.C-24_Ehc3.png

.. image:: ../../../share/validation/benzene/trajectory-BendAHarm-C.C.C-25_Ehc3.png
    :width: 400
    :target: ../../../share/validation/benzene/trajectory-BendAHarm-C.C.C-25_Ehc3.png

.. image:: ../../../share/validation/benzene/trajectory-BendAHarm-C.C.C-26_Ehc3.png
    :width: 400
    :target: ../../../share/validation/benzene/trajectory-BendAHarm-C.C.C-26_Ehc3.png

.. image:: ../../../share/validation/benzene/trajectory-BendAHarm-C.C.C-27_Ehc3.png
    :width: 400
    :target: ../../../share/validation/benzene/trajectory-BendAHarm-C.C.C-27_Ehc3.png

.. image:: ../../../share/validation/benzene/trajectory-BendAHarm-C.C.C-28_Ehc3.png
    :width: 400
    :target: ../../../share/validation/benzene/trajectory-BendAHarm-C.C.C-28_Ehc3.png

.. image:: ../../../share/validation/benzene/trajectory-BendAHarm-C.C.C-29_Ehc3.png
    :width: 400
    :target: ../../../share/validation/benzene/trajectory-BendAHarm-C.C.C-29_Ehc3.png


* **C-C-H bends**:

.. image:: ../../../share/validation/benzene/trajectory-BendAHarm-C.C.H-12_Ehc3.png
    :width: 400
    :target: ../../../share/validation/benzene/trajectory-BendAHarm-C.C.H-12_Ehc3.png

.. image:: ../../../share/validation/benzene/trajectory-BendAHarm-C.C.H-13_Ehc3.png
    :width: 400
    :target: ../../../share/validation/benzene/trajectory-BendAHarm-C.C.H-13_Ehc3.png

.. image:: ../../../share/validation/benzene/trajectory-BendAHarm-C.C.H-14_Ehc3.png
    :width: 400
    :target: ../../../share/validation/benzene/trajectory-BendAHarm-C.C.H-14_Ehc3.png

.. image:: ../../../share/validation/benzene/trajectory-BendAHarm-C.C.H-15_Ehc3.png
    :width: 400
    :target: ../../../share/validation/benzene/trajectory-BendAHarm-C.C.H-15_Ehc3.png

.. image:: ../../../share/validation/benzene/trajectory-BendAHarm-C.C.H-16_Ehc3.png
    :width: 400
    :target: ../../../share/validation/benzene/trajectory-BendAHarm-C.C.H-16_Ehc3.png

.. image:: ../../../share/validation/benzene/trajectory-BendAHarm-C.C.H-17_Ehc3.png
    :width: 400
    :target: ../../../share/validation/benzene/trajectory-BendAHarm-C.C.H-17_Ehc3.png


* **C-C-H-C out-of-plane distances**:

.. image:: ../../../share/validation/benzene/trajectory-Oopdist-C.C.H.C-54_Ehc3.png
    :width: 400
    :target: ../../../share/validation/benzene/trajectory-Oopdist-C.C.H.C-54_Ehc3.png

.. image:: ../../../share/validation/benzene/trajectory-Oopdist-C.C.H.C-55_Ehc3.png
    :width: 400
    :target: ../../../share/validation/benzene/trajectory-Oopdist-C.C.H.C-55_Ehc3.png

.. image:: ../../../share/validation/benzene/trajectory-Oopdist-C.C.H.C-56_Ehc3.png
    :width: 400
    :target: ../../../share/validation/benzene/trajectory-Oopdist-C.C.H.C-56_Ehc3.png

.. image:: ../../../share/validation/benzene/trajectory-Oopdist-C.C.H.C-57_Ehc3.png
    :width: 400
    :target: ../../../share/validation/benzene/trajectory-Oopdist-C.C.H.C-57_Ehc3.png

.. image:: ../../../share/validation/benzene/trajectory-Oopdist-C.C.H.C-58_Ehc3.png
    :width: 400
    :target: ../../../share/validation/benzene/trajectory-Oopdist-C.C.H.C-58_Ehc3.png

.. image:: ../../../share/validation/benzene/trajectory-Oopdist-C.C.H.C-59_Ehc3.png
    :width: 400
    :target: ../../../share/validation/benzene/trajectory-Oopdist-C.C.H.C-59_Ehc3.png


Geometry and frequencies
========================

.. program-output:: python simulate.py
  :cwd: ../share/validation/benzene
