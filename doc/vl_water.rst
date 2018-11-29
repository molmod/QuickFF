.. _seclab_vl_water:

Validation 1 - Water
###################

Force field construction
************************

.. container:: toggle

   .. container:: header

      **Click to show/hide the output of qff-input-ei.py**
      
   .. program-output:: qff-input-ei.py -v --ffatypes=low --gaussian gaussian.fchk gaussian_mbis.h5:charges pars_ei_mbisgauss.txt
      :cwd: ../share/validation/water

|

.. container:: toggle

   .. container:: header

      **Click to show/hide the Yaff parameter file of the electrostatic contribution**
        
   .. program-output:: cat pars_ei_mbisgauss.txt
      :cwd: ../share/validation/water

|

.. container:: toggle

   .. container:: header

      **Click to show/hide the QuickFF log file**
        
   .. program-output:: qff.py -c ../config.txt gaussian.fchk
      :cwd: ../share/validation/water

|

.. container:: toggle

   .. container:: header

      **Click to show/hide the Yaff parameter file of the covalent contribution**
        
   .. program-output:: cat pars_cov.txt
      :cwd: ../share/validation/water

Force field evaluation
**********************

Perturbation trajectories
=========================

The figures below illustrate the reproduction of the *ab initio* energy along 
the constructed perturbation trajectories. For water, there are 3 such
trajectories, two bonds and one bend.


* **First O-H bond**:

.. image:: ../../../share/validation/water/trajectory-BondHarm-H.O-0_Ehc3.png
    :width: 400
    :target: ../../../share/validation/water/trajectory-BondHarm-H.O-0_Ehc3.png    

* **Second O-H bond**:

.. image:: ../../../share/validation/water/trajectory-BondHarm-H.O-1_Ehc3.png
    :width: 400
    :target: ../../../share/validation/water/trajectory-BondHarm-H.O-1_Ehc3.png
    
* **H-O-H bend**:

.. image:: ../../../share/validation/water/trajectory-BendAHarm-H.O.H-2_Ehc3.png
    :width: 400
    :target: ../../../share/validation/water/trajectory-BendAHarm-H.O.H-2_Ehc3.png


Geometry and frequencies
========================

.. program-output:: python simulate.py
  :cwd: ../share/validation/water
