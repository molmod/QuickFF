Dependencies
############

This section describes the dependency libraries of QuickFF.

MolMod dependency
=================

`MolMod <http://molmod.github.com/molmod/>`_ is a Python library used by most
Python programs developed at the CMM. It must be installed before QuickFF can
be used or installed. Installation and download instructions can be found in the
`molmod documentation <http://molmod.github.com/molmod/tutorial/install.html>`_.
The instructions below only work if the MolMod package is installed.

Yaff dependency
=================

`Yaff <http://molmod.github.com/yaff/>`_ is a pythonic force-field (FF)
code used at the Center for Molecular Modeling (CMM) to test-drive new FF models.
Installation and download instructions can be found in the
`yaff documentation <http://molmod.github.io/yaff/ug_install.html>`_.
The instructions below only work if the Yaff package is installed.

External dependencies
=====================

Some software packages should be installed before QuickFF can be installed or
used. It is recommended to use the software package management of your Linux
distribution to install these dependencies.

The following software must be installed:

* Python >= 2.7 (including the development files): http://www.python.org/
* Numpy >= 1.0: http://numpy.scipy.org/
* h5py >= 2.0.0: http://code.google.com/p/h5py/
* matplotlib >= 1.0.0: http://matplotlib.sourceforge.net/

To run the python tests, it is necessary to install the `Nose 
<https://nose.readthedocs.org/en/latest/>`_ package. Most Linux distributions 
can install this software with just a single terminal command.

* Ubuntu 12.4 and later::

    sudo apt-get install python-numpy python-h5py python-matplotlib python-nose

* Fedora 17 and later::

    sudo yum install numpy h5py python-matplotlib python-nose

Finally, to be able to run QuickFF in parallel on multiple cores, the Python
package `Scoop <https://code.google.com/p/scoop/>`_ needs to be installed. To
download and install Scoop (version 0.6.0.final), run the following commands

* Ubuntu 12.4 and later::

    wget 'https://scoop.googlecode.com/files/scoop-0.6.0.final.tar.gz'
    tar -xvzf scoop-0.6.0.final.tar.gz
    cd scoop-0.6.0.final
    python setup.py install

An optional dependency is ALGLIB, which is used to perform the minimization of
a quadratic cost function with box constraints. If ALGLIB is not installed,
this minimization can also be performed using standar SciPy minimizers.
You can find instructions on how to install the Python wrapper for ALGLIB
`here <http://www.alglib.net/download.php>`_. Following commands might work:

* UNIX::

    wget http://www.alglib.net/translator/re/alglib-3.10.0.cpython.gpl.zip
    unzip alglib-3.10.0.cpython.gpl.zip
    cd cpython/core
    chmod u+x aenv build-unix
    ./build-unix
    cd ../
    chmod u+x setup.py
    ./setup.py install --home=~
