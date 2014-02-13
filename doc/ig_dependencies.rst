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
<https://nose.readthedocs.org/en/latest/>`_ package, and to build the 
documentation of this website yourself, the `Spinx <http://sphinx-doc.org/>`_ 
package needs to be installed. Most Linux distributions can install this 
software with just a single terminal command.

* Ubuntu 12.4 and later::

    sudo apt-get install python-numpy python-h5py python-matplotlib python-nose python-sphinx

* Fedora 17 and later::

    sudo yum install numpy h5py python-matplotlib python-nose sphinx

Finally, to be able to run QuickFF in parallel on multiple cores, the Python
package `Scoop <https://code.google.com/p/scoop/>`_ needs to be installed. To
download and install Scoop (version 0.6.0.final), run the following commands

* Ubuntu 12.4 and later::

    wget 'https://scoop.googlecode.com/files/scoop-0.6.0.final.tar.gz'
    tar -xvzf scoop-0.6.0.final.tar.gz
    cd scoop-0.6.0.final
    python setup.py install
