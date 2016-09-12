Dependencies
############

MolMod dependency
*****************

`MolMod <http://molmod.github.com/molmod/>`_ is a Python library used by most
Python programs developed at the CMM. It must be installed before QuickFF can
be used or installed. Installation and download instructions can be found in the
`molmod documentation <http://molmod.github.com/molmod/tutorial/install.html>`_.
The instructions below only work if the MolMod package is installed.

Yaff dependency
***************

`Yaff <http://molmod.github.com/yaff/>`_ is a pythonic force-field (FF)
code used at the Center for Molecular Modeling (CMM) to test-drive new FF models.
Yaff is required to run QuickFF, the implementation of each force field term is
based on the implementation in Yaff. For this version of QuickFF, it is required
to install the 
`1.0.develop.2.13 <https://github.com/molmod/yaff/tree/1.0.develop.2.13>`_ 
branch of Yaff available at github. It can be installed by cloning the git 
repository using the git command::

    git clone https://github.com/molmod/yaff.git

and switch to the branch *1.0.develop.2.13*. Alternatively, a zip file can be
downloaded from the online repository for the branch
`1.0.develop.2.13 <https://github.com/molmod/yaff/tree/1.0.develop.2.13>`_.
After extracting the zip file, Yaff can be installed by following the 
instructions given in the 
`yaff documentation <http://molmod.github.io/yaff/ug_install.html>`_. The 
instructions below only work if the Yaff package is installed.

External dependencies
*********************

Some software packages should be installed before QuickFF can be installed or
used. It is recommended to use the software package management of your Linux
distribution to install these dependencies.

The following software must be installed:

* Python >= 2.7 (including the development files): http://www.python.org/
* Numpy >= 1.0: http://numpy.scipy.org/
* Scipy >= 0.14: http://www.scipy.org/
* h5py >= 2.0.0: http://code.google.com/p/h5py/
* matplotlib >= 1.0.0: http://matplotlib.sourceforge.net/

To run the python tests, it is necessary to install the `Nose 
<https://nose.readthedocs.org/en/latest/>`_ package. Most Linux distributions 
can install this software with just a single terminal command.

* Ubuntu 12.4 and later::

    sudo apt-get install python-numpy python-scipy python-h5py python-matplotlib python-nose

* Fedora 17 and later::

    sudo yum install numpy scipy h5py python-matplotlib python-nose

Finally, to be able to run QuickFF in parallel on multiple cores, the Python
package `Scoop <https://code.google.com/p/scoop/>`_ needs to be installed. To
download and install Scoop (version 0.6.0.final), run the following commands

* Ubuntu 12.4 and later::

    wget 'https://scoop.googlecode.com/files/scoop-0.6.0.final.tar.gz'
    tar -xvzf scoop-0.6.0.final.tar.gz
    cd scoop-0.6.0.final
    python setup.py install


Downloading QuickFF
###################

Stable release
**************

A stable release of the second generation of QuickFF, which includes all 
features described on this website, can be downloaded from:

    http://users.ugent.be/~lvduyfhu/QuickFF-2.1.2.tar.gz

Previous versions available for download can be found below:

    http://users.ugent.be/~lvduyfhu/QuickFF-2.0.1.tar.gz
    
    http://users.ugent.be/~lvduyfhu/quickff-1.0.1.tar.gz

Choose a suitable directory, e.g. ``~/build``, download and unpack the archive::

    mkdir -p ~/build
    cd ~/build
    wget http://users.ugent.be/~lvduyfhu/QuickFF-2.1.2.tar.gz
    tar -xvzf QuickFF-2.1.2.tar.gz
    cd QuickFF-2.1.2


Latest development version (experts only)
*****************************************

The latest development version of QuickFF can only be downloaded using Git.
This also allows you to upload your own changes to QuickFF. Git is free and
open-source distributed revision control system to easily handle programming
projects shared between several people. Further information about git (including
downloads and tutorials) can be found `here <http://git-scm.com/>`_. The
official git URL for QuickFF is http://github.com/molmod/QuickFF . To `clone` 
the QuickFF repository, go to your favorite directory for source code, e.g. 
``~/build``, and execute the following commands::

    git clone git://github.com/molmod/QuickFF.git
    cd quickff

The source code can be updated with the latest patches with the following
command::

    git pull

This will also update the version history so that the progress can easily be
tracked. This version history can be visualized using the `gitk` program. This
program can be downloaded with the following command (Ubuntu)::

    sudo apt-get install gitk

Once, `gitk` is installed, the version history can be visualized by going to the
directory containing the quickff repository and running the command::

    gitk


Installing QuickFF
##################

Once you downloaded the source code and installed all required packages, QuickFF
can be installed. To install QuickFF, simply go to the directory in which you
downloaded the QuickFF source directory and run the following command::

    python setup.py install

After the installation is complete, you can test your installation by running
the following commands::

    cd test
    nosetests

Once all tests are succesfull, you are ready to use QuickFF.
