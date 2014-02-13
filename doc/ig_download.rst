Download QuickFF
################

Stable release
==============

A stable version of QuickFF, which includes all features described on this 
website, can be downloaded here:

    http://users.ugent.be/~lvduyfhu/quickff/quickff-1.0.tar.gz.

Choose a suitable directory, e.g. ``~/build``, download and unpack the archive::

    mkdir -p ~/build
    cd ~/build
    wget http://users.ugent.be/~lvduyfhu/quickff/quickff-1.0.tar.gz
    tar -xvzf quickff-1.0.tar.gz
    cd quickff-1.0


Latest development version (experts only)
=========================================

The latest development version of QuickFF can only be downloaded using Git.
This also allows you to upload your own changes to QuickFF. Git is free and
open-source distributed revision control system to easily handle programming
projects shared between several people. Further information about git (including
downloads and tutorials) can be found `here <http://git-scm.com/>`_. The
official git URL for QuickFF is `http://github.com/molmod/QuickFF/ 
<http://github.com/molmod/QuickFF>`_. To `clone` the QuickFF repository, go to 
your favorite directory for source code, e.g. ``~/build``, and execute the 
following commands::

    git clone git://github.com/molmod/QuickFF/quickff.git
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
