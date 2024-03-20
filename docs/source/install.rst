.. include:: numbatdefs.txt

.. _chap-install-label:

*******************
Installing |NUMBAT|
*******************

This chapter provides instructions on installing |NUMBAT| on each platform.
Please email |NUMBAT_EMAIL| to let us know of any difficulties you encounter, or suggestions for improvements to the install procedure on any platform, but especially for MacOS or Windows.

Information for all platforms
================================
While |NUMBAT| is developed on Linux, it can also be built on MacOS X as a
native command-line application, and under Windows using a virtual machine running Linux.

In all cases, the current source code for |NUMBAT| is hosted `here on Github <https://github.com/michaeljsteel/|NUMBAT|>`_. Please always download the latest release from the github page.


Install locations
---------------------

There is no need to install |NUMBAT| in a central location
such as ``/usr/local/`` or ``/opt/local/`` though you may certainly choose to do so.


Here and throughout this documentation, we use the string ``<NumBAT>`` to indicate the root |NUMBAT| install directory (e.g. ``/usr/local/NumBAT``, ``/home/mike/NumBAT``, ``/home/myuserid/research/NumBAT``).



Installing on Linux
================================


|NUMBAT| has been developed and tested on Ubuntu 23.04 with the following
package versions: Python 3.11.4, Numpy 1.24.2, Arpack-NG, Suitesparse 7.1.0,
and Gmsh 4.8.4.  |NUMBAT| also depends on the BLAS and LAPACK libraries. We
strongly recommend linking |NUMBAT| against optimised versions, such as the MKL
library provided in the  free Intel OneAPI library.

|NUMBAT| has also been successfully installed by users on Debian and
RedHat/Fedora, and with different versions of packages, but these installations
have not been as thoroughly documented so may require user testing.  In
general, any relatively current Linux system should work without trouble.

|NUMBAT| building and installation is easiest if you have root access, but it
is not required.  See the section below if you do not have root access (or the
ability to run ``sudo``) on your machine.

The following steps use package syntax for Ubuntu/Debian systems. For other
Linux flavours, you may need to use different package manager syntax and/or
slightly different package names.

The code depends critically on matrix-vector operations provided by Lapack and
Blas libraries. We strongly recommend using an optimised library such as the
Intel OneAPI library (for Intel CPUs) or the AMD Optimizing CPU Libraries
(AOCL) for AMD CPUs.  The steps below demonstrate the Intel OneAPI approach.

#. Before installing, ensure your system is up to date ::

    $ sudo apt-get update
    $ sudo apt-get upgrade

#. Create a python virtual environment for working with |NUMBAT|.
   You can use any name and location for your tree. 
   To specify a virtual environment tree called `nbpy3` in your home directory, enter ::

    $ cd ~
    $ python3 -m virtualenv nbpy3


#. Activate the new python virtual environment ::

    $ source ~/nbpy3/bin/activate

#. Install necessary python libraries ::

    $ pip3 install numpy matplotlib scipy psutils

#. Create a working directory for your |NUMBAT| work and move into it.

#. To download the current version from the git repository and install any missing library dependencies, use ::

    $ git clone https://github.com/michaeljsteel/NumBAT.git
    $ cd NumBAT


#. Open the file ``<NumBAT>/backend/fortran/Makefile`` in a text editor and check the settings associated with the variables ``PLAT`` that control the preferred math library. The default setting is to use the Intel OneAPI library which is the recommended configuration.

#. Now at last, we can build |NUMBAT| by running the following in the root ``<NumBAT>`` directory. ::

   $ make build

#. If this completes without error, you are ready to proceed to the next chapter to begin using |NUMBAT|.

#. If you hit a compile error you can't resolve, please get in touch at |NUMBAT_EMAIL|.


Other build configurations
--------------------------

The Fortran components (|NUMBAT| source code and libraries) have been successfully compiled with Intel's ``ifortran`` as well as GCC's open-source ``gfortran``. In this documentation we use ``gfortran``, but this can be easily adjusted in ``|NUMBAT|/backend/fortran/Makefile``


Other build configurations
--------------------------

The Fortran components (NumBAT source code and libraries) have been successfully compiled with Intel's ``ifortran`` as well as GCC's open-source ``gfortran``. In this documentation we use ``gfortran``, but this can be easily adjusted in ``NumBAT/backend/fortran/Makefile``

Installing without root access
----------------------------------------------------
Compiling and installing |NUMBAT| itself does not rely on having root access to your machine. However, installing the supporting libraries such as SuiteSparse and Arpack
is certainly simpler if you have root or the assistance of your system admin.

If this is not possible, you can certainly proceed by building and installing all the required libraries into your own tree within your home directory.
It may be helpful to create a tree like the following so that the relevant paths mirror those of a system install::

   $HOME/
    |---my_sys/
         |---usr/
              |---include/
              |---lib/
              |---bin/


Installing on MacOS
================================
|NUMBAT| can also be installed on MacOS, though this is currently somewhat experimental and has only been performed on certain versions of MacOS.  Any comments on difficulties and solutions will be appreciated.

The following steps have worked for us:

#. Open a terminal window on your desktop.

#. Ensure you have the Xcode Command Line Tools installed. This is the basic package for command line development on MacOS. Enter the following command and then follow the prompts.

   $ xcode-select --install

   **Note** that there is a different version of the Xcode tools for each major release of MacOS. If you have upgraded your OS, say from Ventura to Sonoma, you must install the corresponding version of Xcode.

   If the installer says Xcode is installed but an upgrade exists, you almost certainly want to apply that upgrade. 

#. Make a folder for |NUMBAT| studies and clone the github repository::

   $ mkdir numbat
   $ cd numbat
   $ git clone https://github.com/michaeljsteel/NumBAT.git
   $ cd NumBAT

#. If it is not already on your system, install the MacPorts package manager at `this page <https://macports.org/install.php>`_.

#. Install the Gmsh mesh generation tool at `this page <https://gmsh.info>`_.
   Just the main Gmsh installer is fine. The SDK and other features are not required.

   **Note:** After the installer has run, you must move the Gmsh application into your Applications 
   folder by dragging the Gmsh icon into Applications.


#. Install a current gcc (we used gcc)::

   $ sudo port install gcc13

#. Install the Lapack and Blas linear algebra libraries::

   $ sudo port install lapack

#. Install the Arpack eigensolver::

   $ sudo port install arpack

#. Install the SuiteSparse matrix algebra suite::

   $ sudo port install suitesparse

#. Install a current python (we used python 3.12):

   Use the standard installer at `<https://www.python.org/downloads/macos/>`_.

   (Note that this will install everything in `/Library/Frameworks` and **not** override
   the System python in `/System/Library/Frameworks.`)

#. Install python `virtualenv` package

   $ cd /Library/Frameworks/Python.framework/Versions/3.12/bin/
   $ ./python3.12 -m pip install --upgrade pip
   $ ./pip3 install virtualenv

#. Create a |NUMBAT|-specific python virtual environment in `~/nbpy3`

    $ cd /Library/Frameworks/Python.framework/Versions/3.12/bin/
    $ ./python3 -m virtualenv ~/nbpy3


#. Activate the new python virtual environment (note the leading fullstop) ::

    $ . ~/nbpy3/bin/activate

#. Install necessary python libraries ::

    $ pip3 install numpy matplotlib scipy psutil

#. Check that the python installs work and create a matplotlib .config directory::

   $ python3.12
   $ import matplotlib
   $ import numpy
   $ import scipy
   $ import psutil

#. Install the |NUMBAT| matplotlib style file::

   $ mkdir -pR $HOME/.matplotlib/stylelib/
   $ cp <NumBAT>/backend/|NUMBAT|style.mplstyle $HOME/.matplotlib/stylelib

#. Move to the |NUMBAT| fortran directory::

   $ cd backend/fortran

#. Open the file `Makefile` in your preferred text editor and edit the lines at the top  of the file so that: 
   
     - The line `PLAT=MacOS` is active and the others are commented out with a leading `#` symbol. 
     - The value of `MYPYENV` matches the folder of your python virtual environment set up above.
     - The value of `PYVERMAJMIN` and `SYSTEMPYINC` are set appropriately.  

#. Now at last, we can build |NUMBAT| ::

   $ make

#. If this completes without error, you are ready to proceed to the next chapter to begin using |NUMBAT|.

#. If you hit a compile error you can't resolve, please get in touch at |NUMBAT_EMAIL|.


Installing on Windows
================================

The easiest way to run |NUMBAT| on  Windows is usually by installing Ubuntu as
a virtual machine using either `Microsoft Hyper-V
<https://wiki.ubuntu.com/Hyper-V>`_ or `Oracle Virtual Box
<https://ubuntu.com/tutorials/how-to-run-ubuntu-desktop-on-a-virtual-machine-using-virtualbox#1-overview>`_.

Then |NUMBAT| can be installed using exactly the same procedure as described
above for standard Linux installations.  It is also possible to build |NUMBAT|
using the `Windows Subsystem for Linux
<https://msdn.microsoft.com/en-au/commandline/wsl/install_guide>`_, but dealing
with installing the additional required packages may be quite painful.


There is also an outdated Docker package for Windows that could be updated to support the current version of |NUMBAT| if there is demand. Let us know.

.. On non-ubuntu OSes you may also need to compile a local version of Suitesparse, which is described in the next section.

.. Manual installation of SuiteSparse
.. ----------------------------------

.. The FEM routine used in |NUMBAT| makes use of the highly optimised `UMFPACK <https://www.cise.ufl.edu/research/sparse/umfpack/>`_ (Unsymmetric MultiFrontal Package) direct solver for sparse matrices developed by Prof. Timothy A. Davis. This is distributed as part of the  SuiteSparse libraries under a GPL license. It can be downloaded from `https://www.cise.ufl.edu/research/sparse/SuiteSparse/ <https://www.cise.ufl.edu/research/sparse/SuiteSparse/>`_

.. This is the process we have used in the past, however this was some years ago and may need to be modified.

.. Unpack SuiteSparse into ``|NUMBAT|/backend/fortran/``, it should create a directory there; ``SuiteSparse/``
.. Make a directory where you want SuiteSparse installed, in my case SS_installed ::

    $ mkdir SS_installed/

.. edit SuiteSparse/SuiteSparse\_config/SuiteSparse\_config.mk for consistency across the whole build; i.e. if using intel fortran compiler ::

    line 75 F77 = gfortran --> ifort

.. set path to install folder::

    line 85 INSTALL_LIB = /$Path_to_EMustack/|NUMBAT|/backend/fortran/SS_installed/lib
    line 86 INSTALL_INCLUDE = /$Path_to_EMustack/|NUMBAT|/backend/fortran/SS_installed/include

.. line 290ish commenting out all other references to these::

    F77 = ifort
    CC = icc
    BLAS   = -L/apps/intel-ct/12.1.9.293/mkl/lib/intel64 -lmkl_rt
    LAPACK = -L/apps/intel-ct/12.1.9.293/mkl/lib/intel64 -lmkl_rt

.. Now make new directories for the paths you gave 2 steps back::

    $ mkdir SS_installed/lib SS_installed/include

.. Download `metis-4.0 <http://glaros.dtc.umn.edu/gkhome/fsroot/sw/metis/OLD>`_ and unpack metis into SuiteSparse/ Now move to the metis directory::

    $ cd SuiteSparse/metis-4.0

.. Optionally edit ``metis-4.0/Makefile.in`` as per ``SuiteSparse/README.txt`` plus with ``-fPIC``::

    CC = gcc
    or
    CC = icc
    OPTFLAGS = -O3 -fPIC

.. Now make ``metis`` (still in SuiteSparse/metis-4.0/)::

    $ make

.. Now move back to ``|NUMBAT|/backend/fortran/`` ::

    $ cp SuiteSparse/metis-4.0/libmetis.a SS_installed/lib/

.. and then move to ``SuiteSparse/`` and execute the following::

    $ make library
    $ make install
    $ cd SuiteSparse/UMFPACK/Demo
    $ make fortran64
    $ cp SuiteSparse/UMFPACK/Demo/umf4_f77zwrapper64.o into SS_installed/lib/

.. Copy the libraries into ``|NUMBAT|/backend/fortran/Lib/`` so that ``|NUMBAT|/`` is a complete package that can be moved across machine without alteration. This will override the pre-compiled libraries from the release (you may wish to save these somewhere).::

    $ cp SS_installed/lib/*.a |NUMBAT|/backend/fortran/Lib/
    $ cp SS_installed/lib/umf4_f77zwrapper64.o |NUMBAT|/backend/fortran/Lib/


.. Edit ``|NUMBAT|/backend/fortran/Makefile`` to reflect what compiler you are using and how you installed the libraries. The Makefile has further details.

.. Then finally run the setup.sh script!

#. To build the pdf documentation you are currently reading, use ::

    $ make docs

Note however most of the figures will only be available after you have run all the example problems in the  ``tutorial`` ``lit_ex`` and ``JOSAB_tutorial`` directories.
This can be done by running ``make`` in each of those directories. Be aware that  some of these problems are quite large and may require some time to complete depending on your computer's performance.


Building the documentation
===============================

You can rebuild the documentation you are currently reading by moving into the ``<NumBAT>/docs`` directory and running either ``make html`` or ``make latexpdf``.

In each case, the output is placed in the ``<NumBAT>/docs/build`` directory.


Note however most of the figures will only be available after you have run all the example problems in the  ``tutorial`` ``lit_ex`` and ``JOSAB_tutorial`` directories.
This can be done by running ``make`` in each of those directories. Be aware that  some of these problems are quite large and may require some time to complete depending on your computer's performance.





