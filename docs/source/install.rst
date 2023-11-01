.. include:: numbatdefs.txt

*******************
Installing |NUMBAT| 
*******************

Introduction 
================================
While |NUMBAT| is developed on Linux, it can also be built on MacOS X as a 
native command-line application, and under Windows using a virtual machine running Linux.

In all cases, the current source code for NumBAT is hosted `here on Github <https://github.com/michaeljsteel/NumBAT>`_. Please always download the latest release from the github page.

The following sections provide instructions on installing for each platform.
Please let us know of any difficulties you encounter, or suggestions for improvements to the install procedure for MacOS or Windows, by sending an email to 
michael.steel@mq.edu.au.

Installing on Linux
================================


NumBAT has been developed and tested on Ubuntu 23.04 with the following package versions: Python 3.11.4, Numpy 1.24.2, Arpack-NG, Suitesparse 7.1.0, and Gmsh 4.8.4.  NumBAT also depends on the BLAS and LAPACK libraries. We strongly recommend linking NumBAT against optimised versions, such as the MKL library provided in the  free Intel OneAPI library.

|NUMBAT| has also been successfully installed by users on Debian and RedHat/Fedora, and with different versions of packages, but these installations have not been as thoroughly documented so may require user testing.  In general, any relatively current Linux system should work without trouble.


Before installing, you may wish to update your system ::


    $ sudo apt-get update
    $ sudo apt-get upgrade

There is no need to install NumBAT in a central location such as `/usr/local` though you may choose to do so.

To download the current version from the git repository and install any missing library dependencies, use
 ::

    $ git clone https://github.com/michaeljsteel/NumBAT.git
    $ cd NumBAT/
    $ sudo make installdeps 

In this documentation, we indicate the root NumBAT directory (e.g. ``/usr/local/NumBAT``, ``/home/mike/NumBAT``) with ``<NumBAT>``. 

Before attempting to build the NumBAT module, open the file ``<NumBAT>/backend/fortran/Makefile`` in a text editor and check the settings associated with the variables ``PLAT`` that control the preferred math library. The starting makefile is setup to look for the Intel OneAPI library which is the recommended configuration.

You may now build the NumBAT module by running the following in the ``<NumBAT>`` directory. ::

    $ make build

If this completes without error, you may run a short test suite with  ::

    $ make tests

To build the pdf documentation you are currently reading, use ::

    $ make docs 

Note however most of the figures will only be available after you have run all the example problems in the  ``tutorial`` ``lit_ex`` and ``JOSAB_tutorial`` directories.
This can be done by running ``make`` in each of those directories. Be aware that  some of these problems are quite large and may require some time to complete depending on your computer's performance.


Other build configurations 
--------------------------

The Fortran components (NumBAT source code and libraries) have been successfully compiled with Intel's ``ifortran`` as well as GCC's open-source ``gfortran``. In this documentation we use ``gfortran``, but this can be easily adjusted in ``NumBAT/backend/fortran/Makefile``


Installing on MacOS
================================
|NUMBAT| can also be installed on MacOS, though this is currently somewhat experimental and has only been performed on certain versions of MacOS.  Any comments on difficulties and solutions will be appreciated.

The following steps have worked for us:

#. Open a terminal window on your desktop.

#. Make a folder for |NUMBAT| studies and clone the github repository::

   $ mkdir numbat
   $ cd numbat 
   $ git clone https://github.com/michaeljsteel/NumBAT.git
   $ cd numbat

#. If it is not already on your system, install the MacPorts package manager at `this page <https://macports.org/install.php>`_.
   
#. Install the Gmsh mesh generation tool at `this page <https://gmsh.info>`_.
   Just the main Gmsh installer is fine. The SDK and other features are not required.

   Install the Gmsh application into your Applications folder by dragging the Gmsh icon into Applications.


#. Install a current gcc (we used gcc)::

   $ sudo port install gcc13

#. Install the Lapack and Blas linear algebra libraries::

   $ sudo port install lapack

   $ sudo port install blas

#. Install the Arpack eigensolver::

   $ sudo port install arpack

#. Install the SuiteSparse matrix algebra suite::

   $ sudo port install suitesparse

#. Install a current python (we used python 3.12):

   Use the standard installer at `<https://www.python.org/downloads/macos/>`_.

   (Note that this will install everything in `/Library/Frameworks` and **not** override 
   the System python in `/System/Library/Frameworks.`)

#. Install some python packages not included by default::
   
   $ pip3 install numpy scipy matplotlib psutil

#. Check that the python installs work and create a matplotlib .config directory::

   $ python3.12
   $ import matplotlib
   $ import numpy
   $ import scipy
   $ import psutil

#. Install the |NUMBAT| matplotlib style file::

   $ mkdir -pR $HOME/.matplotlib/stylelib/
   $ cp <NumBAT>/backend/NumBATstyle.mplstyle $HOME/.matplotlib/stylelib

#. Move to the |NUMBAT| fortran directory::

   $ cd backend/fortran

#. Open the Makefile in a text editor and edit the lines at the top  of the file so that the line `PLAT=MacOS` is active and the others are commented out. 

#. Now we can build |NUMBAT| ::
   
   $ make - C Makefile

Installing on Windows
================================

The easiest way to run |NUMBAT| on  Windows is usually by installing Ubuntu as a virtual machine using either `Microsoft Hyper-V <https://wiki.ubuntu.com/Hyper-V>`_ or `Oracle Virtual Box <https://ubuntu.com/tutorials/how-to-run-ubuntu-desktop-on-a-virtual-machine-using-virtualbox#1-overview>`_.

Then |NUMBAT| can be installed using exactly the same procedure as described above for standard Linux installations.
It is also possible to build |NUMBAT| using the `Windows Subsystem for Linux <https://msdn.microsoft.com/en-au/commandline/wsl/install_guide>`_, but dealing with installing the additional required packages may be quite painful.



.. On non-ubuntu OSes you may also need to compile a local version of Suitesparse, which is described in the next section.

.. Manual installation of SuiteSparse
.. ----------------------------------

.. The FEM routine used in NumBAT makes use of the highly optimised `UMFPACK <https://www.cise.ufl.edu/research/sparse/umfpack/>`_ (Unsymmetric MultiFrontal Package) direct solver for sparse matrices developed by Prof. Timothy A. Davis. This is distributed as part of the  SuiteSparse libraries under a GPL license. It can be downloaded from `https://www.cise.ufl.edu/research/sparse/SuiteSparse/ <https://www.cise.ufl.edu/research/sparse/SuiteSparse/>`_

.. This is the process we have used in the past, however this was some years ago and may need to be modified.

.. Unpack SuiteSparse into ``NumBAT/backend/fortran/``, it should create a directory there; ``SuiteSparse/``
.. Make a directory where you want SuiteSparse installed, in my case SS_installed ::

    $ mkdir SS_installed/

.. edit SuiteSparse/SuiteSparse\_config/SuiteSparse\_config.mk for consistency across the whole build; i.e. if using intel fortran compiler ::

    line 75 F77 = gfortran --> ifort

.. set path to install folder::

    line 85 INSTALL_LIB = /$Path_to_EMustack/NumBAT/backend/fortran/SS_installed/lib
    line 86 INSTALL_INCLUDE = /$Path_to_EMustack/NumBAT/backend/fortran/SS_installed/include

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

.. Now move back to ``NumBAT/backend/fortran/`` ::

    $ cp SuiteSparse/metis-4.0/libmetis.a SS_installed/lib/

.. and then move to ``SuiteSparse/`` and execute the following::

    $ make library
    $ make install
    $ cd SuiteSparse/UMFPACK/Demo
    $ make fortran64
    $ cp SuiteSparse/UMFPACK/Demo/umf4_f77zwrapper64.o into SS_installed/lib/

.. Copy the libraries into ``NumBAT/backend/fortran/Lib/`` so that ``NumBAT/`` is a complete package that can be moved across machine without alteration. This will override the pre-compiled libraries from the release (you may wish to save these somewhere).::

    $ cp SS_installed/lib/*.a NumBAT/backend/fortran/Lib/
    $ cp SS_installed/lib/umf4_f77zwrapper64.o NumBAT/backend/fortran/Lib/


.. NumBAT Makefile

.. Edit ``NumBAT/backend/fortran/Makefile`` to reflect what compiler you are using and how you installed the libraries. The Makefile has further details.

.. Then finally run the setup.sh script!

.. _sec-contribute-label:

Contributing to NumBAT
----------------------------------

NumBAT is open source software licensed under the GPL with all source and documentation available
at `github.com <https://github.com/michaeljsteel/NumBAT.git>`_. We welcome additions to NumBAT code, documentation and the materials library. Interested users should fork the standard release from github and make a pull request when ready.  For major changes, we strongly suggest contacting the NumBAT team before starting work at ``michael.steel@mq.edu.au``.
