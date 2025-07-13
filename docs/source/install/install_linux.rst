.. include:: ../numbatdefs.txt

.. _sec-linuxinstall-label:


Installing on Linux
================================


Install locations
---------------------

There is no need to install |NUMBAT| in a central location
such as ``/usr/local/`` or ``/opt/local/`` though you may certainly choose to do so.


Here and throughout this documentation, we use the string ``<NumBAT>`` to indicate the root |NUMBAT| install directory (e.g. ``/usr/local/NumBAT``, ``/home/mike/NumBAT``, ``/home/myuserid/research/NumBAT``).


Requirements
---------------------
|NUMBAT| is developed and tested using relatively recent operating system, compiler and libraries. You should not need the very latest releases, but in general compilation will be smoother on an up-to-date system. In particular, we recommend using:

  - an OS release from 2023 or later (eg Ubuntu 24.04/24.10)
  - gcc compiler of version 13.0 or later
  - python 3.11 or later

|NUMBAT| is currently developed and tested on Ubuntu 25.04 with the following
package versions: Python 3.13, Numpy 2.0, Arpack-NG, Suitesparse 7.1.0,
and Gmsh 4.8.4.  |NUMBAT| also depends on the BLAS linear algebra library. We
strongly recommend linking |NUMBAT| against an optimised version, such as the MKL
library provided in the free Intel OneAPI library (for Intel CPUs) or the AMD Optimizing CPU Libraries (AOCL) for AMD CPUs.  The steps below demonstrate the Intel OneAPI approach.


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

Required libraries
--------------------------

#. Before installing, ensure your system is up to date ::

    $ sudo apt-get update
    $ sudo apt-get upgrade

#. Install the required libraries using your distribution's package manager.

   On Ubuntu, perform the following ::

    $ sudo apt-get install gcc gfortran make gmsh python3-venv python3-dev meson pkg-config ninja-build

    $ sudo add-apt-repository universe

    $ sudo apt-get install libarpack2-dev libparpack2-dev libatlas-base-dev libblas-dev liblapack-dev  libsuitesparse-dev

#. If you wish to use the Intel OneAPI math libraries, you need both of the following:

    - *Intel OneAPI Base Toolkit*:
        This is the main Intel developer environment including C/C++ compiler and many high performance math libraries.

        Download and run the `installer <https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html>`_ accepting all defaults.


    - *Intel OneAPI HPC Toolkit*
        This adds the Intel Fortran compiler amongst other HPC tools.

        Download and run the `installer <https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html>`__ accepting all defaults.


#. If you using the Intel OneAPI math libraries, you should add the library path
   ``/opt/intel/oneapi/<release>/lib`` to your ``LD_LIBRARY_PATH`` variable in one of your shell startup files (eg. ``~/.bashrc``).  Replace ``<release>`` with the correct string ``2024.1`` or similar depending on your installed version of OneAPI.



Building |NUMBAT| itself
--------------------------

#. Create a python virtual environment for working with |NUMBAT|.
   You can use any name and location for your tree.
   To specify a virtual environment tree called `nbpy3` in your home directory, enter ::

    $ cd ~
    $ python3 -m venv nbpy3


#. Activate the new python virtual environment ::

      $ source ~/nbpy3/bin/activate

#. Install necessary python libraries ::

      $ pip3 install numpy matplotlib scipy psutils



#. Create a working directory for your |NUMBAT| work and move into it. From now, we will refer to this location as ``<NumBAT>``.

#. To download the current version from the git repository and install any missing library dependencies, use ::

    $ git clone https://github.com/michaeljsteel/NumBAT.git
    $ cd NumBAT

#.  Move to the ``backend\fortran`` directory.

    #. To build with the ``gcc`` compilers, run::

        $ make gcc

    #. To build with the Intel compilers, edit the file ``nb-linuxintel-native-file.ini`` adjusting the variables to point the correct location of the Intel compilers. Then run::

        $ make intel


#. If all is well, this will run to completion. If you encounter errors, please check that all the instructions above have been followed accurately. If you are still stuck, see :ref:`sec-troubleshooting-linux-label` for further ideas.

#. If you hit a compile error you can't resolve, please see the instructions at :ref:`sec-helpinstall-label` on how to seek help.

#. Once the build has apparently succeeded, it is time to test the installation with a short script that tests whether required applications and libraries can be found and loaded. Perform the following commands::

    $ cd <NumBAT>/backend
    $ python ./nb_install_tester.py

#. If this program runs without error, congratulations! You are now ready  to proceed to the next chapter to begin using |NUMBAT|

#. Once again, if you run into trouble, please don't hesitate to get in touch for help using the instructions at :ref:`sec-helpinstall-label`. Please do send all the requested information, as it usually allows us to solve your problem faster.



Using the Intel Fortran compiler
----------------------------------

The default compiler for Linux is GCC's ``gfortran``.

It is also possible to build |NUMBAT| with the ``ifx`` compiler from Intel's free OneAPI HPC toolkit. You may find some modest performance improvements.

To use the Intel compiler,

#. If you have not already done so, install the Intel OneAPI Base and HPC Toolkits as described above.
#. Adjust your LD_LIBRARY_PATH variable in your ``~/.bashrc`` or equivalent to include ``/opt/intel/oneapi/<release>/lib``.  (Replace ``<release>`` with the correct string ``2024.1`` or similar depending on your installed version of OneAPI.)
#. In ``<NumBAT>/backend/fortran``, repeat all the earlier instructions for the standard GCC build but rather than plain ``make gcc``, please use::

    $ make intel


Installing without root access
----------------------------------------------------
Compiling and installing |NUMBAT| itself does not rely on having root access to your machine. However, installing the supporting libraries such as SuiteSparse and Arpack
is certainly simpler if you have root or the assistance of your system admin.

If this is not possible, you can certainly proceed by building and installing all the required libraries into your own tree within your home directory.
It may be helpful to create a tree like the following so that the relevant paths mirror those of a system install ::

   $HOME/
    |---my_sys/
         |---usr/
              |---include/
              |---lib/
              |---bin/




.. _sec-troubleshooting-linux-label:

Troubleshooting Linux installs
-------------------------------------------
Performing a full build of |NUMBAT| and all its libraries from scratch is a non-trivial task and it's possible you will hit a few stumbles.
Here are a few traps to watch out for:

#. Please ensure to use relatively recent libraries for all the Python components. This includes using

   - Python: 3.11 or later
   - ``matplotlib``: 3.9.0 or later
   - ``scipy``:  1.13.0 or later
   - ``numpy``:  2.0 or later


#. Be sure to follow the instructions above about setting up the virtual environment for |NUMBAT| excusively. This will help prevent incompatible Python modules being added over time.


#. If you encounter an error about "missing symbols" in the NumBAT fortran module, there are usually two possibilities:

   - A shared library (a file ending in ``.so``) is not being loaded correctly because it can't be found in the standard search path. To detect this, run ``ldd nb_fortran.so`` in the ``backend/fortran`` directory and look for any lines containing ``not found``. (On MacOS, use ``otools -L nb_fortran.so`` instead of ``ldd``.)


     You may need to add the directory containing the relevant libraries to your ``LD_LIBRARY_PATH`` in your shell setup files (eg. ``~/.bashrc`` or equivalent).

   - You may have actually encountered a bug in the |NUMBAT| build process. Contact us for assistance as described in the introduction.

#. If |NUMBAT| crashes during execution with a ``Segmentation fault``, you have quite possibly found a bug that should be reported to us. Useful information about a crash can be obtained from the GNU debugger ``gdb`` as follows:

    #. Make sure that core dumps are enabled on your system. This `article <https://medium.com/@sourabhedake/core-dumps-how-to-enable-them-73856a437711>`_ provides an excellent guide on how to do so.

    #. Ensure that ``gdb`` is installed.

    #. Rerun the script that causes the crash. You should now have a core dump in the directory determined in step 1.

    #. Execute ``gdb`` as follows::

        $ gdb <path_to_numbat_python_env>  <path_to_core file>

    #. In gdb, enter ``bt`` for *backtrace* and try to identify the point in the code at which the crash has occurred.


