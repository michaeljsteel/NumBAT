
.. include:: numbatdefs.txt

.. _chap-install-label:

*******************
Installing NumBAT
*******************


This chapter provides instructions on installing |NUMBAT| on each platform.
Please email |NUMBAT_EMAIL| to let us know of any difficulties you encounter, or suggestions for improvements to the install procedure on any platform, but especially for MacOS or Windows.

Information for all platforms
================================
While |NUMBAT| is developed on Linux, it can also be built natively on both MacOS and Windows.
The Linux builds can also be run under virtual machines on MacOS and Windows if desired.

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
and Gmsh 4.8.4.  |NUMBAT| also depends on the BLAS and :eq: libraries. We
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

Required libraries
--------------------------

#. Before installing, ensure your system is up to date ::

    $ sudo apt-get update
    $ sudo apt-get upgrade

#. Install the required libraries using your distribution's package manager.

   On Ubuntu, perform the following::

    $ sudo apt-get install gcc gfortran make gmsh python3-venv python3-dev

    $ sudo apt-get install meson pkg-config ninja-build

    $ sudo add-apt-repository universe

    $ sudo apt-get install libarpack2-dev libparpack2-dev

    $ sudo apt-get install libatlas-base-dev libblas-dev liblapack-dev

    $ sudo apt-get install libsuitesparse-dev

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


   Ensure that your ``numpy`` version  is from the 1.26.x and not the new 2.0.0 line.

#. If you wish to be able to rebuild the documentation, we need some additional modules ::
     $ pip3 install sphinx sphinxcontrib-bibtex setuptools

#. Create a working directory for your |NUMBAT| work and move into it. From now, we will refer to this location as ``<NumBAT>``.

#. To download the current version from the git repository and install any missing library dependencies, use ::

    $ git clone https://github.com/michaeljsteel/NumBAT.git
    $ cd NumBAT

#.  Move to the ``backend\fortran`` directory.

    #. To build with the ``gcc`` compilers, run::

        $ make gcc

    #. To build with the Intel compilers, edit the file ``nb-linuxintel-native-file.ini`` adjusting the variables to point the correct location of the Intel compilers. Then run::

        $ make intel


#. If all is well, this will run to completion. If you encounter errors, please check that all the instructions above have been followed accurately. If you are still stuck, see :ref:`sec-troubleshooting-label` for further ideas.

#. If you hit a compile error you can't resolve, please see the instructions at :ref:`sec-helpinstall-label` on how to seek help.

#. Once the build has apparently succeeded, it is time to test the installation with a short script that tests whether required applications and libraries can be found and loaded. Perform the following commands::

    $ cd <NumBAT>/backend
    $ python ./nb_install_tester.py

#. If this program runs without error, congratulations! You are now ready  to proceed to the next chapter to begin using |NUMBAT|




Other build configurations
--------------------------

The default compiler for Linux is GCC's ``gfortran``.

It is also possible to build |NUMBAT| with the ``ifx`` compiler from Intel's free OneAPI HPC toolkit.

To do so,

#. Install the Intel OneAPI Base and HPC Toolkits.
#. Adjust your LD_LIBRARY_PATH variable in your ``~/.bashrc`` or equivalent to include ``/opt/intel/oneapi/<release>/lib``.  (Replace ``<release>`` with the correct string ``2024.1`` or similar depending on your installed version of OneAPI.)
#. In ``<NumBAT>/backend/fortran``, repeat all the earlier instructions for the standard GCC build but rather than plain ``make gcc``, please use::

    $ make intel


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


.. _sec-troubleshooting-label:


Troubleshooting Linux and MacOS installs
-------------------------------------------
Performing a full build of |NUMBAT| and all its libraries from scratch is a non-trivial task and it's possible you will hit a few stumbles.
Here are a few traps to watch out for:

#. Please ensure to use relatively recent libraries for all the Python components. This includes using

   - Python: 3.10 or later
   - ``matplotlib``: 3.9.0 or later
   - ``scipy``:  1.13.0 or later
   - ``numpy``:  1.26.2 or later

#. But try not to use very recently released major upgrades.

   Notably the 2.0.0 series of ``numpy``, which was only released in mid-June 2024 includes major changes to ``numpy`` architecture and is  not yet supported.

#. Be sure to follow the instructions above about setting up the virtual environment for |NUMBAT| excusively. This will help prevent incompatible Python modules being added over time.

#. In general, the GCC build is more tested and forgiving than the build with the Intel compilers and we recommend the GCC option.  However, we do recommend using the Intel OneAPI math libraries as described above. This is the easiest way to get very high performance LAPACK and BLAS libraries with a well-designed directory tree.

#. If you encounter an error about "missing symbols" in the NumBAT fortran module, there are usually two possibilities:

   - A shared library (a file ending in ``.so``) is not being loaded correctly because it can't be found in the standard search path. To detect this, run ``ldd nb_fortran.so`` in the ``backend/fortran`` directory and look for any lines containing ``not found``.


     You may need to add the directory containing the relevant libraries to your ``LD_LIBRARY_PATH`` in your shell setup files (eg. ``~/.bashrc`` or equivalent).

   - You may have actually encountered a bug in the |NUMBAT| build process. Contact us for assistance as described in the introduction.

#. If |NUMBAT| crashes during execution with a ``Segmentation fault``, useful information can be obtained from the GNU debugger ``gdb`` as follows:

    #. Make sure that core dumps are enabled on your system. This `article <https://medium.com/@sourabhedake/core-dumps-how-to-enable-them-73856a437711>`_ provides an excellent guide on how to do so.

    #. Ensure that ``gdb`` is installed.

    #. Rerun the script that causes the crash. You should now have a core dump in the directory determined in step 1.

    #. Execute ``gdb`` as follows::

        $ gdb <path_to_numbat_python_env>  <path_to_core file>

    #. In gdb, enter ``bt`` for *backtrace* and try to identify the point in the code at which the crash has occurred.

Installing on MacOS
================================

|NUMBAT| can also be installed on MacOS, though this is currently somewhat experimental and has only been performed on certain versions of MacOS.  Any comments on difficulties and solutions will be appreciated.

The following steps have worked for us:

#. Open a terminal window on your desktop.

#. Ensure you have the Xcode Command Line Tools installed. This is the basic package for command line development on MacOS. If you do not or are not sure, enter the following command and then follow the prompts::

   $ xcode-select --install

   **Note** that there is a different version of the Xcode tools for each major release of MacOS. If you have upgraded your OS, say from Ventura to Sonoma, you must install the corresponding version of Xcode.

   If the installer says Xcode is installed but an upgrade exists, you almost certainly want to apply that upgrade.

#. Make a folder for your |NUMBAT| work in a suitable location in your folder tree. Then clone the github repository::

   $ mkdir numbat
   $ cd numbat
   $ git clone https://github.com/michaeljsteel/NumBAT.git
   $ cd NumBAT

This new ``NumBAT`` folder location is referred to as ``<NumBAT>`` in the following.


#. If it is not already on your system, install the  `MacPorts package manager <https://macports.org/install.php>`_ using the appropriate installer for your version of MacOS at that page.

#. Install the  `Gmsh  <https://gmsh.info>`_ mesh generation tool.
   Just the main Gmsh installer is fine. The SDK and other features are not required.

   **Note:** After the installer has run, you must move the Gmsh application into your Applications
   folder by dragging the Gmsh icon into Applications.


#. Install a current gcc (we used gcc13)::

   $ sudo /opt/local/bin/port install gcc-devel

#. Install the Lapack and Blas linear algebra libraries::

   $ sudo /opt/local/bin/port install lapack

#. Install the Arpack eigensolver::

   $ sudo /opt/local/bin/port install arpack

#. Install the SuiteSparse matrix algebra suite::

   $ sudo /opt/local/bin/port install suitesparse

#. Install a current python (we used python 3.12):

   Use the standard installer at `<https://www.python.org/downloads/macos/>`_.

   (Note that this will install everything in `/Library/Frameworks` and **not** override
   the System python in `/System/Library/Frameworks.`)

#. Install python `virtualenv` package ::

   $ cd /Library/Frameworks/Python.framework/Versions/3.12/bin/
   $ ./python3.12 -m pip install --upgrade pip
   $ ./pip3 install virtualenv

#. Create a |NUMBAT|-specific python virtual environment in `~/nbpy3` ::

    $ cd /Library/Frameworks/Python.framework/Versions/3.12/bin/
    $ ./python3 -m virtualenv ~/nbpy3


#. Activate the new python virtual environment (note the leading fullstop) ::

    $ . ~/nbpy3/bin/activate

#. Install necessary python libraries ::

    $ pip3 install numpy matplotlib scipy psutil meson ninja

#. Check that the python installs work.

   $ python3.12
   >>> import matplotlib
   >>> import numpy
   >>> import scipy
   >>> import psutil

#. Install the |NUMBAT| matplotlib style file::

   $ mkdir -p $HOME/.matplotlib/stylelib/
   $ cp <NumBAT>/backend/|NUMBAT|style.mplstyle $HOME/.matplotlib/stylelib

#. Move to the |NUMBAT| fortran directory::

   $ cd backend/fortran

#. Move to the ``<NumBAT>/backend/fortran/`` directory and open the file ``meson.options`` in a text editor. Check the values of the options in the ``MacOS`` section and change any of the paths in the ``value`` fields as required.

#. To start the build, enter::

   $ make mac

#. If all is well, this will run to completion. If you encounter errors, please check that all the instructions above have been followed accurately. If you are still stuck, see :ref:`sec-troubleshooting-label` for further ideas.

#. If you hit a compile error you can't resolve, please see the instructions at :ref:`sec-helpinstall-label` on how to seek help.

#. Once the build has apparently succeeded, it is time to test the installation with a short script that tests whether required applications and libraries can be found and loaded. Perform the following commands::

    $ cd <NumBAT>/backend
    $ python ./nb_install_tester.py

#. If this program runs without error, congratulations! You are now ready  to proceed to the next chapter to begin using |NUMBAT|


Installing on Windows
================================

There are several approaches for installing |NUMBAT| on Windows. You can build the entire system from scratch, as with the Linux and MacOS installs, or you can use a pre-built installer available from the github page.

In both cases, you need to set up the python environment.  This is described first below, followed by the two alternative
methods of installing |NUMBAT| itself.

Setting up Python on Windows
-------------------------------------

Python installers can be downloaded from the `Python website <https://www.python.org/downloads/windows/>`_.


  #. If you do not have a current Python, download and run the installer for a recent version
  from the `Python website <https://www.python.org/downloads/windows/>`_.

  By default, this will install Python in the directory ``%HOMEPATH%\AppData\Local\Programs\Python\<PythonVer>``,
  say ``%HOMEPATH%\AppData\Local\Programs\Python\Python312``.


  #. Create a python virtual environment for working with |NUMBAT|.
      You can use any name and location for your environment.

    To specify a virtual environment tree called `nbpy3`, open a command prompt (or your favourite Windows terminal app)
    from the Start Menu and  enter ::

        $ %HOMEPATH%\AppData\Local\Programs\Python\Python312\python.exe -m venv nbpy3


  #. Activate the new python virtual environment ::

       $ %HOMEPATH%\npby3\Scripts\activate

  #. Install the necessary python tools and libraries ::

     $ pip3 install numpy==1.26.4 matplotlib scipy psutil ninja meson==1.4.1

   Note that at last check, the most recent meson (1.5.0) is broken and we specify the earlier 1.4.1 version.

   Similarly we specify a version of ``numpy`` from the 1.26 series as the new 2.0 version is not yet
   supported by some other packages we use.

  #. Finally, we will also need the Gnu ``make`` tool. This can be installed by typing

     $ winget install ezwinports.make

    and then starting a new terminal (so that the PATH variable is updated to find ``make.exe``.)

  #. Now you can proceed to install |NUMBAT| using either Method 1,  building fully from source, or Method 2, using the pre-built installer from the github page.


Installing |NUMBAT| Method 1: full build from source
----------------------------------------------------------


The Windows version of |NUMBAT| is built using the native Windows toolchain including Visual Studio and the Intel Fortran compiler. There is a substantial number of steps and tools required, but it should go relatively smoothly.


Windows build tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following tools are all free but will use several GB of disk space.

    - *Visual Studio*:
        This is the primary Microsoft development environment.

        To install the free Community 2022 edition, download the `main installer <https://visualstudio.microsoft.com/vs/community/>`_ and follow the instructions.

    - *Intel OneAPI Base Toolkit*:
        This is the main Intel developer environment including C/C++ compiler and many high performance math libraries.

        Download and run the `installer <https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html>`_ accepting all defaults.


    - *Intel OneAPI HPC Toolkit*
        This adds the Intel Fortran compiler amongst other HPC tools.

        Download and run the `installer <https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html>`__ accepting all defaults.


    - *Git*
        This is a source control that we use to download |NUMBAT| and some other tools.


        Download and run the  `latest Git for Windows release <https://git-scm.com/download/win>`_, accepting all defaults.

        Some users may prefer to use a graphical interface such as `GitHub Desktop <https://desktop.github.com/>`_. This is fine too.

    - *CMake*
        This is a cross-platform build tool we will need for building some of the libraries.

        Download and run the  `latest release <https://cmake.org/download/>`_ accepting all defaults.


NumBAT code and libraries
^^^^^^^^^^^^^^^^^^^^^^^^^
We can now build the supporting libraries, and then |NUMBAT| itself.

#. Choose a location for the base directory for building |NUMBAT| and supporting libraries,
    say ``c:\Users\<myname>\numbat``, which we will refer to as ``<NumBAT_BASE>``.

#. Use the Start Menu to open the *Intel OneAPI Command Prompt for Intel 64 for Visual Studio 2022*.
   This is simply a Windows terminal with some Intel compiler environment variables pre-defined.

#. In the terminal window, change to the ``<NumBAT_BASE>`` directory, then execute the following commands::

    $ mkdir nb_releases
    $ mkdir usr_local
    $ mkdir usr_local\include
    $ mkdir usr_local\lib
    $ mkdir usr_local\packages
    $ cd usr_local\packages
    $ git clone https://github.com/opencollab/arpack-ng.git arpack-ng
    $ git clone https://github.com/jlblancoc/suitesparse-metis-for-windows.git suitesparse-metis
    $ cd ..\..\nb_releases
    $ git clone https://github.com/michaeljsteel/NumBAT.git nb_latest

#. Download the `Windows build of Gmsh <https://gmsh.info>`_ and unzip the tree into ``usr_local\packages\gmsh``.  The Gmsh executable should now be at ``<NumBAT>\usr_local\packages\gmsh\gmsh.exe``.



#. Your ``<NumBAT_BASE>`` tree should now look like this::

    %HOME%
        |---numbat
            |---nb_releases
            |---usr_local
                |---include
                |---lib
                |---packages



Building SuiteSparse
""""""""""""""""""""""
This library performs sparse matrix algebra, used in the eigensolving routines of |NUMBAT|.

1. In the Intel command terminal, cd to ``<NumBAT_BASE>\usr_local\packages\suitesparse-metis``.

2. Enter the following command. It may take a minute or two to complete::

     $ cmake -D WITH_MKL=ON -B build .

   Make sure to include the fullstop after ``build``.

3. If that completes correctly, use Windows Explorer to open ``<NumBAT_BASE>\usr_local\packages\suitesparse-metis\build\SuiteSparseProject.sln`` with Visual Studio 2022.

4. In the pull-down menu in the ribbon, select the *Release* build. Then from the *Build* menu, select the *Build Solution* item to commence the build.  This will take a couple of minutes.

5. Return to the command terminal and  cd to ``<NumBAT_BASE>\usr_local``. Then execute the following commands::

    $ copy packages\suitesparse-metis\build\lib\Release\*.dll lib
    $ copy packages\suitesparse-metis\build\lib\Release\*.lib lib
    $ copy packages\suitesparse-metis\SuiteSparse\AMD\Include\*.h include
    $ copy packages\suitesparse-metis\SuiteSparse\UMFPACK\Include\*.h include
    $ copy packages\suitesparse-metis\SuiteSparse\SuiteSparse_config\*.h include


Building Arpack-ng
""""""""""""""""""""""
This library performs an iterative algorithm for finding matrix eigensolutions.

1. In the Intel command terminal, cd to ``<NumBAT_BASE>\usr_local\packages\arpack-ng``.

2. Enter the following command. It may take a minute or two to complete::

      $ cmake -B build -T "fortran=ifx" -D CMAKE_BUILD_TYPE=Release -D BUILD_SHARED_LIBS=OFF .

   Note the final fullstop!

3. If that completes correctly, use Windows Explorer to open ``<NumBAT_BASE>\usr_local\packages\arpack-ng\build\arpack.sln`` with Visual Studio 2022.

4. In the pull-down menu in the ribbon, select the *Release* build. Then from the *Build* menu select
the *Build solution* option. This will take a few minutes.

5. Return to the command terminal and  cd to ``<NumBAT_BASE>\usr_local``. Then execute the following commands::

    $ copy packages\arpack-ng\build\Release\* lib
    $ copy packages\arpack-ng\ICB\*.h include


Building |NUMBAT|
""""""""""""""""""""""
At long last, we are ready to build |NUMBAT| itself.



#. Move to your root ``<NumBAT_BASE>`` directory and then to the |NUMBAT| folder itself::

    $ cd <NumBAT_BASE>
    $ cd nb_releases\nb_latest

   From this point, we refer to the current directory as ``<NumBAT>``.  In other words, ``<NumBAT> = <NumBAT_BASE>\nb_releases\nb_latest``.

#. Setup the environment variables for the Intel compiler::

    $  "c:\Program Files (x86)\Intel\oneAPI\setvars.bat"


#. Move to the ``<NumBAT>\backend\fortran`` directory and open the file ``meson.options`` in a text editor. Check the values of the options in the ``Windows`` section, particularly the value for ``windows_dir_nb_usrlocal`` and change any of the paths in the ``value`` fields as required.

#. To initiate the build, enter ::

     $ make win

   This should take 2 to 3 minutes.

.. #. Move to the ``<NumBAT>\backend\fortran\`` directory and open the file ``<NumBAT>\backend\fortran\Makefile.win`` in a text editor and change any absolute paths that involve your username. Now at last, we can build |NUMBAT| by running the following in the root ``<NumBAT>`` directory. ::

..    $ cd backend\fortran
..    $ make -f Makefile.win



#. If all is well, this will run to completion. If you encounter errors, please check that all the instructions above have been followed accurately. If you are still stuck, see :ref:`sec-troubleshooting-label` for further ideas.

#. If you hit a compile error you can't resolve, please see the instructions at :ref:`sec-helpinstall-label` on how to seek help.

#. Copy the .dlls built earlier to this directory::

    $ copy ..\..\..\..\usr_local\lib\*.dll .

#. Copy Intel oneAPI .dlls to this directory::

    $ copy "c:\Program Files (x86)\Intel\oneAPI\mkl\latest\bin\mkl_rt.2.dll" .
    $ copy "c:\Program Files (x86)\Intel\oneAPI\mkl\latest\bin\mkl_intel_thread.2.dll" .
    $ copy "c:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin\svml_dispmd.dll" .
    $ copy "c:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin\libmmd.dll" .
    $ copy "c:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin\libifcoremd.dll" .
    $ copy "c:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin\libifportMD.dll" .

#. At this point, we are ready to test the installation with a short script that checks whether required applications and libraries can be found and loaded. Perform the following commands::

    $ cd <NumBAT>/backend
    $ python .\nb_install_tester.py

#. If this program runs without error, congratulations! You are now ready  to proceed to the next chapter to begin using |NUMBAT|.  If not, please see the suggestions at :ref:`sec-troubleshooting-windows-label`.


Installing |NUMBAT| Method 2: pre-built Windows installer
-------------------------------------------------------------
*Note*: The installer will allow you to run |NUMBAT| problems and make changes to the code on the Python side only. (For most people, this is fine.)
The installer is not updated as frequently as the main source tree.


#. Choose a location for the base directory for building |NUMBAT| and supporting libraries,
    say ``c:\Users\<myname>\numbat``, which we will refer to as ``<NumBAT_BASE>``.

#. Use the Start Menu to open a Windows terminal.

#. In the terminal window, change to the ``<NumBAT_BASE>`` directory, then execute the following commands::

    $ mkdir nb_releases
    $ mkdir usr_local
    $ mkdir usr_local\packages

#. Download the `Windows build of Gmsh <https://gmsh.info>`_ and unzip the tree into ``usr_local\packages\gmsh``.  The Gmsh executable should now be at ``<NumBAT>\usr_local\packages\gmsh\gmsh.exe``.

#. Download the `Windows installer for |NUMBAT| <https://github.com/michaeljsteel/NumBAT>`_ from the |NUMBAT| github page. The link to the installer can be found at the bottom of the *Readme* section and also under the *Releases* heading in the right-hand column of the page.

   Run the installer choosing an install directory in the ``<NumBAT_BASE>\nb_releases`` folder.

#. In your Windows terminal, change directory to the newly installed |NUMBAT| folder.

#. Activate your python environment and then move to the |NUMBAT| tutorials directory::

      $ %HOMEPATH%\npby3\Scripts\activate
      $ cd tutorials

#. You should now be able to run a |NUMBAT| calculation::

      $ make tut01

#. If this program runs without error, congratulations! You are now ready  to proceed to the next chapter to begin using |NUMBAT|.  If not, please check the instructions above again, and if still stuck, follow the instructions under :ref:`sec-helpinstall-label` to seek assistance.



Creating a self-contained command terminal
---------------------------------------------

If you plan to build the fortran code frequently,
both the python and Intel oneAPI paths need to be set up in your terminal.
Doing this manually requires typing::

    $ %HOMEPATH%\npby3\Scripts\activate
    $ c:\Program Files (x86)\Intel\oneAPI\setvars.bat

This quickly becomes tedious. To automatically activate your python environment and ensure all other necessary paths are correctly setup, it is helpful to create a dedicated launcher for the desktop that executes the required commands on first opening the terminal.

Here is a procedure for doing this

  #. Copy the launcher file  ``numbat_cmd.bat`` to your |NUMBAT| root directory::

       $ copy <NumBAT>\share\numbat_cmd.bat <NumBAT_BASE>

  #. Create a desktop shortcut to the Windows command terminal by using File Explorer to open the folder ``c:\Windows\System32``, typing ``cmd.exe`` in the search box at top right, and then right-clicking *Send to Desktop*.

  #. Right click on the new shortcut and open its *Properties* dialog.

  #. Select the *General* tab and change the name field at the top to *NumBAT Terminal*.

  #. Select the *Shortcut* tab and change the *Target* field to ``%windir%\System32\cmd.exe "/K" %HOMEPATH%\numbat\numbat_cmd.bat``

  #. Click the *Change Icon* button and select the |NUMBAT| icon at ``<NumBAT>\docs\source\numbat.ico``.



.. _sec-troubleshooting-windows-label:


Troubleshooting a Windows installation
-------------------------------------------

#. My build of |NUMBAT| completes but the  ``nb_install_tester.py`` program complains the |NUMBAT| fortran ``nb_fortran.pyd`` dynamically linked library (DLL) can't be loaded.

  This is usually due to another DLL required by |NUMBAT| not being found, either because it is in an unexpected location or missing altogether.  This can be a little painful to diagnose. The following procedure is relatively straightforward.

  #. Download the *Dependencies* tool available as a zip file install from  `github <https://github.com/lucasg/Dependencies>`_. This tool displays all the DLL dependencies of a given file and whether or not they have been located in the file system. Extract the zip file to a folder named ``dependencies`` in ``<NumBAT_BASE>\usr_local\packages``.

  #. Now we can apply the tool to the |NUMBAT| python dll.

    Start the ``DependenciesGui.exe`` tool::

      $ <NUMBAT_BASE>\usr_local\packages\dependencies\DependenciesGUI.exe

    Browse to your |NUMBAT| folder and open ``backend\fortran\nb_fortran.pyd``.


    Examine the output and note any red highlighted entries. These indicate required DLLs that have not been found.  If one or more such lines appear, read through the install instructions again and ensure that any commands to copy DLLs to particular locations have been executed.

  #. Alternatively, you can try the command line version of this tool::

      $ <NUMBAT_BASE>\usr_local\packages\dependencies\Dependencies.exe -depth 1 -modules nb_fortran.pyd



Installing the Linux version via a Virtual Machine
------------------------------------------------------

Another way to run |NUMBAT| on  Windows or MacOS is by installing Ubuntu as
a virtual machine using either `Microsoft Hyper-V
<https://wiki.ubuntu.com/Hyper-V>`_ or `Oracle Virtual Box
<https://ubuntu.com/tutorials/how-to-run-ubuntu-desktop-on-a-virtual-machine-using-virtualbox#1-overview>`_, or a similar tool on MacOS.

Then |NUMBAT| can be installed using exactly the same procedure as described
above for standard Linux installations.  It is also possible to build |NUMBAT|
using the `Windows Subsystem for Linux
<https://msdn.microsoft.com/en-au/commandline/wsl/install_guide>`_, but dealing
with installing the additional required packages may be quite painful.



Building the documentation
===============================

You can rebuild the documentation you are currently reading by moving into the ``<NumBAT>/docs`` directory and running either ``make html`` or ``make latexpdf``.

In each case, the output is placed in the ``<NumBAT>/docs/build`` directory.


Note however most of the figures will only be available after you have run all the example problems in the  ``tutorial`` ``lit_ex`` and ``JOSAB_tutorial`` directories.
This can be done by running ``make`` in each of those directories. Be aware that  some of these problems are quite large and may require some time to complete depending on your computer's performance.



.. _sec-helpinstall-label:

Seeking help with building |NUMBAT|
=====================================

If you are having trouble building |NUMBAT| or are experiencing crashes, we will do our best to assist you.

Before writing for help (see contact details in :ref:`chap-intro-label`), please do the following:

    #. Download the latest version of |NUMBAT| from github and ensure the problem remains.

    #. If on Linux or MacOS, run the script ``./backend/nb_runconfig_test.sh`` from the main |NUMBAT| directory.

       This will create the file ``./nb_buildconfig.txt`` in the main directory with useful details about your build configuration and environment.

       **If on Windows**, do the same using the file ``.\backend\nb_runconfig_test.bat`` instead.

    #. In your email, indicate the date on which you last downloaded |NUMBAT| from github.

    #. Attach the ``nb_buildconfig.txt`` file.

    #. If the problem is a crash while running a |NUMBAT| script, please attach the python script file and a description of how to observe the problem.

    #. If on Linux or MacOS, follow the instructions under :ref:`sec-troubleshooting-label` to generate a debugging trace with GDB and send a screen shot of the trace.
