.. include:: ../numbatdefs.txt


.. _sec-wininstall-label:

Installing on Windows
================================

There are three methods for installing |NUMBAT| on Windows:

#. Use the pre-built binary installer.

   This is usually the easiest method to get going quickly. However, the pre-built installer is updated less frequently than the core github source code tree.

#. Build the entire system from scratch including all the mathematical support libraries.

   This method provides the most flexibility but does involve a significant number of steps. This is recommended for users with some experience of building code.

#. Build the core NumBAT code but obtain the supporting mathematical libraries with a pre-built installer.

   This is a good compromise for most users, allowing immediate access to new features in the github tree.


For all three methods, there are some common steps, including setting up the python environment and installing the open-source GMsh tool.
These common steps are described in the next section and should be performed first. Then proceed to the section corresponding to your preferred installation mode.

Common steps - Setting up the |NUMBAT| Python environment
-----------------------------------------------------------------

|NUMBAT| is entirely driven from Python, either as scripts or Jupyter notebooks, so a working Python installation is a basic requirement.

You can use an existing Python installation if you have a sufficiently recent version (>3.11.0) installed, or download a new Python installer. Note that |NUMBAT| is a self-contained tree and will not add any files to your Python installation.

**Note: if you are using the pre-built NumBAT installer, for binary compatibility you *must* use Python version 3.13 and numpy version 2.2.**

#.  If you do not have a suitable current Python, download and run the installer for a recent version from the `Python website <https://www.python.org/downloads/windows/>`_.

    By default, this will install Python in the directory ``%HOMEPATH%\AppData\Local\Programs\Python\<PythonVer>``, say ``%HOMEPATH%\AppData\Local\Programs\Python\Python312``.


#.  Create a python virtual environment for working with |NUMBAT|.

    You can use any name and location for your environment.

    To specify a virtual environment tree called `nbpy3`, open a command prompt (or your favourite Windows terminal app) from the Start Menu and  enter ::

        $ %HOMEPATH%\AppData\Local\Programs\Python\Python312\python.exe -m venv nbpy3


#.  Activate the new python virtual environment ::

        $ %HOMEPATH%\npby3\Scripts\activate

#.  Install the necessary python tools and libraries ::

        $ pip3 install numpy  matplotlib scipy psutil ninja meson==1.4.1

    Note that at last check, the most recent meson (1.5.0) is broken and we specify the earlier 1.4.1 version.


#.  Finally, we will also need the Gnu ``make`` tool. This can be installed by typing ::

       $ winget install ezwinports.make

    and then starting a new terminal (so that the PATH variable is updated to find ``make.exe``.)


Now you can proceed to install |NUMBAT| using one of the following methods:

#. Method 1: Using the complete pre-built installer in :ref:`sec-wininstallmeth1-label`
#. Method 2: Building all components from source in :ref:`sec-wininstallmeth2-label`
#. Method 3: Building core |NUMBAT| with pre-built mathematical libraries in :ref:`sec-wininstallmeth3-label`


.. _sec-wininstallmeth1-label:

Installing |NUMBAT| Method 1: pre-built Windows installer
-----------------------------------------------------------------
*Note*: The installer will allow you to run |NUMBAT| problems and make changes to the code on the Python side only.
You will not be able to make changes to the Fortran finite element code. (For most people, this is completely fine.)
The installer is not updated as frequently as the main source tree.

**Note: the pre-built installer will only work with Python version 3.13. If you do not have this version installed, please follow the instructions for installing Python earlier in this section.**

We will actually run two installers: one for the free GMsh mesh generation tool, and one for |NUMBAT| itself.

.. _sec-wininstallfoldergmsh-label:

Setting up NumBAT's folder structure and installing gmsh
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


#. Choose a location for the base folder for building |NUMBAT| and supporting libraries, say ``c:\Users\<myname>\numbat``, which we will refer to as ``<NumBAT_BASE>``.

   Create this folder using either Windows Explorer or a command prompt.

#. Use the Start Menu to open a Windows terminal.

#. In the terminal window, change to the ``<NumBAT_BASE>`` directory, then execute the following commands ::

      $ mkdir nb_releases
      $ mkdir usr_local
      $ mkdir usr_local\packages

   Your |NUMBAT| trees will be stored inside the `nb_releases` directory.


#. Download the `Windows build of Gmsh <https://gmsh.info>`_ and unzip the tree into ``<NumBAT_Base>\usr_local\packages\gmsh``.  The Gmsh executable should now be at ``<NumBAT_Base>\usr_local\packages\gmsh\gmsh.exe``.


Installing the NumBAT binary installer
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#.  Download the `Windows installer <https://github.com/michaeljsteel/NumBAT>`_ from the |NUMBAT| github page. The link to the installer can be found at the bottom of the *Readme* section and also under the *Releases* heading in the right-hand column of the page.

    Be sure to download the `numbat_complete_installer_win64.exe` file and *not* the mathlibs installer.

#.  Run the installer specifying a *new* install directory inside`` the ``<NumBAT_BASE>\nb_releases`` folder.

    We recommend including the |NUMBAT| release number for the install directory name, for example ``<NumBAT_BASE>\nb_releases\NumBAT-2.1.3``.


Testing the installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#.  In your Windows terminal, change directory to the newly installed |NUMBAT| folder.

#.  Activate your python environment and then move to the |NUMBAT| examples/tutorials directory ::

      $ %HOMEPATH%\npby3\Scripts\activate
      $ cd examples\tutorials

    Remember, you must activate a python environment that uses Python 3.13.

#.  You should now be able to run a |NUMBAT| calculation ::

      $ make tut01

#.  If this program runs without error, congratulations! You are now ready  to proceed to the next chapter to begin using |NUMBAT|.

    If not, please check the instructions above again, and if still stuck, please read the :ref:`troubleshooting <sec-troubleshooting-windows-label>` section to attempt to diagnose the problem. Then follow the instructions at :ref:`sec-helpinstall-label` to seek assistance.


.. _sec-wininstallmeth2-label:

Method 2: Building all components from source
----------------------------------------------------------

The Windows version of |NUMBAT| is built using the native Windows toolchain including Visual Studio and the Intel Fortran compiler. Since a standard Windows installation is not set up as a development environment, there are a number of steps and tools required. If you are new to building software,
you might like to begin with Method 3.


.. _sec-winbuildtoolsmeth2-label:

Windows build tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We need to install a number of compilers and math libraries.
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
    This is a widely-used source control tool that we use to download |NUMBAT| and some other tools.


    Download and run the  `latest Git for Windows release <https://git-scm.com/download/win>`_, accepting all defaults.

    Some users may prefer to use a graphical interface such as `GitHub Desktop <https://desktop.github.com/>`_. This is fine too.

#. *GNU make*
      We will need the Gnu ``make`` tool. This can be installed by typing ::

         $ winget install ezwinports.make

      and then starting a new terminal (so that the PATH variable is updated to find ``make.exe``.)


Building the supporting math libraries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
We can now build the supporting libraries, before proceeding to build |NUMBAT| itself.

#. To build Arpack and SuiteSparse libraries, we need the *CMake* tool.

    This is a cross-platform build tool we will need for building some of the libraries.

    Download and run the  `latest release <https://cmake.org/download/>`_ accepting all defaults.

#. Choose a location for the base directory for building |NUMBAT| and supporting libraries, say ``c:\Users\<myname>\numbat``, which we will refer to as ``<NumBAT_BASE>``.

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

4. In the pull-down menu in the ribbon, select the *Release* build. Then from the *Build* menu select the *Build solution* option. This will take a few minutes.

5. Return to the command terminal and  cd to ``<NumBAT_BASE>\usr_local``. Then execute the following commands::

    $ copy packages\arpack-ng\build\Release\* lib
    $ copy packages\arpack-ng\ICB\*.h include


At long last, we are ready to build |NUMBAT| itself.


.. _sec-wininstall-numbatitself-label:

Building |NUMBAT|
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
We are now ready to build and test |NUMBAT| itself.

#. Your command prompt needs to have your |NUMBAT| python environment activated.

   If you have opened a new prompt since setting up the python environment, enter ::

       $ %HOMEPATH%\npby3\Scripts\activate

#. We need two additional python modules to run the build process ::

     $ pip3 install ninja meson==1.4.1

   Note that at last check, the most recent meson (1.5.0) is broken and we specify the earlier 1.4.1 version.

#. In your command prmopt, move to your root ``<NumBAT_BASE>`` directory and then to the |NUMBAT| folder itself::

    $ cd <NumBAT_BASE>
    $ cd nb_releases\nb_latest

   From this point, we refer to the current directory as ``<NumBAT>``.

   In other words, ``<NumBAT> = <NumBAT_BASE>\nb_releases\nb_latest``.

#. Setup the environment variables for the Intel compiler::

    $  "c:\Program Files (x86)\Intel\oneAPI\setvars.bat"


#. Move to the ``<NumBAT>\backend\fortran`` directory and open the file ``meson.options`` in a text editor. Check the values of the options in the ``Windows`` section, particularly the value for ``windows_dir_nb_usrlocal`` and change any of the paths in the ``value`` fields as required.

#. To initiate the build, enter ::

     $ make win

   This should take 2 to 3 minutes.


#. If all is well, this will run to completion. If you encounter errors, please check that all the instructions above have been followed accurately. If you are still stuck, see :ref:`sec-troubleshooting-windows-label` for further ideas.

#. If you hit a compile error you can't easily resolve, please see the instructions at :ref:`sec-helpinstall-label` on how to seek help.

#. Copy the .dlls generated earlier to this directory::

    $ copy ..\..\..\..\usr_local\lib\*.dll .

#. Also, copy the following Intel oneAPI .dlls to this directory::

    $ copy "c:\Program Files (x86)\Intel\oneAPI\mkl\latest\bin\mkl_rt.2.dll" .
    $ copy "c:\Program Files (x86)\Intel\oneAPI\mkl\latest\bin\mkl_intel_thread.2.dll" .
    $ copy "c:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin\svml_dispmd.dll" .
    $ copy "c:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin\libmmd.dll" .
    $ copy "c:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin\libifcoremd.dll" .
    $ copy "c:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin\libifportMD.dll" .

#. At this point, we are ready to test the installation with a short script that checks whether required applications and libraries can be found and loaded. Perform the following commands::

    $ cd <NumBAT>/backend
    $ python .\nb_install_tester.py

#. If this program runs without error, congratulations! You are now ready  to proceed to the next chapter to begin using |NUMBAT|.  If not, please see the suggestions at :ref:`sec-troubleshooting-windows-label`. But before moving on, please read the  section :ref:`sec-winbuild-cmdterminal-label` on creating a specalised |NUMBAT| command terminal.

#. Once again, if you run into trouble, please don't hesitate to get in touch for help using the instructions at :ref:`sec-helpinstall-label`. Please do send all the requested information, as it usually allows us to solve your problem faster.



.. _sec-wininstallmeth3-label:

Method 3: Building core |NUMBAT| with pre-built mathematical libraries
------------------------------------------------------------------------

#.  Set up the Windows build environment.

    Follow the instructions in section :ref:`sec-winbuildtoolsmeth2-label` to install *Visual Studio*, the *Intel OneAPI Base and HPC Toolkits*, *Git* and *GNU Make*.


#.  Download the `Windows maths support libraries for NumBAT <https://github.com/michaeljsteel/NumBAT>`_ installer from the |NUMBAT| github page. The link to the installer can be found at the bottom of the *Readme* section and also under the *Releases* heading in the right-hand column of the page.

    Be sure to download the `numbat_mathlibs_installer_win64.exe` file and *not* the complete installer.

    Run the installer choosing a base install directory of your choice, say ``c:\Users\<myname>\numbat``, which we will refer to as the ``<NumBAT_BASE>`` folder.

#.  Download the `Windows build of Gmsh <https://gmsh.info>`_ and unzip the tree into ``<NumBAT_Base>\usr_local\packages\gmsh``.  The Gmsh executable should now be at ``<NumBAT_Base>\usr_local\packages\gmsh\gmsh.exe``.

#.  Carry out the instructions in the section :ref:`sec-wininstall-numbatitself-label` to build and test the core |NUMBAT| code.

#.  You are now ready to start working with |NUMBAT|.

    Before moving on, to make updating |NUMBAT| simpler, please read the  section :ref:`sec-winbuild-cmdterminal-label` on creating a specalised |NUMBAT| command terminal.




.. _sec-winbuild-cmdterminal-label:

Creating a self-contained command terminal
----------------------------------------------------------

If you plan to build the fortran code frequently to keep up to date with changes in the source tree,
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

#.  My build of |NUMBAT| completes but the  ``nb_install_tester.py`` program complains the |NUMBAT| fortran ``nb_fortran.pyd`` dynamically linked library (DLL) can't be loaded.

    This is usually due to another DLL required by |NUMBAT| not being found, either because it is in an unexpected location or missing altogether.  This can be a little painful to diagnose. The following procedure is relatively straightforward.

    #.  Download the *Dependencies* tool available as a zip file install from  `github <https://github.com/lucasg/Dependencies>`_. This tool displays all the DLL dependencies of a given file and whether or not they have been located in the file system. Extract the zip file to a folder named ``dependencies`` in ``<NumBAT_BASE>\usr_local\packages``.

    #.  Now we can apply the tool to the |NUMBAT| python dll.

    Start the ``DependenciesGui.exe`` tool::

      $ <NUMBAT_BASE>\usr_local\packages\dependencies\DependenciesGUI.exe

    Browse to your |NUMBAT| folder and open ``backend\fortran\nb_fortran.pyd``.


    Examine the output and note any red highlighted entries. These indicate required DLLs that have not been found.  If one or more such lines appear, read through the install instructions again and ensure that any commands to copy DLLs to particular locations have been executed.

    #.  Alternatively, you can try the command line version of this tool ::

           $ <NUMBAT_BASE>\usr_local\packages\dependencies\Dependencies.exe -depth 5 -modules nb_fortran.pyd

  If you find a missing DLL by one of these methods, please :ref:`let us know <sec-helpinstall-label>`. It may suggest a problem with the pre-built installer or the documentation instructions.



Installing the Linux version via a Virtual Machine
------------------------------------------------------

Yet another way to run |NUMBAT| on  Windows or MacOS is by installing Ubuntu as
a virtual machine using either `Microsoft Hyper-V
<https://wiki.ubuntu.com/Hyper-V>`_ or `Oracle Virtual Box
<https://ubuntu.com/examples/tutorials/how-to-run-ubuntu-desktop-on-a-virtual-machine-using-virtualbox#1-overview>`_, or a similar tool on MacOS.

Then |NUMBAT| can be installed using exactly the same procedure as described
above for standard Linux installations.

As the native installation methods now work well on all platforms, this is not recommended for normal use.
