.. include:: ../numbatdefs.txt

.. _sec-macosinstall-label:

Installing on MacOS
================================

|NUMBAT| can also be installed on MacOS, though this is less tested than other versions.  We are keen to increase our support of the MacOS build. Any comments on difficulties and solutions will be appreciated, so please don't hesitate to get in touch using the directions at :ref:`sec-helpinstall-label`.

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

    $ pip3 install numpy matplotlib scipy psutil meson ninja gitpython

#. Check that the python installs work. ::

     $ python3.12
     >>> import matplotlib
     >>> import numpy
     >>> import scipy
     >>> import psutil

#. Move to the |NUMBAT| fortran directory::

   $ cd backend/fortran

#. Move to the ``<NumBAT>/backend/fortran/`` directory and open the file ``meson.options`` in a text editor. Check the values of the options in the ``MacOS`` section and change any of the paths in the ``value`` fields as required.

#. To start the build, enter::

   $ make mac

#. If all is well, this will run to completion. If you encounter errors, please check that all the instructions above have been followed accurately. If you are still stuck, see :ref:`sec-troubleshooting-linuxmacos-label` for further ideas.

#. If you hit a compile error you can't resolve, please see the instructions at :ref:`sec-helpinstall-label` on how to seek help.

#. Once the build has apparently succeeded, it is time to test the installation with a short script that tests whether required applications and libraries can be found and loaded. Perform the following commands::

    $ cd <NumBAT>/backend
    $ python ./nb_install_tester.py

#. If this program runs without error, congratulations! You are now ready  to proceed to the next chapter to begin using |NUMBAT|.

#. Once again, if you run into trouble, please don't hesitate to get in touch for help using the instructions at :ref:`sec-helpinstall-label`. Please do send all the requested information, as it usually allows us to solve your problem faster.


.. _sec-troubleshooting-macos-label:

Troubleshooting MacOS installs
-------------------------------------------
Performing a full build of |NUMBAT| and all its libraries from scratch is a non-trivial task and it's possible you will hit a few stumbles.
Here are a few traps to watch out for:

#. Please ensure to use relatively recent libraries for all the Python components. This includes using

   - Python: 3.11 or later
   - ``matplotlib``: 3.9.0 or later
   - ``scipy``:  1.13.0 or later
   - ``numpy``:  2.0 or later


#. Be sure to follow the instructions above about setting up the virtual environment for |NUMBAT| excusively. This will help prevent incompatible Python modules being added over time.

#. In general, the GCC build is more tested and forgiving than the build with the Intel compilers and we recommend the GCC option.  However, we do recommend using the Intel OneAPI math libraries as described above. This is the easiest way to get very high performance LAPACK and BLAS libraries with a well-designed directory tree.

#. If you encounter an error about "missing symbols" in the NumBAT fortran module, there are usually two possibilities:

   - A shared library (a file ending in ``.so``) is not being loaded correctly because it can't be found in the standard search path. To detect this, run ``ldd nb_fortran.so`` in the ``backend/fortran`` directory and look for any lines containing ``not found``. (On MacOS, use ``otools -L nb_fortran.so`` instead of ``ldd``.)


     You may need to add the directory containing the relevant libraries to your ``LD_LIBRARY_PATH`` in your shell setup files (eg. ``~/.bashrc`` or equivalent).

   - You may have actually encountered a bug in the |NUMBAT| build process. Contact us for assistance as described in the introduction.

