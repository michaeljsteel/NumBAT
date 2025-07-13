
.. include:: numbatdefs.txt

.. _chap-install-label:

*******************
Installing NumBAT
*******************


This chapter provides instructions on installing |NUMBAT| on each platform.
Please email |NUMBAT_EMAIL| to let us know of any difficulties you encounter, or suggestions for improvements to the install procedure on any platform, but especially for MacOS or Windows.


While |NUMBAT| is developed on Linux, it can also be built natively on both MacOS and Windows. A pre-built executable version is also available for Windows, though it is updated less frequently than the main source code tree.
The Linux builds can also be run under virtual machines on MacOS and Windows if desired.

In all cases, the current source code for |NUMBAT| is hosted `here on Github <https://github.com/michaeljsteel/NumBAT>`_. Please always download the latest release from the github page.

You can now skip forward to the section for your operating system: :ref:`Linux <sec-linuxinstall-label>`, :ref:`MacOS <sec-macosinstall-label>` or :ref:`Windows <sec-wininstall-label>`.


**Installation instructions for your OS**


.. toctree::
    :maxdepth: 1

    install/install_linux
    install/install_macos
    install/install_windows


.. _sec-helpinstall-label:

Seeking help with building |NUMBAT|
=====================================

If you are having trouble building |NUMBAT| or are experiencing crashes, we will do our best to assist you.

Before writing for help (see contact details in :ref:`chap-intro-label`), please do the following:

    #. Download the latest version of |NUMBAT| from github and ensure the problem remains.

    #. If on Linux or MacOS, run the script ``./backend/nb_runconfig_test.sh`` from the main |NUMBAT| directory.

       This will create the file ``./nb_buildconfig.txt`` in the main directory with useful details about your build configuration and environment.

       *If on Windows*, do the same using the file ``.\backend\nb_runconfig_test.bat`` instead.

    #. In your email, indicate the date on which you last downloaded |NUMBAT| from github.

    #. Attach the ``nb_buildconfig.txt`` file.

    #. If the problem is a crash while running a |NUMBAT| script, please attach the python script file and a description of how to observe the problem.

    #. If you are on Linux and encounter a segfault or "Segmentation violation", it is especially helpful if you able to follow the instructions under :ref:`sec-troubleshooting-linux-label` to generate a debugging trace with GDB and send a screen shot of the trace.




Building the documentation
===============================

If you should want to, you can rebuild the documentation you are currently reading. You will need a working LaTeX installation to build the pdf version.


Steps for building the documentation
----------------------------------------

*   Ensure that all examples have been built.


    Most of the figures will only be available after you have run all the problems in the  ``tutorial``, ``lit_ex`` and ``josab_tutorial`` sub-directories of the ``examples`` directory have been run.
    This can be done by running ``make`` in each of those directories. Be aware that  some of these problems are quite large and may require some time to complete depending on your computer's performance.

*  Install documentation related python modules.

   Activate your normal |NUMBAT| python environment and then type ::

         $ pip3 install sphinx nbsphinx sphinx_subfigure sphinxcontrib-bibtex setuptools pandoc ipython pygments

*  Install Pandoc.

    You may also need to install the Pandoc package  for your distribution using your package manager or from `<https://pandoc.org/>`_ .

*  Build the docs (choosing whichever version you require) ::

      $ cd <numbat_base>/docs
      $ make html
      $ make latexpdf

*  The output will be found in the `<numbat_base>\docs\build` directory.


