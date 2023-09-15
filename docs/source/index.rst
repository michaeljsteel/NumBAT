.. NumBAT documentation master file, created by
   sphinx-quickstart on Wed Oct 19 15:23:46 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. include:: numbatdefs.txt


Contents
====================================

.. toctree::
   :maxdepth: 2
   :numbered: 


==================
Introduction
==================

.. toctree::
    :maxdepth: 4

    intro


.. _chap-install-label:

==================
Installation
==================

.. toctree::
    :maxdepth: 4

    install

====================================
Basic Usage 
====================================

.. toctree::
    :maxdepth: 4

    usage

====================================
Tutorial
====================================

.. toctree::
    :maxdepth: 4

    tutorial


.. _chap-josab-label:

====================================
JOSA-B Tutorial Paper
====================================

.. toctree::
    :maxdepth: 4

    josab

.. _chap-literature-label:

====================================
Literature Examples
====================================

.. toctree::
    :maxdepth: 4

    lit_ex

.. _chap-pythonbackend-label:

==================
Python Interface
==================

.. toctree::
    :maxdepth: 4

    pythonapi

==================
Fortran Backend
==================

The intention of |NUMBAT| is that the Fortran FEM routines are essentially black boxes. 
They are called from ``mode_calcs.py`` and return the modal fields. 
However, there are a few important things to know about the workings of these routines.

.. toctree::
    :maxdepth: 4

    fem





Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

