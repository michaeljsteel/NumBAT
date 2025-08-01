
.. include:: numbatdefs.txt

.. _chap-pythonapi-label:


*********************
Python Interface API
*********************



This chapter provides an auto-generated summary of the NumBAT Python API.

The API consists of several core modules\:

  - ``numbat``, for creating the top-level |NUMBAT| application;
  - ``materials``, for defining waveguide materials and their properties;
  - ``voigt``, for manipulating tensor quantities.
  - ``structure``, for constructing waveguides from materials;
  - ``modecalcs``, for the core calculation of electromagnetic and acoustic modes;
  - ``integration``, for performing calculations relating to SBS gain;
  - ``plotting``, for creating output plots of modes and gain functions.

For the most part, users should only encounter the ``numat`` and ``integration`` modules.

numbat module
=============

The ``numbat`` module contains the :class:`NumBATApp` through which the bulk of the |NUMBAT| API is accessed.

Creating a :class:`NumBATApp` object is normally the first main step in a |NUMBAT| script.

.. automodule:: numbat
    :members:
    :undoc-members:
    :show-inheritance:

materials module
================

The ``materials`` module provides functions for specifying all relevant optical and elastic properties of waveguide materials.

The primary class is :class:`materials.Material` however users will rarely use this class directly. Instead, we generally specify material properties by writing new  ``.json`` files stored in the folder ``backend/materials_data``. Materials are then loaded to build waveguide structures using the function :func:`materials.make_material`.

.. automodule:: materials
    :members:
    :undoc-members:
    :show-inheritance:

voigt module
================

The ``voigt`` module provides functions for performing a number of tensor operations with regular and Voigt style tensors.

.. automodule:: voigt
    :members:
    :undoc-members:
    :show-inheritance:

modes module
================

The ``modes`` module defines the main classes for displaying and interrogating optical and elastic modes.

.. automodule:: modes
    :members:
    :undoc-members:
    :show-inheritance:

structure module
==============

The ``structure`` module provides functions for defining and constructing waveguides.
The diagrams in Chapter 2 can be used to identify which parameters (``slab_a_x, slab_c_y, material_d`` etc) correspond to each region.

.. automodule:: structure
    :members:
    :undoc-members:
    :show-inheritance:


modecalcs module
=================

The ``modecalcs`` module is responsible for the core engine to construct and solve the optical and elastic finite-element problems.

.. automodule:: modecalcs
    :members:
    :undoc-members:
    :show-inheritance:

integration module
==================

The ``integration`` module is responsible for calculating gain and loss information from existing mode data.

.. automodule:: integration
    :members:
    :undoc-members:
    :show-inheritance:

plotting module
=================

The ``plotting`` module is responsible for generating all standard graphs.

.. automodule:: plotting
    :members:
    :undoc-members:
    :show-inheritance:



