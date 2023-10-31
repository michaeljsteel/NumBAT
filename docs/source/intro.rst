.. include:: numbatdefs.txt

.. _chap-intro-label:

****************************
Introduction to |NUMBAT|
****************************


.. role:: raw-math(raw)
    :format: latex html

.. figure:: NumBAT_logo.png
   :scale: 40 %

Introduction
================

|NUMBAT|, the Numerical Brillouin Analysis Tool, integrates electromagnetic and acoustic mode solvers to calculate the interactions of optical and acoustic waves in waveguides.

Goals
================
|NUMBAT| is designed primarily to calculate the optical gain response from
stimulated Brillouin scattering (SBS) in integrated waveguides. It uses finite element algorithms
to solve the electromagnetic and acoustic modes of a wide range of 2D waveguide structures. It
can account for photoelastic/electrostriction and moving boundary/radiation pressure effects, as well as
arbitrary acoustic anisotropy. 

|NUMBAT| also supports user-defined material properties and we hope its creation will drive a community-driven
set of standard properties and geometries which will allow all groups to test and validate each other's
work.

A full description of the |NUMBAT| physics and numerical algorithms  is available in the article B.C.P Sturmberg at al., "Finite element analysis of stimulated Brillouin scattering in integrated photonic waveguides", *J. Lightwave Technol.*  **37**, 3791-3804 (2019), 
available at `<https://dx.doi.org/10.1109/JLT.2019.2920844>`_.


|NUMBAT| is open-source software and the authors welcome additions to the code.  Details for how
to contribute are available in :ref:`sec-contribute-label`.

Citing |NumBAT|
===============
If you use |NumBAT| in published work, we would appreciate a citation 
to B.C.P Sturmberg at al.,
"Finite element analysis of stimulated Brillouin scattering in integrated
photonic waveguides", *J. Lightwave Technol.*  **37**, 3791-3804 (2019), 
available at `<https://dx.doi.org/10.1109/JLT.2019.2920844>`_ 
and `<https://arxiv.org/abs/1811.10219>`_,
and a link to the github page at `<https://github.com/michaeljsteel/NumBAT>`_.

Development team
================

|NUMBAT| was developed by Bjorn Sturmberg, Kokou Dossou, Blair Morrison, Chris
Poulton and Michael Steel in a collaboration between Macquarie University, the
University of Technology Sydney, and the University of Sydney.

We thank Christian Wolff, Mike Smith and Mikolaj Schmidt for contributions.

Contributing to NumBAT
================================
NumBAT is open source software licensed under the GPL with all source and documentation available
at `github.com <https://github.com/michaeljsteel/NumBAT.git>`_. We welcome additions to NumBAT code, documentation and the materials library. Interested users should fork the standard release from github and make a pull request when ready.  For major changes, we strongly suggest contacting the NumBAT team before starting work at ``michael.steel@mq.edu.au``.


Support
=============
Development of |NumBAT| has been supported in part by the
Australian Research Council under Discovery Projects 
DP130100832, DP160101691, DP200101893 and DP220100488.

Release notes
=============

Version 2.0
-----------

A number of API changes have been made in |NUMBAT| 2.0 to tidy up the interface and make plotting and analysis simpler and more powerful. 
You will need to make some changes to existing files to run in |NUMBAT| 2.0.  Your best guide to new capabilities and API changes is to look through the code in the tutorial examples.

Some key changes you will need to make are as follows:
 * The waveguide class ``Struct`` has been renamed to ``Structure`` 
 * The interface for creating materials has changed. You now call the ``materials. make_material(`` *name* ``)`` function. For example ``material_a = materials.make_material('Vacuum')``
 * To access an existing material in an  existing ``Struture`` object (usually in a variable called ``wguide`` use ``wguide.get_material(`` *label* ``)`` For example, ``mat_a = wguide.get_material('b')`` where the allowed labels are ``bkg`` and the letters ``a`` to ``r``.
 * The member name for refractive index in a ``Material`` object has changed from ``n`` to ``refindex_n``.
 * The member name for density in a ``Material`` object has changed from ``n`` to ``rho``.
 * Due to a change in parameters, the function ``plotting.gain_spectra`` is deprecated and replaced by ``plotting.plot_gain_spectra`` with the following changes:
      * The frequency arguments ``freq_min`` and ``freq_max`` should now be passed in units of Hz, not GHz.
      * The argument ``k_AC`` has been removed.
 * In all functions the parameter ``prefix_str`` has been renamed to ``prefix`` for brevity.


What does |NUMBAT| actually calculate? 
=======================================

|NUMBAT| performs three main types of calculations given a particular waveguide design: 
 * solve the electromagnetic modal problem  using the finite element method (FEM).
 * solve the elastic modal problem  using FEM.
 * calculate Brillouin gain coefficients and linewidths for a given triplet of two optical and one elastic mode, and use this to generate gain spectra.


Here we specify the precise mathematical problems been solved.
For further details, see  the |NUMBAT| paper in the Journal of Lightwave Technolgoy at
at `<https://dx.doi.org/10.1109/JLT.2019.2920844>`_.

1. Electromagnetic modal problem
----------------------------------


WRITE ME

2. Elastic  modal problem
----------------------------------

WRITE ME

3. SBS gain calculation  modal problem
-----------------------------------------

WRITE ME
