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

.. _sec-contribute-label:

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
 * On Linux, the fortran Makefile is now designed to work with a virtual environment python to avoid dependencies on your system python.
 * There is a new core |NUMBAT| module ``numbat`` that should be imported before any other |NUMBAT| modules.
 * It should no longer be necessary to import the ``object`` or ``Numbat`` (note different case) modules.
 * The first call to any |NUMBAT| code should be to create a |NUMBAT| application object by calling ``nbapp = numbat.NumBATApp()``.
 * The default output prefix can now be set as an argument to ``numbat.NumBATApp()``. All output can be directed to a sub-folder of the starting directory with a second argument: ``nbapp = numbat.NumBATApp('tmp', 'tmpdir')``.
 * The waveguide class ``Struct`` has been renamed to ``Structure``.
 * A waveguide is now constructed using ``nbapp.make_waveguide`` rather than ``object.Structure``.
 * The interface for creating materials has changed. You now call the ``materials. make_material(`` *name* ``)`` function. For example ``material_a = materials.make_material('Vacuum')``
 * To access an existing material in an  existing ``Structure`` object (say, in a variable called ``wguide``) use ``wguide.get_material(`` *label* ``)`` For example, ``mat_a = wguide.get_material('b')`` where the allowed labels are ``bkg`` and the letters ``a`` to ``r``.
 * The member name for refractive index in a ``Material`` object has changed from ``n`` to ``refindex_n``.
 * The member name for density in a ``Material`` object has changed from ``n`` to ``rho``.
 * Due to a change in parameters, the function ``plotting.gain_spectra`` is deprecated and replaced by ``plotting.plot_gain_spectra`` with the following changes:
      * The frequency arguments ``freq_min`` and ``freq_max`` should now be passed in units of Hz, not GHz.
      * The argument ``k_AC`` has been removed.
 * In all functions the parameter ``prefix_str`` has been renamed to ``prefix`` for brevity. Using the default output settings in ``NumBATApp()``, these should be rarely needed.
 * All waveguides are now specified as individual plugin classes in the files ``backend/msh/user_waveguides.json`` and ``backend/msh/user_meshes.py``.  These files provide useful examples of how to design and load new waveguide templates. See the following chapter for more details.


What does |NUMBAT| actually calculate?
=======================================

|NUMBAT| performs three main types of calculations given a particular waveguide design:
 * solve the electromagnetic modal problem  using the finite element method (FEM).
 * solve the elastic modal problem  using FEM.
 * calculate Brillouin gain coefficients and linewidths for a given triplet of two optical and one elastic mode, and use this to generate gain spectra.


Here we specify the precise mathematical problems been solved.
For further details, see  the |NUMBAT| paper in the Journal of Lightwave Technology at
at `<https://dx.doi.org/10.1109/JLT.2019.2920844>`_.

Electromagnetic modal problem
----------------------------------
The electromagnetic wave problem is defined by the vector wave equation

.. math::

   - \nabla \times (\nabla \times {\vec E}) + \omega^2 \epsilon_0 \, \epsilon_r(x,y) \vec E =0,

where the *real-valued* electric field has the following form for modal propagation along :math:`z`:

CHECK THE FACTOR of HALF here.

.. math::
   :nowrap:

   \begin{align*}
   \vec E(x,y,z,t) = & \frac{1}{2} \left ( {\vec \cal E}(\vec r) e^{- i  \omega t } + {\vec \cal E}^* (\vec r) e^{- i  \omega t } \right) \\
                   = & \frac{1}{2} \left ( a(z) \vec e(x,y) e^{i (kz-\omega t) } + a^*(z) \vec e(x,y) e^{-i (kz-\omega t) } \right),
   \end{align*}


in terms of the complex field amplitude :math:`\mathcal{E}(\vec r)`, the mode profile :math:`\vec e(x,y)` and the complex slowly-varying envelope function :math:`a(z)`.

By Faraday's law the complex magnetic field amplitude is given by

.. math::

   {\vec {\cal B} } =  \frac{1}{i \omega} \nabla \times {\vec {\cal E}}.



Elastic  modal problem
----------------------------------

The elastic modal problem is defined by the wave equation

.. math::

   \nabla \cdot \bar{T} + \omega^2 \rho(x,y) \vec U = 0,

where :math:`\vec u` is the elastic displacement and :math:`\bar{T}=\mathbf{c}(x,y) \bar{S}` is the stress tensor,
defined in terms of the stiffness :math:`\mathbf{c}` and the strain tensor
:math:`\bar{S}=S_{ij} = \frac{1}{2}(\frac{\partial U_i}{\partial r_j} + \frac{\partial U_j}{\partial r_i})`.

The displacement has the modal propagation form

.. math::

   \vec U = b(z) \, \vec u(x,y) e^{i (qz-\Omega t) } + b^*(z) \, \vec u(x,y)^* e^{-i (qz-\Omega t) } ,


For details on how these problems are framed as finite element problems, we refer to `<https://dx.doi.org/10.1109/JLT.2019.2920844>`_.

Modal properties
-----------------------

For propagation in a given mode :math:`\vec e_n` or :math:`\vec U_n`, the optical (:math:`o`)
and elastic (:math:`a`) energy fluxes in Watts and linear energy densities in (J/m) are given by the
following expressions

.. math::

   {\cal P}_n^{(o)}        & = 2 \mathrm{Re} \int \mathrm{d}^2 r \, \hat{z} \cdot (\vec e_n^*(x,y) \times \vec h_n(x,y)), \\
   {\cal E}_n^{(o)} & = 2 \epsilon_0 \int \mathrm{d}^2 r \, \epsilon_r(x,y) |\vec e_n(x,y)|^2 \\
   {\cal P}_n^{(a)} & =  \mathrm{Re} \int \mathrm{d}^2 r \, (-2i\Omega) \sum_{jkl} c_{zjkl}(x,y) u^*_{mj}(x,y) \partial_k u_{ml}(x,y) \\
   {\cal E}_n^{(a)} & = 2 \Omega^2 \int_A \mathrm{d}^2 r \, \rho(x,y) |\vec u_n(x,y)|^2.

For fields with slowly-varying optical and elastic amplitudes :math:`a_m(z)` and :math:`b_m(z)`, the total
carried powers are

.. math::
   P^{(o)}(z)  & = \sum_m |a_n(z)|^2 {\cal P}_n^{(o)} \\
   P^{(a)}(z)  & = \sum_m |b_n(z)|^2 {\cal P}_n^{(a)}.

Note that in this convention, the amplitude functions :math:`a_m(z)` and :math:`b_m(z)`
are dimensionless and the units of the fields live in the modal functions
:math:`\vec e_n, \vec h_n, \vec U_n`.


SBS gain calculation  modal problem
-----------------------------------------

The photoelastic and moving boundary couplings in J/m are given by

.. math::

   Q^{\mathrm{(PE)}} & = - \epsilon \int_A  \mathrm{d}^2 r \, \sum_{ijkl} \epsilon_r^2  \,
   e_i^{(s)*}  \, e_j^{(p)}  \, p_{ijkl}   \, \partial_k u_l^* \\
   Q^{\mathrm{(MB))}} & = \int_{\cal C}  \mathrm{d} {\vec r} \, (\vec u^* \cdot \hat{n})
   \times \\
   & ~~~~
   \left [
   (\epsilon_a - \epsilon_b) \epsilon_0 (\hat{n} \times \vec u^{(s)})^* \cdot (\hat{n} \times \vec e^{(p)})
   -
   (\epsilon_a^{-1} - \epsilon_b^{-1}) \epsilon_0^{-1} (\hat{n} \cdot \vec d^{(s)})^*
          \cdot (\hat{n} \cdot \vec d^{(p)})
   \right]


Then the peak SBS gain :math:`\Gamma` is given by

.. math::
   \Gamma = \frac{2\omega \Omega}{\alpha_t} \frac{|Q_\mathrm{tot}|^2}{P^{(s)}P^{(p)}{\cal E}^{(a)}},

where the total SBS coupling is :math:`Q_\mathrm{tot} = Q^{(\mathrm{PE})} + Q^{(\mathrm{MB})}`.

Here :math:`\alpha_t` is the temporal elastic loss coefficient in :math:`\mathrm{s}^{-1}`.
It is related to the spatial attenuation coefficient by
:math:`\alpha_s = \alpha_t /v_{\mathrm{p}}^{(\mathrm{a})}` with :math:`v_{\mathrm{p}}^{(\mathrm{a})}` being the elastic phase velocity.

In a backward SBS problem, where there is genuine gain in Stokes optical field propagating in the negative
:math:`z` direction, its power evolves as

.. math::

   P^{\mathrm{(s)}}(z) = P_\mathrm{in}^{\mathrm{(s)}} e^{-\Gamma z}.


SBS equations of motion
-----------------------------------------

With the above conventions, the dynamical equations for the slowly-varying amplitudes are

..  math::

    \frac{1}{v_g^{p}} \frac{\partial }{\partial t} a_p +   \frac{\partial }{\partial z} a_p & = -i \frac{\omega_p Q_{\mathrm{tot}}}{{\cal P}_o} a_s b     \\
    \frac{1}{v_g^{s}} \frac{\partial }{\partial t} a_s -    \frac{\partial }{\partial z} a_s & = i  \frac{\omega_s Q_{\mathrm{tot}}^*}{{\cal P}_o} a_p b^*   \\
    \frac{1}{v_g^{a}} \frac{\partial }{\partial t} b  +    \frac{\partial }{\partial z} b
    + \frac{\alpha_s}{2} b & =  -i  \frac{\Omega Q_{\mathrm{tot} }^*}{{\cal P}_a}  a_p a_s^*



Here we've chosen the group velocities to be positive and included the propagation direction explicitly.


Connect these to output quantities from code
-------------------------------------------------


Equivalent forms of equations
-----------------------------------

TODO: show forms without the normalisation energies and with cubic style effective area.

Compare to some fiber literature and the hydrodynamic reprn.

