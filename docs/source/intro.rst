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


|NUMBAT|, the *Numerical Brillouin Analysis Tool*, is a software tool integrating electromagnetic and acoustic mode solvers to calculate the interactions of optical and acoustic waves in waveguides.
Most notably, this includes Stimulated Brillouin Scattering (SBS) frequency shifts and optical gains.

This chapter provides some background on the capabilities and techniques used in |NUMBAT|.
If you would like to get straight to computations,

jump ahead to the installation and setup instructions in :ref:`chap-install-label`.

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

A full description of |NUMBAT|'s physical and numerical algorithms is available in the article B.C.P Sturmberg at al., "Finite element analysis of stimulated Brillouin scattering in integrated photonic waveguides", *J. Lightwave Technol.*  **37**, 3791-3804 (2019),
available at `<https://dx.doi.org/10.1109/JLT.2019.2920844>`_.


|NUMBAT| is open-source software and the authors welcome additions to the code.  Details for how
to contribute are available in :ref:`sec-contribute-label`.

Development team
================

|NUMBAT| was developed by Bjorn Sturmberg, Kokou Dossou, Blair Morrison, Chris
Poulton and Michael Steel in a collaboration between Macquarie University, the
University of Technology Sydney, and the University of Sydney.

We thank Christian Wolff, Mike Smith and Mikolaj Schmidt for contributions.


Citing |NumBAT|
===============
If you use |NumBAT| in published work, we would appreciate a citation
to B.C.P Sturmberg at al.,
"Finite element analysis of stimulated Brillouin scattering in integrated
photonic waveguides", *J. Lightwave Technol.*  **37**, 3791-3804 (2019),
available at `<https://dx.doi.org/10.1109/JLT.2019.2920844>`_
and `<https://arxiv.org/abs/1811.10219>`_,
and a link to the github page at `<https://github.com/michaeljsteel/NumBAT>`_.


.. _sec-contribute-label:

Contributing to NumBAT
================================
NumBAT is open source software licensed under the GPL with all source and documentation available
at `github.com <https://github.com/michaeljsteel/NumBAT.git>`_. We welcome additions to NumBAT code, documentation and the materials library. Interested users should fork the standard release from github and make a pull request when ready.  For major changes, we strongly suggest contacting the NumBAT team before starting work at ``michael.steel@mq.edu.au``.

About our mascot
================================
The **numbat** (*Myrmecobius fasciatus*) is a delightful insect-eating marsupial

from Western Australia, of which it is the official state animal.
It has two other common names, *noombat* in the Nyungar language,

and *walpurti* in the Pitjantjatjara language.

As a carnivorous marsiupial, they belong to the order Dasyuromorphia, closely related to quolls and
the famed thylacines which had similar markings on their lower back.
Once found across southern Australia, numbats are now confined to small local groups

in Western Australia and the species has Endangered status.

.. figure:: ./numbat_face.jpg
   :scale: 40 %

   A numbat at Perth zoo in 2010. `(Creative commons) <https://commons.wikimedia.org/wiki/File:Numbat_Face.jpg>`_.

Apart from the distinctive striped back (which we like to think of as an acoustic wave made flesh),
numbats have a number of unique properties. They are the only fully diurnal marsupial.

They are insectivores and eat exclusively termites, perhaps 20000 each day!

To find out how you can support the care and revitalisation of this beautiful animal, check out the work at

`projectnumbat <numbat.org.au>`_ and the `Australian Wildlife Conservancy <https://www.australianwildlife.org/wildlife/numbat>`_.


Acknowledgements
===================

Development of |NumBAT| has been supported in part by the
Australian Research Council under Discovery Projects
DP130100832, DP160101691, DP200101893 and DP220100488.


Suggested reading on SBS and opto-elastic interactions in nanophotonics
===============================================================================

A very extensive literature on SBS and related effects in nanophotonics has
arisen over the period since 2010. Here we provide a few suggestions for
entering and navigating that literature.


Books
--------------------
The centenary of Brillouin scattering was marked with the publication of a two-volume book
featuring contributions from many of the leading researchers in the field.
These books provide detailed background and the history, theory, observation and application
of Brillouin scattering in guided wave systems.

#. *Brillouin Scattering, Parts 1 and 2*, eds: B.J. Eggleton, M.J. Steel, C.G. Poulton, (Academic, 2022). https://doi.org/10.1016/S0080-8784(22)00024-2

Reviews
--------------------
There are several excellent reviews covering the theory and experiment of Brillouin photonics.


#. A. Kobyakov, M. Sauer, and D. Chowdhury, "Stimulated Brillouin scattering in optical fibers," *Adv. Opt. Photon.* **2**, 1-59 (2010). https://doi.org/10.1364/AOP.2.000001

#. B.J. Eggleton, C.G. Poulton, P.T. Rakich, M.J. Steel and G. P. Bahl, "Brillouin integrated photonics," *Nat. Photonics* **13**, 664–677 (2019). https://doi.org/10.1038/s41566-019-0498-z

#. G.S. Wiederhecker, P. Dainese, T.P. Mayer Alegre, "Brillouin optomechanics in nanophotonic structures," *APL Photonics* **4**, 071101 (2019).
   https://doi.org/10.1063/1.5088169


Theoretical development
-------------------------
The following papers feature more depth on the theory of SBS in waveguides.
Chapters 2 and 3 of the *Brillouin Scattering* book listed above are also thorough
resources for this material.

#. P.T. Rakich, C. Reinke, R. Camacho, P. Davids, and Z. Wang, "Giant Enhancement of Stimulated Brillouin Scattering in the Subwavelength Limit,"
   *Phys. Rev. X* **2**, 011008 (2012).  https://doi.org/10.1103/PhysRevX.2.011008

#. C. Wolff, M.J. Steel, B.J. Eggleton, and C.G. Poulton "Stimulated Brillouin scattering in integrated photonic waveguides: Forces, scattering mechanisms, and coupled-mode analysis," *Phys. Rev. A* **92**, 013836 (2015). https://doi.org/10.1103/PhysRevA.92.013836

#. J.E. Sipe and M.J. Steel, "A Hamiltonian treatment of stimulated Brillouin scattering in nanoscale integrated waveguides," *New J. Phys.* **18**, 045004 (2016).  https://doi.org/10.1088/1367-2630/18/4/045004

#. B.C.P Sturmberg at al., "Finite element analysis of stimulated Brillouin scattering in integrated photonic waveguides", *J. Lightwave Technol.*  **37**, 3791-3804 (2019). https://dx.doi.org/10.1109/JLT.2019.2920844

#. C. Wolff, M.J.A. Smith, B. Stiller, and C.G. Poulton, "Brillouin scattering—theory and experiment: tutorial," *J. Opt. Soc. Am. B* **38**, 1243-1269 (2021).  https://doi.org/10.1364/JOSAB.416747


What does |NUMBAT| actually calculate?
=======================================
|NUMBAT| performs three main types of calculations given a particular waveguide design:
 * solve the electromagnetic modal problem  using the finite element method (FEM).
 * solve the elastic modal problem  using FEM.
 * calculate Brillouin gain coefficients and linewidths for a given triplet of two optical and one elastic mode, and use this to generate gain spectra.


Here we specify the precise mathematical problems been solved.
For further details, see the |NUMBAT| paper :cite:p:`Sturmberg:2019` in the Journal of Lightwave Technology  at `<https://dx.doi.org/10.1109/JLT.2019.2920844>`_.

Electromagnetic modal problem
----------------------------------
The electromagnetic wave problem is defined by the vector wave equation

.. math::

   - \nabla \times (\nabla \times {\vec E}) + \omega^2 \epsilon_0 \, \epsilon_r(x,y) \vec E =0,

where the *real-valued* electric field has the following form for modal propagation along :math:`z`:


.. math::
   :nowrap:

   \begin{align*}
   \vec E(x,y,z,t) = &
   \left(
   {\vec {\mathcal{E}}}   (\vecr) e^{- i  \omega t } +
   {\vec {\mathcal{E}}}^* (\vecr) e^{ i  \omega t }
   \right) \\
      = &  \left(
            a(z) \vece(x,y) e^{i (kz-\omega t) } + a^*(z) \vece^*(x,y) e^{-i (kz-\omega t) }
            \right),
   \end{align*}


in terms of the complex field amplitude :math:`\vcalE(\vecr)`, the mode profile :math:`\vece(x,y)` and the complex slowly-varying envelope function :math:`a(z)`. Note that some authors include a factor of :math:`\frac{1}{2}` in these definitions which leads to slightly different expressions for energy and power flow below.

By Faraday's law the complex magnetic field amplitude is given by

.. math::

   {\vec {\cal B} } =  \frac{1}{i \omega} \nabla \times {\vec {\cal E}}.



Elastic  modal problem
----------------------------------

The elastic modal problem is defined by the wave equation

.. math::

   \nabla \cdot \bar{T} + \Omega^2 \rho(x,y) \vec U = 0,

where :math:`\vec u` is the elastic displacement and :math:`\bar{T}=\mathbf{c}(x,y) \bar{S}` is the stress tensor,
defined in terms of the stiffness :math:`\mathbf{c}` and the strain tensor
:math:`\bar{S}=S_{ij} = \frac{1}{2}(\frac{\partial U_i}{\partial r_j} + \frac{\partial U_j}{\partial r_i})`.

The displacement has the modal propagation form

.. math::

   \vec U = b(z) \, \vec u(x,y) e^{i (qz-\Omega t) } + b^*(z) \, \vec u(x,y)^* e^{-i (qz-\Omega t) } ,


For details on how these problems are framed as finite element problems, we refer to Ref. :cite:p:`Sturmberg:2019`.

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
   P^{(o)}(z)  & = \sum_n |a_n(z)|^2 {\cal P}_n^{(o)} \\
   P^{(a)}(z)  & = \sum_n |b_n(z)|^2 {\cal P}_n^{(a)}.

Note that in this convention, the amplitude functions :math:`a_n(z)` and :math:`b_n(z)`
are dimensionless and the dimensionality of the fields lives in the modal functions
:math:`\vec e_n, \vec h_n, \vec U_n`, which are taken to have their conventional SI units.


SBS gain calculation  modal problem
-----------------------------------------

The photoelastic and moving boundary couplings in J/m are given by

.. math::

   Q^{\mathrm{(PE)}} & = - \epsilon \int_A  \mathrm{d}^2 r \, \sum_{ijkl} \epsilon_r^2  \,
   e_i^{(s)*}  \, e_j^{(p)}  \, p_{ijkl}   \, \partial_k u_l^* \\
   Q^{\mathrm{(MB)}} & = \int_{\cal C}  \mathrm{d} {\vec r} \, (\vec u^* \cdot \hat{n})
   \times \\
   & ~~~~
   \left [
   (\epsilon_a - \epsilon_b) \epsilon_0 (\hat{n} \times \vec u^{(s)})^* \cdot (\hat{n} \times \vec e^{(p)})
   -
   (\epsilon_a^{-1} - \epsilon_b^{-1}) \epsilon_0^{-1} (\hat{n} \cdot \vec d^{(s)})^*
          \cdot (\hat{n} \cdot \vec d^{(p)})
   \right]

Note that in general these functions are complex, rather than purely real or imaginary. The equations of motion in the next section show how this is consistent with energy conservation requirements.

Then, at least for backward SBS, the peak SBS gain of the Stokes wave :math:`\Gamma` is given by

.. math::
   \Gamma = \frac{2\omega \Omega}{\alpha_t} \frac{|Q_\mathrm{tot}|^2}{P^{(s)}P^{(p)}{\cal E}^{(a)}},

where the total SBS coupling is :math:`Q_\mathrm{tot} = Q^{(\mathrm{PE})} + Q^{(\mathrm{MB})}`.

Here :math:`\alpha_t` is the temporal elastic loss coefficient in :math:`\mathrm{s}^{-1}`.
It is related to the spatial attenuation coefficient by
:math:`\alpha_s = \alpha_t /v_{\mathrm{p}}^{(\mathrm{a})}` with :math:`v_{\mathrm{p}}^{(\mathrm{a})}` being the elastic phase velocity.

In a backward SBS problem, where there is genuine gain in Stokes optical field propagating in the negative
:math:`z` direction, its optical power evolves as

.. math::

   P^{\mathrm{(s)}}(z) = P_\mathrm{in}^{\mathrm{(s)}} e^{-\Gamma z}.

In forward Brillouin scattering, the same equations for the couplings :math:`Q_i` apply, but it is generally
more helpful to think in terms of the spectral processing of the Stokes field rather than a
simple gain. For this reason, we prefer the general term "forward Brillouin scattering" to
"forward SBS", which you may also encounter.  See refs. :cite:p:`Wolff:2017` :cite:p:`Wolff:2022` for detailed discussion of this issue.

SBS equations of motion
-----------------------------------------

With the above conventions, the dynamical equations for the slowly-varying amplitudes are

..  math::

    \frac{1}{v_g^{p}} \frac{\partial }{\partial t} a_p +   \frac{\partial }{\partial z} a_p & = -i \frac{\omega_p Q_{\mathrm{tot}}}{{\cal P}_o} a_s b     \\
    \frac{1}{v_g^{s}} \frac{\partial }{\partial t} a_s -    \frac{\partial }{\partial z} a_s & = i  \frac{\omega_s Q_{\mathrm{tot}}^*}{{\cal P}_o} a_p b^*   \\
    \frac{1}{v_g^{a}} \frac{\partial }{\partial t} b  +    \frac{\partial }{\partial z} b
    + \frac{\alpha_s}{2} b & =  -i  \frac{\Omega Q_{\mathrm{tot} }^*}{{\cal P}_a}  a_p a_s^*



Here we've chosen the group velocities to be positive and included the propagation direction explicitly.
Note that the coupling :math:`Q_\text{tot}` enters both in native and complex conjugate form. This ensures the total energy behaves appropriately.

Finite element formulation
-------------------------------
|NUMBAT| solves both the electromagnetic and elastic modal properties using finite element method (FEM) algorithms.
For details see refs. :cite:p:`Sturmberg:2019` and :cite:p:`hladkyhennion:1997`.

Here we give a brief outline of the essentials.


Electromagnetic problem
^^^^^^^^^^^^^^^^^^^^^^^^^

Elastic problem
^^^^^^^^^^^^^^^
To formulate a FEM algorithm, the elastic eigenvalue problem is expressed in the so-called *weak form*.
Over the elastic domain :math:`A`, we seek eigenpairs :math:`(\vecu_n(\vecr), \Omega_n)` such that for all test functions :math:`\vecv(\vecr)`, we have

..  math::

    \int_A \vecv^*(\vecr) \cdot \left( \nabla \cdot \bar{T} + \Omega^2 \rho \vecu(\vecr)\right) \, \dA = 0.

The functions :math:`\vecu(\vecr)` and :math:`\vecv(\vecr)` are expanded in a finite set of :math:`N` basis functions :math:`\vecg_m`:

 ..  math::

     \vecu = \sum_{m=1}^N u_m \vecg_m(\vecr),

for some set of coefficients :math:`\vecu_h = (u_1,u_2,\dots u_N)`. In |NUMBAT|, the :math:`\vecg_m` are chosen as piecewise quadratic polynomials on a triangular grid.

As shown in :cite:p:`Sturmberg:2019`, this problem statement ultimately leads to the linear generalised eigenproblem

  .. math::

     \mathrm{K} \vecu_h = \Omega^2 \mathrm{M} \vecu_h ,

where the *stiffness* and *mass* matrices are defined as

.. math::

   K_{lm} = & \int_A (\nabla_s \vecg_l)^* : (\bar{c} : \nabla_s \vecg_m) \, \dA

   M_{lm} = & \int_A \rho  \vecg_l^* \cdot  \vecg_m \, \dA

This problem is then solved using standard linear algebra numerical libraries including ARPACK-NG, UMFPACK on top of the LAPACK and BLAS libraries .

Connecting these to output quantities from code
-------------------------------------------------


Equivalent forms of equations
-----------------------------------

TODO: show forms without the normalisation energies and with cubic style effective area.

Compare to some fiber literature and the hydrodynamic reprn.

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

