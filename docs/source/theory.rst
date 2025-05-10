.. include:: numbatdefs.txt

.. _chap-theory-label:

****************************
Background theory
****************************






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

where :math:`\omega=2\pi/\lambda` is the optical angular frequency with :math:`\lambda` the free space wavelength. The waveguide properties are defined by the spatially varying relative dielectric constant :math:`\epsilon_r(x,y)`.

The *real-valued* electric field :math:`\vec E(x,y,z,t)` has the following form for modal propagation along :math:`z`:

.. math::
   :nowrap:

   \begin{align*}
   \vec E(x,y,z,t) = &
   \left(
   {\vec {\mathcal{E}}}   (\vecr) e^{- i  \omega t } +
   {\vec {\mathcal{E}}}^* (\vecr) e^{ i  \omega t }
   \right) \\
      = &  \left(
            a \vece(x,y) e^{i (kz-\omega t) } + a^* \vece^*(x,y) e^{-i (kz-\omega t) }
            \right),
   \end{align*}

in terms of the wavenumber :math:`k`, complex field amplitude :math:`\vcalE(\vecr)`, mode profile :math:`\vece(x,y)` and the mode amplitude :math:`a`. Note that some authors include a factor of :math:`\frac{1}{2}` in these definitions which leads to slightly different expressions for the energy and power flow below.

By Faraday's law the complex magnetic field amplitude is given by

.. math::

   {\vec {\cal B} } =  \frac{1}{i \omega} \nabla \times {\vec {\cal E}}.

Dependent and independent variables in the modal dispersion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In |NUMBAT|, the electromagnetic mode problem is formulated with the optical angular frequency :math:`\omega` as the independent variable and the wavenumber or *propagation constant* :math:`k(\omega)` as the eigenvalue.  For a given frequency, |NUMBAT| thus finds the eigenvalues :math:`k_n(\omega)` and eigenmodes :math:`\vece_n(x,y)`. The algorithm finds the most-bound or "lowest" eigenvalues first. Thus the values :math:`k_n(\omega)` *decrease* with :math:`n` and there are only a finite number of guided mode solutions for a given structure and frequency :math:`\omega`. Note that the propagation constant is also often denoted by :math:`\beta_n(\omega)`.

It is thus straightforward to generate dispersion relations :math:`k_n(\omega)` on an equal-spaced frequency grid. If you require the wavenumber on an equal-spaced grid, you will need to interpolate the :math:`k_n(\omega)` functions.

Elastic  modal problem
----------------------------------

The elastic modal problem is defined by the wave equation

.. math::

   \nabla \cdot \bar{T} + \Omega^2 \rho(x,y) \vec U = 0,

where :math:`\Omega` is the elastic angular frequency and `:math:`\vec U` is the elastic displacement.
Further, :math:`\bar{T}=\mathbf{c}(x,y) \bar{S}` is the rank-2 stress tensor,
defined in terms of the rank-4 stiffness tensor :math:`\mathbf{c}` and the rank-2 strain tensor
:math:`\bar{S}=S_{ij} = \frac{1}{2}(\frac{\partial U_i}{\partial r_j} + \frac{\partial U_j}{\partial r_i})`.

The waveguide properties are defined by the spatially varying density and stiffness functions.

The displacement field has the modal propagation form

.. math::

   \vec U = b(z) \, \vec u(x,y) e^{i (qz-\Omega t) } + b^*(z) \, \vec u(x,y)^* e^{-i (qz-\Omega t) } .

where :math:`q` is the elastic wavenumber or propagation constant.

Dependent and independent variables in the modal dispersion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In |NUMBAT|, the elastic mode problem is formulated in the *opposite sense* to the electromagnetic problem.
In this case, the elastic wavenumber :math:`q` is the independent variable, and |NUMBAT| solves for the corresponding angular frequency eigenvalues :math:`\Omega_n(q)` and eigenmodes :math:`\vecu_n(x,y)`. With, the wavenumber as the independent variable, the eigenvalues :math:`q_n(\Omega)` *increase* with the mode index :math:`n`. In such a formulation, there is no upper bound to the number of physical eigenstates, other than the numerical resolution of the FEM grid. In practice, above frequencies of :math:`\Omega/(2\pi)` of a few tens of GHz, the elastic losses become so large that the propagation is restricted to distances of only a few microns.

It is thus straightforward to generate dispersion relations :math:`\Omega_n(q)` on an equal-spaced frequency grid. If you require the wavenumber on an equal-spaced grid, you will need to interpolate the :math:`\Omega_n(q)` functions.



For details on how these problems are framed as finite element problems, we refer to Ref. :cite:p:`Sturmberg:2019`, though some details are provided in :ref:`chap-fem-formulation-label`.

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





Connecting these to output quantities from code
-------------------------------------------------


Equivalent forms of equations
-----------------------------------

TODO: show forms without the normalisation energies and with cubic style effective area.

Compare to some fiber literature and the hydrodynamic reprn.




The method of solution in |NUMBAT|
=========================================

|NUMBAT| solves both the electromagnetic and elastic modal properties using *finite element method* (FEM) algorithms.
Chapter :ref:`chap-fem-formulation-label` provides a summary of the finite element formulation.

For further details, see refs. :cite:p:`Sturmberg:2019`, :cite:p:`Dossou:2005` and :cite:p:`hladkyhennion:1997`.





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

