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





Connecting these to output quantities from code
-------------------------------------------------


Equivalent forms of equations
-----------------------------------

TODO: show forms without the normalisation energies and with cubic style effective area.

Compare to some fiber literature and the hydrodynamic reprn.




The finite element method in |NUMBAT|
=========================================

|NUMBAT| solves both the electromagnetic and elastic modal properties using finite element method (FEM) algorithms.
For details see refs. :cite:p:`Sturmberg:2019`, :cite:p:`Dossou:2005` and :cite:p:`hladkyhennion:1997`.

Here we give a brief outline of the essentials.


Electromagnetic problem
-----------------------------
We closely follow the exposition in :cite:p:`Dossou:2005`.

Expressed in the modal form :math:`\vecE(\vecr)= [\vecE_t, E_z]e^{i \beta z}`, the wave equation becomes the pair of equations:

.. math::

   \nabla_t \times \left(\frac{1}{\mu}(\nabla_t \times \vecE_t) \right)
      - \omega^2 \epsilon \vecE_t = & \, \beta^2 \frac{1}{\mu} ( \nabla \hE_z- \vecE_t)

   - \nabla_t \cdot (\frac{1}{\mu} \vecE_t) + \nabla_t \cdot(\frac{1}{\mu} \nabla_t \hE_z) + \omega^2 \epsilon \hE_z =& \, 0 ,

where for convenience we take :math:`E_z = -\beta \hE_z`.

In the so-called *weak-formulation*, we seek pairs :math:`(\vecE, \beta)` so that for all test functions :math:`\vecF(\vecr)`, the following equations hold

.. math::

   \left(\frac{1}{\mu} (\nabla_t \times \vecE_t), \nabla_t \times \vecF_t\right)
      - \omega^2 \left(\epsilon \vecE_t, \vecF_t\right)
           =  \, \beta^2 \left(\frac{1}{\mu} ( \nabla \hE_z- \vecE_t), \vecF_t \right)

   \left(\frac{1}{\mu} \vecE_t,\nabla_t F_z  \right) -
      \left(\frac{1}{\mu} \nabla_t \hE_z, \nabla_t F_z \right)
      + \omega^2 \left (\epsilon \hE_z, F_z \right ) = \, 0 ,

To build the FEM formulation, we introduce transverse vector basis functions :math:`\vphi_i(x,y)` and scalar longitudinal functions :math:`\psi_i(x,y)`  so that a general field has the form

..  math::

     \vecE = & \vecE_t + \unitz E_z \\
         = & \sum_{i=1}^{N} e_{t,i} \, \vphi_i(\vecr)
                + \unitz \sum_{i=1}^{N} \he_{z,i} \, \psi_i(\vecr) \\
                = & [\vphi_1, \vphi_2, \ldots, \vphi_N ,\psi_1, \psi_2, \ldots, \psi_N, ] \begin{bmatrix} e_{t,1} \\ e_{t,2}   \\ \ldots \\ e_{t,N} \\
                \he_{z,1} \\ \he_{z,2} \\ \ldots \\ \he_{z,N} \end{bmatrix} \\
                = & [\bvphi ; \unitz \bpsi]
                 \begin{bmatrix} \bfe_t \\ \bfhez \end{bmatrix}

An arbitrary inner product :math:`(\calL_1 E, \calL_2 F)` with :math:`\calF = a \vphi_m+ b\psi_m \unitz` where :math:`a` and :math:`b` are zero or one, expands  to

..  math::

     (\calL_1 \vecE, \calL_2 \vecF) = & \int (\calL_2 \vecF)^* \cdot (\calL_1 \vecE) \, \dA \\
      = & \int (\calL_2 a \vphi_m + b \psi_m \unitz)^* \cdot
      ( \calL_1
         [\bvphi ; \unitz \bpsi]
                 \begin{bmatrix} \bfe_t \\ \mathbf{\hat{e}}_z \end{bmatrix} ) \\
      = & a \int (\calL_2 \vphi_m )^* \cdot ( \calL_1 \vphi_n ) \, \dA \, e_{t,n}
        + b \int (\calL_2 \psi_m \unitz )^* \cdot ( \calL_1 \vphi_n ) \, \dA\,  e_{t,n} \\
       ~ & + a\int  (\calL_2 \vphi_m )^* \cdot ( \calL_1 \psi_n \unitz ) \, \dA\,\he_{z,n}
        + b \int (\calL_2 \psi_m \unitz )^* \cdot ( \calL_1 \psi_n \unitz ) \, \dA\, \he_{z,n}
      \\
      = & \calL_{tt} \bfe_{t} +  \calL_{zt} \bfe_{z}
           + \calL_{tz}\bfhez +  \calL_{zz} \bfhez.


With these definitions we can identify the matrices

..  math::

    K_{mn} = & \int \left [ (\nabla_t \times \vphi_m )^*  \cdot   \frac{1}{\mu} (\nabla_t \times \vphi_n) - \omega^2 \vphi_m ^* \cdot ( \epsilon \vphi_n) \right ] \, \dA  \\
    K_{zt} = & \int (\nabla_t \psi_m )^*  \cdot  \frac{1}{\mu}   \vphi_n \, \dA  \\
    K_{zz} = & \omega^2 \int   \psi_m ^*  \cdot  \epsilon   \psi_n \, \dA  \\
    M_{tt} = & - \int   \vphi_m ^*  \cdot    \vphi_n \, \dA  \\


Then the eigenproblem to be solved is the generalised linear  problem

..  math::

    \begin{bmatrix} [K_{tt}] &[ 0] \\ [K_{zt}] & [K_{zz}] \end{bmatrix}
    \begin{bmatrix} e_{t,h} \\ e_{z,h} \end{bmatrix}
    = \beta^2 \begin{bmatrix} [M_{tt}] &[K_{zt}]^T \\ [0]  & [0] \end{bmatrix}
      \begin{bmatrix} e_{t,h} \\ e_{z,h} \end{bmatrix} .

By swapping the sides of the second row, the two matrices involved become symmetric:

..  math::

    \begin{bmatrix} [K_{tt}] &[ 0] \\ [0] & [0] \end{bmatrix}
    \begin{bmatrix} e_{t,h} \\ e_{z,h} \end{bmatrix}
    = \beta^2 \begin{bmatrix} [M_{tt}] &[K_{zt}]^T \\ [K_{zt}]  & [K_{zz}] \end{bmatrix}
      \begin{bmatrix} e_{t,h} \\ e_{z,h} \end{bmatrix},

which is ideally posed to solve using standard numerical libraries.


Elastic problem
--------------------
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

