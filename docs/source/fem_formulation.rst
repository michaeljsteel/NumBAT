.. include:: numbatdefs.txt

.. _chap-fem-formulation-label:

*******************************************************
Numerical formulation of the finite element problems
*******************************************************

Introduction
=============

Here we provide more detail on how the finite element problems defined in chapter :ref:`chap-theory-label` is expressed as numerical problems that can be solved by standard matrix libraries.



We start with the elastic problem which is somewhat easier to solve since it is formulated with a simpler of elements.


Elastic modal problem
======================
Recall from :ref:`chap-theory-label` that the elastic wave equation has the form

.. math::

   \nabla \cdot \bar{T} + \Omega^2 \rho(x,y) \vec u = 0,

where we seek eigensolutions :math:`(\vecu_n(\vecr), \Omega_n)`. See :ref:`chap-theory-label` for the complete definition of the density :math:`\rho(x,y)`, the rank-2 stress tensor :math:`\bar{T}` and its connection with the stiffness tensor :math:`\mathbf{c}`.

In the so-called *weak-formulation* of this eigenproblem, we seek pairs :math:`(\vecu_n(\vecr), \Omega_n)`  so that for all "test" functions :math:`\vecF(\vecr)`, the following integral equation holds

.. math::

    \int_A \vecv^*(\vecr) \cdot \left( \nabla \cdot \bar{T} + \Omega^2 \rho \vecu(\vecr)\right) \, \dA = 0.




Building the FEM formulation
------------------------------

The weak form equation is a general mathematical statement.

To develop a finite-element numerical algorithm,
the functions :math:`\vecu(\vecr)` and :math:`\vecv(\vecr)` are expanded in a finite set of :math:`N` basis functions :math:`g_m(\vecr_t)`:

 ..  math::

     \vecu(\vecr_t) = \sum_{m=1}^N (u^x_m, u^y_m, u^z_m) \vecg_m(\vecr),

for some set of coefficients :math:`\vecu_h = (u^x_1,u^y_1,u^z_1,u^x_2,u^y_2,u^z_2,\dots u_N)`. In |NUMBAT|, the :math:`g_m(\vecr_t)` are chosen as piecewise quadratic polynomials defined on each *element* of an irregularly shaped triangular grid. These basis functions are sometimes known as Lagrange P2 polynomials.

The numerical task is to find eigenvector solutions for the vector of coefficients or *degrees of freedom* :math:`\vecu_h`.

On inserting the field expansion into the weak form integral equation, the eigen-problem ultimately leads to the linear generalised linear eigenvalue equation (see :cite:p:`Sturmberg:2019`),

  .. math::

     \mathrm{K} \vecu_h = \Omega^2 \mathrm{M} \vecu_h ,

where the *stiffness* and *mass* matrices are defined as

.. math::

   K_{lm} = & \int_A (\nabla_s \vecg_l)^* : (\bar{c} : \nabla_s \vecg_m) \, \dA

   M_{lm} = & \int_A \rho  \vecg_l^* \cdot  \vecg_m \, \dA

This problem is then solved using standard linear algebra numerical libraries including ARPACK-NG, UMFPACK on top of the LAPACK and BLAS libraries.

Once the coefficients :math:`\vecu_h` are known for a particular mode, they can be interpolated onto a rectangular grid for plotting. Note that evaluation of the key SBS coupling integrals is performed on the original FEM grid for accuracy.





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

