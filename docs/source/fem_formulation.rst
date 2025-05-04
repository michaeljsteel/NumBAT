.. include:: numbatdefs.txt

.. _chap-theory-label:

*******************************************************
Numerical formulation of the finite element problems
*******************************************************

Introduction
=============

Here we provide more detail on how the finite element problems defined in chapter :ref:`chap-theory-label` is expressed as numerical problems that can be solved by standard matrix libraries.



We start with the elastic problem which is somewhat easier to solve since it is formulated with a simpler of elements.


Elastic FEM problem
======================

Recall from :ref:`chap-theory-label` that the elastic finite element problem can be expressed in weak form as

.. math::

    \int_A \vecv^*(\vecr) \cdot \left( \nabla \cdot \bar{T} + \Omega^2 \rho \vecu(\vecr)\right) \, \dA = 0.

where we seek the solution eigenpairs  :math:`(\vecu_n(\vecr), \Omega_n)` that satisfy the equation for all test functions :math:`\vecv(\vecr)`.