
.. include:: numbatdefs.txt

.. _chap-tutorial-label:

***********
Tutorial
***********

This chapter provides a sequence of graded tutorials for learning |NUMBAT|,
exploring its applications and validating it against literature results and
analytic solutions where possible. Before attempting your own calculations with
|NUMBAT|, we strongly advise working through the sequence of tutorial exercises
which are largely based on literature results.

We will meet a significant number of |NUMBAT| functions in these tutorials, though certainly not all.  The full Python interface is documented in the section :ref:`chap-pythonapi-label`.

You may then choose to explore relevant examples drawn from a recent
tutorial paper by Dr Mike Smith and colleagues, and a range of other literature
studies,  which are provided in the following two chapters,
:ref:`chap-josab-label` and :ref:`chap-literature-label`.



Some Key Symbols
=====================
As far as practical we use consistent notation and symbols in the tutorial
files. The following list introduces a few commonly encountered ones. Note that
with the exception of the free-space wavelength :math:`\lambda` and the spatial
dimensions of waveguide structures, which are both specified in nanometres
(nm), all quantities in |NUMBAT| should be expressed in the standard SI units.
For example, elastic frequencies :math:`\nu` are expressed in Hz, not GHz.

``lambda_nm``
    This is the *free-space* optical wavelength :math:`\lambda` satisfying
    :math:`\lambda = 2\pi c/\omega`, where :math:`c` is the speed of light and :math:`\omega` is the
    angular frequency. **For convenience, this parameter is specified in nm.**

    For most examples, we use the conventional value :math:`\lambda` =1550 nm.

``omega, omega_EM, om_EM``
    This is the electromagnetic *angular* frequency :math:`\omega = 2 \pi c/\lambda`  specified in
    :math:`\mathrm{rad.s}^{-1}`.

``k, beta, k_EM``
    This is the electromagnetic *wavenumber* or *propagation constant* :math:`k`
    or :math:`\beta`, specified in :math:`\mathrm{m}^{-1}`.

``neff, n_eff``
    This is the electromagnetic modal *effective index* :math:`\bar{n}=ck/\omega`, which is dimensionless.

``nu, nu_AC``
    This is the acoustic frequency :math:`\nu` specified in Hz.

``Omega, Omega_AC, Om_AC``
    This is the acoustic *angular* frequency :math:`\Omega = 2 \pi \nu`  specified in :math:`\mathrm{rad.s}^{-1}`.

``q, q_AC``
    This is the acoustic *wavenumber* or *propagation constant* :math:`q=v_{ac} \Omega`,
    where :math:`v_{ac}` is the phase speed of the wave.
    The acoustic wavenumber is specified in :math:`\mathrm{m}^{-1}`.

``m``
    This is an integer(8) corresponding to the mode number :math:`m` of an electromagnetic
    mode :math:`\vec E_m(\vec r)` or an acoustic mode :math:`\vec u_m(\vec r)`.

    For both electromagnetic and acoustic modes, counting of modes begins with ``m=0``
    and are ordered by decreasing effective index and increasing frequency respectively.

    For the electromagnetic problem in which frequency/free-space wavelength is the
    independent variable, the :math:`m=0` mode has the *highest* effective index
    :math:`\bar{n}` and *highest* wavenumber :math:`k` of any mode for a given angular frequency
    :math:`\omega`.

    For the acoustic problem, the wavenumber :math:`q` is the
    independent variable and we solve for frequency :math:`\nu=\Omega/(2\pi)`.
    The :math:`m=0` mode has the *lowest* frequency
    :math:`\nu` of any mode for a given wavenumber :math:`q`.

    The integer(8) :math:`m` therefore has no particular correspondence to
    the conventional two index mode indices for fibre or rectangular waveguides.


``inc_a_x, inc_a_y, inc_b_x, inc_b_y, slab_a_x, slab_a_y,``  ...  etc
   These are dimensional parameters specifying the lengths of different aspects
   of a given structure: rib height, fibre radius etc.
   **For convenience, these parameters are specified in nm.**


.. raw:: latex

    \clearpage





Elementary Tutorials
=======================

We now walk through a number of simple simulations that demonstrate the basic use of |NUMBAT| located in the ``<NumBAT>/tutorials`` directory.


.. toctree::
   :maxdepth: 1

   tute_elem

   notebooks/jup_09_smf28


.. raw:: latex

    \clearpage


Intermediate tutorials
======================

The next set of tutorials begin to explore some more advanced features.

.. toctree::
   :maxdepth: 1


   notebooks/jup_09a_bulk_anisotropy

.. raw:: latex

    \clearpage

.. toctree::
   :maxdepth: 1

   tute_intermediate



