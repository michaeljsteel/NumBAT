
.. include:: numbatdefs.txt

.. _chap-technical-label:

*******************
Additional details
*******************


User plot preferences
-------------------------------------------
A large number of plotting properties can be controlled using the ``numbat.toml`` file.

A sample file is provided in the root directory of the |NUMBAT| installation.

To provide an adjustable set of global settings copy this into your home directory,
either as ``~/.numbat.toml`` on Linux or MacOS, or ``numbat.toml`` on Windows.

You can also places copies of this file in any working directory to override
your normal defaults for a particular set of calculations.


Specifying and using anistropic materials

-------------------------------------------

WRITE ME

Orientation of the coordinate axes in NumBAT

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Cartesian coordinates in |NUMBAT| are defined so that the waveguide lies
in the :math:`x-y` plane with propagation along :math:`z`.

To obtain a right-handed system, one should think of the propagation as being out
of the screen (though it is rare that this matters) :

  - :math:`x` : Increases to the right. Usually the dominant electric field component for modes labelled as  TE polarisation.
  - :math:`y` : Increases up the screen. Usually the dominant electric field component for modes labelled as  TM polarisation.
  - :math:`z` : Increases out of the screen. Propagation direction.


It is common to encounter different coordinate choices in the literature, in particular

the *bird's-eye view* corresponding to viewing a photonic circuit from above is frequently encountered:

  - :math:`x`: Increases to the right.  Usually the dominant electric field component for modes labelled as TE polarisation.
  - :math:`y`: Increases out of the screen. Propagation direction.
  - :math:`z`: Increases up the screen.  Usually the dominant electric field component for modes labelled as TM polarisation.

When dealing with elastically anisotropic materials, it is important to understand the relationship between the
|NUMBAT| coordinate system and that of any literature source you may be consulting, as it may be necessary to
perform a rotation on the tensor properties of the material in question.  This is discussed below.




Supported crystal classes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Depending on the degree of crystal symmetry of a material, the optical and elastic response tensors
obey constraints that can substantially reduce the number of independent components.

|NUMBAT| currently supports the following crystal classes: Isotropic, Cubic, Trigonal.

Note that this choice merely reflects what has been useful to the authors
and new classes are easily added as needed. Please get in touch.

The number and location of the independent elements of each tensor &mdash;
stiffness  tensor :math:`c_{ijkl}`,

photolelastic  tensor :math:`p_{ijkl}`,

viscosity tensor :math:`p_{ijkl}` &mdash; depends on the particular crystal class.

As all the relevant tensors are symmetric, there can be at most 6 independent elements for rank-2 tensors (rather than 9),
and 36 independent elements in rank-4 tensors (rather than 81).  In particular, tensor subscripts

are all symmetric in pairs as follows:

.. math::

   T_{ij} = T_{ji}
   T_{ijkl} = T_{jikl} = T_{ijlk} = T_{jikl}


where the subscripts range over all values of

:math:`x,y,z`.

Consequently we can use the Voigt notation

in which pairs of indices are represented by a single integer(8) in the range 1 to 6 as follows:


.. math:: \begin{bmatrix}  xx \\ yy \\ zz \\ xz \\ yz \\ zz \end{bmatrix}

   = \begin{bmatrix}  1 \\ 2 \\ 3 \\ 4 \\ 5 \\ 6 \end{bmatrix}


Then a rank two tensor like strain or stress can be represented as the column

.. math:: \bar{S} \equiv S_I \equiv

    \begin{bmatrix}  S_1 \\ S_2 \\ S_3 \\ S_4 \\ S_5 \\ S_6 \end{bmatrix}  .

A fourth rank tensor like the stiffness or photoelastic tensor is represented in the form


.. math:: \bar{c} \equiv c_{IJ} \equiv

    \begin{bmatrix}

   c_{11} & c_{12} & c_{13} & c_{14} & c_{15} & c_{16} \\
   c_{21} & c_{22} & c_{23} & c_{24} & c_{25} & c_{26} \\
   c_{31} & c_{32} & c_{33} & c_{34} & c_{35} & c_{36} \\
   c_{41} & c_{42} & c_{43} & c_{44} & c_{45} & c_{46} \\
   c_{51} & c_{52} & c_{53} & c_{54} & c_{55} & c_{56} \\
   c_{61} & c_{62} & c_{63} & c_{64} & c_{65} & c_{66} \\
   \end{bmatrix}  .



This table summarises the symmetry properties for each class and the required tensor elements that must be specified in a material ``.json`` file.


.. table:: Properties of crystal classes and required elements in |NUMBAT|.

    ==========  =======================================   ==========    ================================   ==========================================================   ==================================================================================
    Structure   Nature                                     Examples      :math:`\epsilon`                   Stiffness elements                                           Photoelastic elements
    ==========  =======================================   ==========    ================================   ==========================================================   ==================================================================================
    Isotropic   Every direction equivalent                 Glass         :math:`\epsilon_1`                         :math:`c_{11}`                                        :math:`p_{11}, p_{12}, p_{14}`

    Cubic       3 equivalent perpendicular directions      Silicon       :math:`\epsilon_1`                         :math:`c_{11}, c_{12}, c_{44}`                        :math:`p_{11}, p_{12}, p_{44}`

    Trigonal    1 threefold rotation axis                  LiNbO3        :math:`\epsilon_1, \epsilon_3`     :math:`c_{11}, c_{12}, c_{13}, c_{14}, c_{33}, c_{44}`        :math:`p_{11}, p_{12}, p_{13}, p_{14}, p_{31}, p_{33}, p_{41}, p_{44}`
    ==========  =======================================   ==========    ================================   ==========================================================   ==================================================================================


The full form of the material tensors  for each crystal class is as follows:


**Isotropic**


.. math::
   \epsilon_{I} =

    \begin{bmatrix}
    \epsilon_1 & 0 & 0 \\
    0 & \epsilon_1 &  0 \\
    0 & 0 & \epsilon_1

    \end{bmatrix}

   \quad
   c_{IJ} =

    \begin{bmatrix}
   c_{11} & c_{12} & c_{12} & 0      & 0      & 0 \\
   c_{12} & c_{11} & c_{12} & 0      & 0      & 0 \\
   c_{12} & c_{12} & c_{11} & 0      & 0      & 0 \\
   0      & 0      & 0      & c_{44} & 0      & 0 \\
   0      & 0      & 0      & 0      & c_{44} & 0 \\
   0      & 0      & 0      & 0      & 0      & c_{44} \\
    \end{bmatrix}

   p_{IJ} =

    \begin{bmatrix}
   p_{11} & p_{12} & p_{12} & 0      & 0      & 0 \\
   p_{12} & p_{11} & p_{12} & 0      & 0      & 0 \\
   p_{12} & p_{12} & p_{11} & 0      & 0      & 0 \\
   0      & 0      & 0      & p_{44} & 0      & 0 \\
   0      & 0      & 0      & 0      & p_{44} & 0 \\
   0      & 0      & 0      & 0      & 0      & p_{44} \\
    \end{bmatrix}

where :math:`c_{44} = (c_{11}-c_{12})/2` and :math:`p_{44} = (p_{11}-p_{12})/2`.
These quantities are related to the *Lame parameters* :math:`\mu=c_{44}` and :math:`\lambda=c_{11}`
which may in turn be expressed in terms of the Young's modulus :math:`E` and Poisson ratio :math:`\nu`:

.. math::
   \mu = \frac{E}{2(1+\nu)}, \qquad \lambda = \frac{E\nu}{(1+\nu)(1-2\nu)}

**Cubic**


The matrix expressions are identical to the isotropic case except that

:math:`c_{44}` and :math:`p_{44}` are now independent quantities that must be specified directly.

**Trigonal**


.. math::
   \epsilon_{I} =

    \begin{bmatrix}
    \epsilon_1 & 0 & 0 \\
    0 & \epsilon_1 &  0 \\
    0 & 0 & \epsilon_3

    \end{bmatrix}

   \quad
   c_{IJ} =

    \begin{bmatrix}
   c_{11} & c_{12} & c_{13} & c_{14} & 0      & 0 \\
   c_{12} & c_{11} & c_{13} & -c_{14}& 0      & 0 \\
   c_{13} & c_{13} & c_{33} & 0      & 0      & 0 \\
   c_{14} &-c_{14} & 0      & c_{44} & 0      & 0 \\
   0      & 0      & 0      & 0      & c_{44} & c_{14} \\
   0      & 0      & 0      & 0      & c_{14} & (c_{11}-c_{12})/2
    \end{bmatrix}

   p_{IJ} =

    \begin{bmatrix}
   p_{11} & p_{12} & p_{13} & p_{14} & 0      & 0 \\
   p_{12} & p_{11} & p_{13} &-p_{14} & 0      & 0 \\
   p_{31} & p_{31} & p_{33} & 0      & 0      & 0 \\
   p_{41} & -p_{41}& 0      & p_{44} & 0      & 0 \\
   0      & 0      & 0      & 0      & p_{44} & p_{41} \\
   0      & 0      & 0      & 0      & p_{14} &  (p_{11}-p_{12})/2
    \end{bmatrix}






WRITE ME

Required tensor components
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

WRITE ME




