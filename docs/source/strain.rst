################################
Acoustic strain post-processing
################################

Acoustic strain is evaluated directly from the native quadratic FEM
displacement solution. It is not calculated by differentiating an
interpolated rectangular displacement grid. The six independent components
are ``Sxx``, ``Syy``, ``Szz``, ``Syz``, ``Sxz``, and ``Sxy``.

2D plots and exports
====================

Use the acoustic simulation result to plot or export selected modes:

.. code-block:: python

   simres_AC.plot_strain_modes(mode_indices=[0])
   simres_AC.write_strain_modes(mode_indices=[0])

Omit ``mode_indices`` (or pass ``"All"``) to process every acoustic mode.

With no ``comps`` argument, ``plot_strain_modes`` writes a six-panel
overview using the display-primary components ``Re[Sxx]``, ``Re[Syy]``,
``Re[Szz]``, ``Im[Syz]``, ``Im[Sxz]``, and ``Re[Sxy]``. Request individual
components with names such as ``comps="xx,yz"``. Plot values and colourbars
are displayed in megastrain (``Mε``), meaning the native dimensionless strain
is scaled by :math:`10^{-6}`. Text exports retain the unscaled values.

Each 2D export produces twelve files:
``AC_strain_mode_{mode:02d}_S{ij}_{re|im}.txt``. They contain the physical
real or imaginary sampled component grid, with coordinate bounds in the file
header.

1D cuts and exports
===================

Native strain cuts use the same ``"x"``, ``"y"``, and ``"line"`` conventions
as displacement cuts:

.. code-block:: python

   simres_AC.plot_strain_modes_1D("x", 0, mode_indices=[0], comps="xx,yz")
   simres_AC.write_strain_modes_1D("x", 0, mode_indices=[0])
   simres_AC.plot_strain_modes_1D(
       "line", (-0.5, 0), (0.5, 0), mode_indices=[0]
   )

An exported cut is named ``AC_strain_mode_{mode:02d}_{cut}cut.txt`` and has
the columns ``x Sxxr Sxxi Syyr Syyi Szzr Szzi Syzr Syzi Sxzr Sxzi Sxyr
Sxyi``. At material interfaces, cuts retain a deterministic one-sided value:
the first matching native local subtriangle is used and values are never
averaged across the interface.
