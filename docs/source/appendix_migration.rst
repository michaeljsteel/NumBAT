.. include:: numbatdefs.txt

.. _chap-migration-label:

*******************************************************
Appendix: Migrating scripts to the ``numbat`` package
*******************************************************

Overview
========

As of |NUMBAT| version 2.4.0, the internal code has been reorganized into a structured package
``numbat_core``. The public interface is now provided by the ``numbat`` package,
which exposes the user-facing pieces under clear names. The tutorials and examples in ``examples/`` now use these new imports. A frozen copy of the original scripts is kept in
``examples_legacy/`` for reference.

Legacy flat imports are no longer supported and raise ``ImportError``. Scripts
must use the public ``numbat`` package imports. A conversion tool is provided
to help rewrite legacy imports in Python scripts and Jupyter notebooks. See the
discussion below for details.

Details
========

EM wavelength units
-------------------

``wl_nm`` is deprecated in EM mode calculations. Use ``wl_um`` in
micrometres instead; for example, replace ``wl_nm=1550`` with
``wl_um=1.55``. Existing nanometre calls continue to work and emit a
``DeprecationWarning``. If both are supplied, ``wl_nm`` takes precedence.

Structure coordinate shifts
---------------------------

``set_xyshift_em`` and ``set_xyshift_ac`` have been removed. Specify a single
``translation=(x, y)`` in nanometres when creating the structure instead.
This moves the mesh and wireframe together, with positive values moving in
the positive coordinate direction. EM and AC calculations share this origin;
separate old plotting offsets cannot in general be converted automatically.
Review the original geometry and cut coordinates when migrating a script.

The independent ``set_r0_offset`` and width-reference APIs remain unchanged.


.. _field-selection-migration-label:

Field selection migration
-------------------------

The legacy ``FieldCode`` and ``FieldTag`` classes remain available during the
2.6 compatibility period. They emit ``DeprecationWarning`` when constructed or
created through a legacy factory. Replace them with immutable typed selections:

========================================================  ======================================================================
Legacy API                                                Replacement
========================================================  ======================================================================
``FieldCode(field_type, force_AC)``                       ``resolve_field_type(field_type, simulation_is_acoustic=...)``
``FieldTag.make_from_field_and_component(ft, "x")``       ``FieldSelection.primary(ft, FieldComponent.X)``
``tag.set_to_minor()``                                    ``selection.with_part(ComplexPart.SECONDARY)``
``tag.component_as_F()``                                  ``storage_key(selection)``
``tag.get_tex_plot_label()``                              ``tex_label(selection)``
``tag.linecolor()``                                       ``line_color(selection)``
``tag.linestyle(comps)``                                  ``line_style(selection, show_magnitude=...)``
========================================================  ======================================================================

For plotting components, prefer ``FieldComponent`` values such as
``FieldComponent.X``. Existing ``"x"``, ``"y"``, ``"z"``, ``"t"``, and
``"a"`` string inputs remain supported during the compatibility period.

Prior to version 2.4.0, |NUMBAT| was used by adding the ``backend`` directory to ``sys.path`` and then importing a collection of flat modules:

.. code-block:: python

    sys.path.append(str(Path('../../backend')))
    import numbat
    import materials
    import integration
    from nbtypes import SI_GHz

The code now lives in a structured internal package ``numbat_core``, and
``numbat`` is a public package that exposes the user-facing pieces under clear
names. The tutorials and examples in ``examples/`` now use these new imports.
A frozen copy of the original scripts is kept in ``examples_legacy/`` for
reference.

**Update existing scripts before running them.** Legacy flat modules and the
``plotting``, ``meshing``, and ``msh.builtin_meshes`` compatibility packages
raise ``ImportError`` with the corresponding public ``numbat`` import path.

What does *not* change:

- You still add the ``backend`` directory to ``sys.path`` (or set ``PYTHONPATH``).
- Apart from the migrations described above, function, class and method names,
  and calculation results, are unchanged.
- ``numbat.NumBATApp``, ``numbat.load_simulation`` and every other name that
  used to be available as ``numbat.<name>`` still is.
- Saved simulation results and meshes load as before.

Import mapping
==============

=================================================  =========================================================
Old import                                         New import
=================================================  =========================================================
``import numbat``                                  ``import numbat`` (unchanged)
``import materials``                               ``from numbat import materials``
``import integration``                             ``from numbat.solvers import sbsgain as integration``
``import modecalcs``                               ``from numbat.solvers import modecalcs``
``import modes``                                   ``from numbat.solvers import modes``
``import simresult``                               ``from numbat.solvers import simresult``
``import structure``                               ``from numbat.geometry import structure``
``from structure import Structure``                ``from numbat.geometry import Structure``
``from usermesh import UserGeometryBase``          ``from numbat.mesh import UserGeometryBase``
``import nbgmsh``                                  ``from numbat.mesh import MailData``
                                                   (and use ``MailData(...)`` in place of ``nbgmsh.MailData(...)``)
``from nbtypes import SI_GHz``                     ``from numbat.types import SI_GHz``
``from nbthread import launch_worker_...``         ``from numbat.parallel import launch_worker_...``
``from numbattools import ...``                    ``from numbat.utils import ...``
``import plotgain``                                ``from numbat import plotting``
                                                   (and use ``plotting.<name>`` in place of ``plotgain.<name>``)
``from plotgain import Decorator``                 ``from numbat.plotting import Decorator``
``from plotting.ipythontools import ...``          ``from numbat.plotting.ipython import ...``
``import plotting``                                ``from numbat import plotting``
``import meshing``                                 ``from numbat import mesh``
``import msh.builtin_meshes``                      ``from numbat import mesh``
=================================================  =========================================================

The new package also offers ``numbat.mesh.FemMesh`` and the simulation and
result classes in ``numbat.solvers`` (``EMSimulation``, ``ACSimulation``,
``Mode``, ``EMSimResult``, ``ACSimResult``).

Notes for Jupyter notebooks
---------------------------

Notebooks use the same imports. If your notebook imported
``mode_calcs`` or the ``plotting`` package directly, replace these with
``numbat.solvers.modecalcs`` and ``numbat.plotting``.

Multiprocessing on Python 3.14
==============================

This is unrelated to the import changes but affects the multicore examples.
Since Python 3.14, Linux uses the ``forkserver`` start method for
``multiprocessing`` by default. Scripts that create a ``Pool`` at module level
and rely on global variables in worker functions (as the tutorials do) must
request the ``fork`` context explicitly, and should select a non-interactive
matplotlib backend before importing ``pyplot``:

.. code-block:: python

    import multiprocessing
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ...
    with multiprocessing.get_context("fork").Pool(num_cores) as pool:
        ...

Automatic conversion
====================

The script ``scripts/nb_convert_imports.py`` rewrites legacy imports in Python
scripts and Jupyter notebooks. Only ``import`` statements are changed.

Preview the changes to a file or directory (nothing is written):

.. code-block:: bash

    python scripts/nb_convert_imports.py my_project/

Apply them, keeping the original of each changed file as ``<name>.orig``:

.. code-block:: bash

    python scripts/nb_convert_imports.py --write my_project/

Use ``--no-backup`` to skip the ``.orig`` copies. The tool is safe to run
repeatedly. Lines it cannot convert safely, such as ``import numbat, materials``,
are reported for you to edit by hand.

Two conversions keep your code working without further edits:

- ``import plotgain`` becomes ``from numbat import plotting as plotgain``.
- ``import nbgmsh`` becomes ``from numbat import mesh as nbgmsh``.

You may then tidy these up by hand to the forms shown in the table above.
