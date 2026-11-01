.. include:: numbatdefs.txt

.. _chap-migration-label:

*******************************************************
Appendix: Migrating scripts to the ``numbat`` package
*******************************************************

Overview
========

Earlier versions of |NUMBAT| were used by adding the ``backend`` directory to
``sys.path`` and then importing a collection of flat modules:

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

**Your existing scripts continue to work.** The old flat imports remain as
compatibility aliases, and each refers to exactly the same objects as the new
names. Importing a legacy flat module or the ``plotting``, ``meshing``, or
``msh.builtin_meshes`` compatibility packages now emits a ``FutureWarning``
with a suggested public ``numbat`` import. The warning does not mean the
import has stopped working; migrate your scripts when convenient. No removal
release is scheduled by this warning alone.

What does *not* change:

- You still add the ``backend`` directory to ``sys.path`` (or set ``PYTHONPATH``).
- All function, class and method names, and all calculation results, are unchanged.
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
``import integration``                             ``from numbat.solvers import integration``
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
