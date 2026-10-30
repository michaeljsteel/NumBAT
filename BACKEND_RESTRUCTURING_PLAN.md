# Backend Restructuring Plan

## Goal

Move NumBAT's backend from its current flat import namespace to a structured
internal package hierarchy without breaking existing tutorials, tests, dynamic
template loading, or user scripts.

The migration should be compatibility-first and dependency-driven. Do not
combine it with behaviour changes, broad formatting, or unrelated refactors.

## Current constraints

- Tutorials commonly add `backend` directly to `sys.path` and then use
  `import numbat`.
- Backend modules use flat absolute imports such as `import reporting` and
  `from modes import ModeEM`.
- The application root (`numbat`) is imported from low-level mesh and plotting
  code to access `NumBATApp`, producing dependency cycles.
- Mesh metadata currently identifies built-in implementations by stable module
  paths such as `msh.builtin_meshes`.
- Dynamic geometry/template loading and pickle identity depend on importable,
  stable module names.

## Proposed destination

Introduce `backend/numbat_core` as the internal package during the transition.
Keep the current top-level modules as compatibility facades until a later,
explicit deprecation release.

```text
backend/
├── numbat.py                     # Legacy public entry point / facade
├── materials.py                  # Legacy facade
├── femmesh.py                    # Legacy facade
├── modes.py                      # Legacy facade
├── ...
│
├── numbat_core/
│   ├── __init__.py
│   ├── runtime/
│   │   ├── app.py
│   │   ├── configuration.py
│   │   └── reporting.py
│   ├── geometry/
│   │   ├── structure.py
│   │   ├── usermesh.py
│   │   └── waveguides.py
│   ├── materials/
│   │   ├── materials.py
│   │   ├── materialprops.py
│   │   ├── bulkprops.py
│   │   ├── crystalsyms.py
│   │   └── elastic_tensors.py
│   ├── mesh/
│   │   ├── fem.py
│   │   ├── gmsh_io.py
│   │   ├── mail.py
│   │   ├── templates.py
│   │   └── resources/
│   ├── solvers/
│   │   ├── modecalcs.py
│   │   ├── modes.py
│   │   ├── integration.py
│   │   └── simresult.py
│   ├── plotting/
│   │   ├── fields.py
│   │   ├── modes.py
│   │   ├── profiles.py
│   │   └── mesh.py
│   ├── analytics/
│   │   └── ...
│   └── utils/
│       ├── numeric.py
│       ├── platform.py
│       └── threading.py
│
├── meshing/                      # Compatibility package during transition
├── plotting/                     # Compatibility package during transition
├── msh/                           # Retained until resource migration is complete
└── nbanalytic/                   # Migrate later
```

The final names may change, but the core rule is that the new hierarchy is
internal first. Existing import paths remain valid throughout the migration.

## Desired dependency direction

```text
utils / runtime / materials
            ↓
geometry + mesh
            ↓
solvers
            ↓
plotting
            ↓
public application API
```

Low-level modules must not import the public application module. In particular,
mesh and plotting code should receive output directories or an output context
instead of importing `numbat` only to reach `NumBATApp`.

## Staged transition

### Phase 0: Establish a safety net

Do not move production files yet.

1. Add import-compatibility tests that run with `PYTHONPATH=backend`.
2. Verify imports used by tutorials and tests, including:

   ```python
   import numbat
   import materials
   import femmesh
   import modecalcs
   import modes
   import structure
   import simresult
   import meshing.templates
   import plotting.profiles
   ```

3. Add focused smoke tests for:
   - built-in waveguide construction;
   - mesh generation;
   - EM simulation object construction;
   - built-in template loading;
   - user geometry serialization and deserialization.
4. Record circular imports and import-time dependencies around the application,
   structure, mesh, solver, and plotting modules.

### Phase 1: Create boundaries before core moves

Create `numbat_core.runtime` and move or extract low-level application-runtime
concerns there.

The immediate objective is to remove imports of `numbat` from low-level mesh
and plotting code. Prefer explicit output paths:

```python
plotter.save(output_dir, prefix)
```

over application-root imports:

```python
numbat.NumBATApp().outdir()
```

A small runtime context helper is acceptable as an intermediate step if passing
the paths explicitly would make the first migration excessively invasive.

### Phase 2: Move low-risk leaf modules

Move independent modules first, leaving top-level facades that re-export the
old public API.

Suggested early candidates:

| Current area                                | Proposed internal area                       |
| ------------------------------------------- | -------------------------------------------- |
| `reporting`, `nbversion`, `nbplatform`      | `numbat_core` top level (see note below)     |
| `nbtypes`                                   | `numbat_core.nbtypes`                        |
| `numbattools`                               | `numbat_core.utils.numeric`                  |
| `bulkprops`, `materialprops`, `crystalsyms` | `numbat_core.materials`                      |
| `voigt`                                     | `numbat_core.materials.elastic_tensors`      |

`numtools` does not exist in this repository and was not moved. `reporting`,
`nbversion`, and `nbplatform` sit at the top level of `numbat_core` rather than
under `runtime/`, because `numbat_core.runtime` is a module, not a package.
`nbtypes` was not in the original table. It moves so that `numbat_core` does
not import flat modules. `voigt_old.py` is left untouched.

For each move:

1. Move the implementation.
2. Update internal consumers to use the internal import path.
3. Retain a top-level facade, for example:

   ```python
   from numbat_core.materials.bulkprops import *
   ```

4. Run import and smoke tests.
5. Do not emit deprecation warnings yet.

### Phase 3: Consolidate mesh infrastructure

The following related components should eventually become `numbat_core.mesh`:

- `meshing/`;
- `msh/`;
- `nbgmsh.py`;
- `femmesh.py`;
- `usermesh.py`.

Suggested internal structure:

```text
numbat_core/mesh/
├── fem.py
├── mail.py
├── geometry.py
├── templates.py
├── loader.py
├── gmsh_io.py
└── builtin/
    ├── __init__.py                 # re-exports the builtin mesh classes
    ├── _common.py
    ├── basic.py
    ├── two_incl.py
    ├── onion.py
    ├── rib.py
    ├── slot.py
    ├── rounded.py
    ├── waveguides.json
    └── geo/
```

The built-in mesh classes have already been split from a single
`msh/builtin_meshes.py` into the `msh/builtin_meshes/` package, grouped by
waveguide family, with `__init__.py` re-exporting all 23 classes. Classes now
report their submodule as `__module__` (for example `msh.builtin_meshes.rib`).
The target layout above keeps the same family split.

Do not change built-in mesh metadata such as:

```json
"module": "msh.builtin_meshes"
```

in the first mesh migration. The `msh.builtin_meshes` package keeps this path
importable, so the JSON metadata stays valid. First make the loader accept both the current and
new module paths. This preserves template registration, user geometry pickle
identity, and compatibility with existing user code. Change the registry's
default module path only after the dual-path loader has shipped and is tested.

### Phase 4: Consolidate plotting

The existing `plotting/` package is already a useful subsystem. Move loose
top-level plotting modules into the internal plotting package:

| Current module         | Proposed internal module         |
| ---------------------- | -------------------------------- |
| `plotmodes.py`         | `numbat_core.plotting.modes`     |
| `plotmoderaw.py`       | `numbat_core.plotting.raw_modes` |
| `plotgain.py`          | `numbat_core.plotting.gain`      |
| `plotting/profiles.py` | `numbat_core.plotting.profiles`  |
| `plotting/gmsh.py`     | `numbat_core.plotting.mesh`      |

Keep mesh rendering separate from Gmsh subprocess and mesh-file I/O. The former
belongs in plotting; the latter belongs in `numbat_core.mesh`.

### Phase 5: Move simulation-domain modules

After mesh and plotting no longer import the application root, move the
simulation modules as a coherent group:

```text
numbat_core/
├── geometry/
│   ├── structure.py
│   └── usermesh.py
├── solvers/
│   ├── modes.py
│   ├── modecalcs.py
│   ├── integration.py
│   └── simresult.py
└── api/
    └── application.py
```

This phase includes the largest and most interconnected modules, so it should
follow rather than precede the low-level and mesh migrations.

### Phase 6: Introduce the modern public package

After the internal hierarchy is stable, introduce the preferred installed
package API, while retaining `import numbat` permanently if practical.

Potential modern imports:

```python
from numbat import Simulation
from numbat.materials import Material
from numbat.mesh import FemMesh
```

Deprecate specialist flat imports only over one or more releases, and update
tutorials incrementally rather than all at once.

## Do not do these things

- Do not require `import backend.numbat`; existing code intentionally treats
  `backend` as an import root.
- Do not move `msh` resources and alter their metadata in the same change.
- Do not move `femmesh.py` first; it currently has application and plotting
  coupling that should be removed beforehand.
- Do not use broad search-and-replace for imports.
- Do not mix this migration with behaviour changes, broad formatting, or
  unrelated typing work.

## Recommended first implementation milestone

The first restructuring branch should contain only:

1. import smoke tests;
2. `numbat_core.runtime`;
3. migration of reporting, platform/version, and numerical utility modules;
4. top-level compatibility facades;
5. conversion of a small number of internal consumers to the new paths.

This validates the migration mechanism before moving large modules such as
`femmesh.py`, `modecalcs.py`, and `materials.py`.
