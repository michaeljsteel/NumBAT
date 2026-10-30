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
    ├── builtin_waveguides.json
    └── geo/
```

The built-in mesh classes have already been split from a single
`msh/builtin_meshes.py` into the `msh/builtin_meshes/` package, grouped by
waveguide family, with `__init__.py` re-exporting all 23 classes. Classes now
report their submodule as `__module__` (for example `msh.builtin_meshes.rib`).
The target layout above keeps the same family split.

The legacy file-based `wg_impl` mechanism, its dynamic loader, and
`backend/msh/user_waveguides.json` have been removed. User waveguides are
normal Python modules in registered out-of-tree directories. Consequently,
this phase does not need a file loader, file-based pickle reconstruction, or
compatibility resources for user models.

Migrate this phase in the following order:

1. Move the implementations of `usermesh.py`, `nbgmsh.py`, `femmesh.py`, and
   `meshing/templates.py` to `numbat_core.mesh`, then replace their top-level
   modules and the `meshing` package with re-exporting facades.
2. Move the built-in family modules to `numbat_core.mesh.builtin`, retaining
   `msh.builtin_meshes` as a facade that re-exports every built-in class.
3. Package the built-in manifest and `.geo` files with the internal mesh
   package. Provide a resource-path helper so built-in template discovery does
   not rely on `backend/msh` being the implementation location.
4. Update internal consumers to import `numbat_core.mesh.*`; preserve the
   existing flat imports and `msh.builtin_meshes.*` imports for tutorials,
   existing scripts, and serialized built-in geometries.
5. Test the built-in manifest, a registered external user-waveguide directory,
   Gmsh template resolution, mesh generation, and geometry pickle round trips.

Do not change the built-in manifest's initial module path:

```json
"module": "msh.builtin_meshes"
```

The facade keeps this metadata valid while implementations move internally.
The template loader now uses ordinary `import_module()` calls, so it already
accepts either the legacy facade path or `numbat_core.mesh.builtin`; test both
paths with equivalent manifests. Change the shipped manifest to the internal
path only in a later, explicit compatibility release.

### Phase 4: Consolidate plotting

Move all plotting implementations to `numbat_core.plotting`:

```text
numbat_core/plotting/
├── axes.py
├── fields.py
├── gain.py
├── ipython.py
├── materials.py
├── mesh.py
├── modes.py
├── preferences.py
├── profiles.py
├── raw_modes.py
├── tools.py
└── resources/
    └── gmsh/
```

The legacy flat modules and `plotting/` package remain compatibility aliases:

| Legacy path              | Internal implementation                |
| ------------------------ | -------------------------------------- |
| `plotmodes.py`           | `numbat_core.plotting.modes`           |
| `plotmoderaw.py`         | `numbat_core.plotting.raw_modes`       |
| `plotgain.py`            | `numbat_core.plotting.gain`            |
| `plotting/plotprefs.py`  | `numbat_core.plotting.preferences`     |
| `plotting/plottools.py`  | `numbat_core.plotting.tools`           |
| `plotting/tidyaxes.py`   | `numbat_core.plotting.axes`            |
| `plotting/profiles.py`   | `numbat_core.plotting.profiles`        |
| `plotting/gmsh.py`       | `numbat_core.plotting.mesh`            |
| `msh/gmsh.py`            | `numbat_core.plotting.mesh`            |

Before moving the modules, remove plotting behavior from `numbat_core.mesh`.
`FEMScalarFieldPlotter` belongs in `numbat_core.plotting.fields`, and MAIL
mesh rendering belongs in `numbat_core.plotting.mesh`; `FemMesh` and
`MailData` remain pure mesh representations. Move MSH file parsing to
`numbat_core.mesh.gmsh_io`, then have plotting consume that API.

Update all implementation imports to use `numbat_core.plotting.*` and retain
no `numbat_core.plotting` imports under `numbat_core.mesh`. Add headless
rendering tests and verify both legacy and internal module paths resolve to
the same module objects.

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
