# Mode Plotting Separation Plan

## Purpose

Separate mode-domain, field-sampling, rendering, and export responsibilities
currently combined in `backend/numbat/solvers/modes.py`.

This work follows completion of
[FIELD_SELECTION_MIGRATION_PLAN.md](./FIELD_SELECTION_MIGRATION_PLAN.md). It
must preserve the public plotting APIs and the field-selection compatibility
guarantees established by that migration.

## Current Responsibilities

`backend/numbat/solvers/modes.py` currently combines:

1. **Field sampling and interpolation**
   - Builds 1D and 2D sampling grids.
   - Interpolates FEM fields.
   - Applies magnetic `Z0` scaling.
   - Produces component arrays including `Fa` and `Fi`.

2. **Mode-domain analysis**
   - Computes polarisation fractions.
   - Computes centre of mass and second-moment widths.
   - Stores mode metadata and analysis results.

3. **Plot-output policy and rendering**
   - Interprets `comps`, including the `e` overview directive.
   - Creates 1D Matplotlib figures directly.
   - Selects overview and individual 2D output.
   - Saves figures.

4. **Field export**
   - Writes 1D and 2D field arrays to text files.

The plotting package already owns lower-level 2D presentation:

- contour and quiver rendering;
- plot labels, colour choices, and line styles;
- multi-panel and individual 2D figure layout;
- plot preferences, axes styling, and figure-save helpers.

## Target Module Boundaries

```text
solvers/
    modes.py                Public Mode, ModeEM, and ModeAC domain façade
    mode_sampling.py        Field interpolation and sampling-grid management

plotting/
    mode_presentation.py    Future shared labels, styles, and display helpers
    modes_2d.py             Future low-level 2D mode drawing primitives
    mode_rendering.py       New 1D/2D mode render orchestration
    modes.py                Compatibility re-export module
    raw_modes.py            Existing raw FEM plotting

io/
    mode_fields.py          New field-data export helpers
```

### `solvers/modes.py`

Retain:

- `Mode`, `ModeEM`, and `ModeAC`;
- mode-type checks and mode-domain metadata;
- analysis methods and analysis results;
- public API compatibility wrappers:
  - `plot_mode`;
  - `plot_mode_1D`;
  - `plot_mode_H`;
  - `write_mode`;
  - `write_mode_1D`;
  - `get_mode_1D`.

It must not import `matplotlib.pyplot` or create figures directly after the
rendering extraction is complete.

### `solvers/mode_sampling.py`

Move `ModeInterpolator` here, retaining its current sampling semantics.
Consider renaming it to `ModeFieldSampler` only with a documented
compatibility alias.

This module is solver/data infrastructure, not plotting code, because its
interpolated arrays are consumed by:

- mode analysis;
- 1D and 2D plotting;
- field export;
- programmatic `get_mode_1D` access.

It owns:

- 1D and 2D grid creation;
- FEM interpolation setup;
- field interpolation;
- field-unit scaling;
- canonical generic field arrays (`Fxr`, `Fxi`, ..., `Fa`, `Fi`).

It must not import Matplotlib or plotting helpers.

### `plotting/mode_rendering.py`

Add a new mode-rendering orchestration module rather than expanding
`plotting/modes.py`.

It owns:

- 1D figure construction, curves, axes, legends, and figure saving;
- 2D overview-versus-individual-output routing;
- conversion from requested plot components into `FieldSelection` values;
- 2D-only `EVERYTHING`/`"e"` output semantics;
- supplied-axis restrictions.

It should depend on field arrays and coordinate data supplied by the solver
layer rather than on interpolation internals.

Suggested rendering boundaries:

```python
render_mode_1d(...)
render_mode_2d(...)
resolve_mode_plot_request(...)
```

The render functions should receive typed field selections, field arrays,
coordinates, `PlotParams`, and only the mode/simulation metadata needed for
titles and overview annotations.

### `plotting/mode_presentation.py`

Eventually extract shared, dimension-independent presentation helpers:

- `FieldSelection` labels;
- line colours and styles;
- field and domain display names;
- field-component-to-selection conversion;
- shared mode-plot filename construction.

This keeps `modes_2d.py` accurately named while supporting 1D rendering
without a dependency on 2D drawing code.

### `plotting/modes_2d.py`

Rename the current low-level 2D drawing implementation to `modes_2d.py` only
after the shared presentation helpers have been extracted. It will own:

- contour and quiver drawing;
- overview layout;
- individual-component drawing.

### `plotting/modes.py`

Retain `modes.py` as a compatibility re-export module after the rename:

- update NumBAT's internal imports to use the explicit new modules;
- preserve existing external `numbat.plotting.modes` imports;
- defer any import-path deprecation decision to a later release.

Do not move 1D orchestration into `modes_2d.py` merely because it is
mode-related; that would make the renamed module another catch-all.

### `io/mode_fields.py`

Move text serialization out of the solver module:

- 2D component writer;
- 1D multi-component writer;
- filename/header construction tied to field data exports.

The export schema must explicitly decide when to include derived fields such
as `Fa`, `Fi`, and `Ft`; it must not change accidentally as a side effect of
plotting refactors.

## Component Specification Policy

The generic component parser remains in `numbat.nbtypes`:

```python
parse_field_components("xyat")
parse_field_components((FieldComponent.X, FieldComponent.Y))
```

Plot-output policy is presentation-specific and moves to
`plotting/mode_rendering.py`:

- no `comps`: combined 2D overview only;
- `comps="xyt"`: individual x, y, and transverse plots only;
- `comps="exyt"`: combined overview plus those individual plots;
- `EVERYTHING`/`"e"`: valid for 2D rendering only;
- `EVERYTHING`/`"e"`: rejected for 1D rendering;
- a supplied Matplotlib axis: exactly one individual component and no
  overview.

`FieldComponent.EVERYTHING` is retained as the user-facing compact grammar.
It is a rendering directive, not a physical field quantity, so it must be
removed before `FieldSelection` construction or generic storage-key lookup.

## Migration Phases

### Phase 0: Characterise Current Behaviour

Before moving code:

- add or retain focused tests for 1D and 2D field selection;
- cover default overview-only output;
- cover individual-only output;
- cover overview plus individual output through `e`;
- cover supplied-axis restrictions;
- record representative filenames and field labels.

### Phase 1: Extract 1D Rendering

**Status: Complete**

Move direct Matplotlib construction from `Mode._plot_mode_1D_gencut` into
`plotting/mode_rendering.py`.

The solver-facing method remains as a compatibility wrapper that:

1. resolves the field type;
2. obtains interpolated field data;
3. delegates to `render_mode_1d`.

Success criterion: `solvers/modes.py` no longer imports `matplotlib.pyplot`.
Keep shared presentation helpers in `plotting/modes.py` during this phase;
their extraction precedes the later `modes_2d.py` rename.

### Phase 2: Extract 2D Render Orchestration

**Status: Complete**

Move `Mode._plot_me_2D` output routing to `render_mode_2d`.

Keep current contour/quiver and overview primitives in `plotting/modes.py`.
The new renderer orchestrates those primitives using the requested overview
flag and individual component selections.

Success criterion: solver code supplies data and requests rendering, but does
not select decorators or call plotting primitives directly.

### Phase 3: Move Plot-Request Policy

**Status: Complete**

Move plot-specific `comps` interpretation from `solvers/modes.py` into
`plotting/mode_rendering.py`.

Retain `parse_field_components` in `nbtypes.py`; do not duplicate parsing or
component validation.

Success criterion: field-component normalization is shared, while overview
semantics are owned by the renderer.

### Phase 4: Extract Field Sampling

**Status: Complete**

Move `ModeInterpolator` to `solvers/mode_sampling.py`.

Preserve:

- interpolation cache behavior;
- 1D and 2D coordinate conventions;
- primary/secondary z convention;
- H-field `Z0` scaling;
- `Fi = Fa**2` after scaling;
- consumers in plotting, analysis, export, and `get_mode_1D`.

Provide a compatibility import or alias if external users can import
`ModeInterpolator` from its current location.

### Phase 5: Extract Field Export

Move `_write_one_component_to_file`, `write_mode`, and `write_mode_1D`
serialization work into `io/mode_fields.py`.

Decide and document the exported derived-field schema, including whether
`Fi` is emitted for 1D and/or 2D exports.

### Phase 6: Clean Up and Document

- remove now-unused solver-to-plotting imports;
- update module and class docstrings;
- document architectural boundaries for contributors;
- deprecate compatibility aliases only in a future release, if warranted.

## Compatibility Requirements

The following must remain supported throughout the migration:

- existing `Mode` and `SimResult` plotting method names and arguments;
- legacy component strings;
- enum component values;
- compact component strings;
- default 2D overview-only behavior;
- explicit `e` overview behaviour;
- existing plot output filenames and labels unless a deliberate release-note
  change is approved;
- existing field storage keys and `FieldSelection` mappings.

Avoid moving the public API at the same time as changing its semantics.
Extract implementation behind compatibility wrappers first.

## Validation Strategy

For each phase:

1. Run field-selection regression tests.
2. Run focused tests for the extracted renderer, sampler, or exporter.
3. Verify output routing with mocked drawing functions:
   - overview only;
   - individual only;
   - overview plus individual;
   - invalid supplied-axis and 1D overview requests.
4. Verify field-data values:
   - primary and secondary component keys;
   - z primary/secondary convention;
   - `Fa`;
   - `Fi = Fa**2`;
   - H-field scaling before intensity calculation.
5. Run syntax/lint checks for touched modules.
6. Run representative example scripts where practical.

## Non-Goals

- Replacing Matplotlib.
- Redesigning field-selection types during this migration.
- Changing numerical interpolation algorithms.
- Changing plot appearance, default filenames, or user preferences without an
  explicit separate decision.
- Combining the rendering and export refactors into one large change.
