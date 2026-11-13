# Field Selection Migration Plan

## Goals

Replace the overlapping `FieldType`, `FieldCode`, and mutable `FieldTag`
responsibilities with:

- a typed field-family selector;
- a typed field component or derived-quantity selector;
- an explicit primary or secondary complex-part selector;
- one conversion layer between those types and internal `Fxr`/`Fyi`/`Fa`
  storage keys;
- plotting-only formatting and styling helpers.

The migration must preserve existing user code such as:

```python
result.plot_modes(field_type="EM_E", comps=("x", "y", "a"))
result.plot_modes(field_type=FieldType.EM_H, comps=("x", "z"))
```

Existing user-facing string forms remain supported throughout the compatibility
period.

## Target Interface

### Core types

Add these alongside the existing types in `backend/numbat/nbtypes.py`:

```python
class FieldComponent(Enum):
    X = "x"
    Y = "y"
    Z = "z"
    TRANSVERSE = "t"
    MAGNITUDE = "a"
    INTENSITY = "i"


class ComplexPart(Enum):
    PRIMARY = "primary"
    SECONDARY = "secondary"
```

Use `PRIMARY` and `SECONDARY` rather than `major` and `minor`. They represent
the current phase convention:

| Component | Primary storage key | Secondary storage key |
|---|---:|---:|
| `X` | `Fxr` | `Fxi` |
| `Y` | `Fyr` | `Fyi` |
| `Z` | `Fzi` | `Fzr` |

Add an immutable selection:

```python
@dataclass(frozen=True)
class FieldSelection:
    field_type: FieldType
    component: FieldComponent
    complex_part: ComplexPart | None = None
```

Validation invariants:

- `X`, `Y`, and `Z` require `complex_part`.
- `TRANSVERSE`, `MAGNITUDE`, and `INTENSITY` require `complex_part=None`.
- Invalid states raise a useful `ValueError` during construction.

Provide concise constructors:

```python
FieldSelection.primary(FieldType.EM_E, FieldComponent.X)
FieldSelection.secondary(FieldType.EM_H, FieldComponent.Z)
FieldSelection.magnitude(FieldType.AC)
FieldSelection.intensity(FieldType.EM_E)
FieldSelection.transverse(FieldType.EM_H)
```

`FieldSelection` must not own TeX formatting, Matplotlib styling, file paths,
NumPy data lookup, or simulation-domain coercion.

## Phase 0: Establish the Compatibility Contract

**Status: Complete**

### Objective

Document precisely what must remain supported before changing implementation.

### Work

1. Audit every public method accepting `field_type`, `comps`, direct
   `FieldTag` values, or raw `Fxr`/`Fyi`-style identifiers.
2. Record the compatibility matrix below as parameterized tests.
3. Keep all existing strings valid at public API boundaries.

| Existing input | Supported during migration | Normalized to |
|---|---:|---|
| `FieldType.EM_E`, `EM_H`, or `AC` | Yes | `FieldType` |
| `"EM_E"`, `"EM_H"`, or `"AC"` | Yes | `FieldType` |
| `"x"`, `"y"`, `"z"`, `"t"`, or `"a"` in `comps` | Yes | `FieldComponent` |
| `"i"` in APIs that support intensity | Yes | `FieldComponent.INTENSITY` |
| Existing `FieldTag` passed to internal plotting helpers | Yes, initially | Compatibility adapter |
| `Fxr`, `Fxi`, and similar data keys | Yes, internal only | Storage-key adapter |

Parsers must reject invalid values explicitly and list accepted values.

### Tests

Extend `tests/test_field_tags.py` with parameterized compatibility tests for
every currently accepted field type and component string.

### Exit criteria

- The compatibility matrix is represented by tests.
- No production behaviour has changed.

## Phase 1: Add Typed Value Objects and Parsers

**Status: Complete**

### Objective

Add the new model without changing public plotting behaviour.

### Work

1. Add `FieldComponent`, `ComplexPart`, and immutable `FieldSelection`.
2. Add:

   ```python
   FieldType.parse(value)
   FieldComponent.parse(value)
   ComplexPart.parse(value)
   ```

3. Preserve existing `FieldType` serialized values:

   ```python
   FieldType.EM_E.value == "EM_E"
   FieldType.EM_H.value == "EM_H"
   FieldType.AC.value == "AC"
   ```

4. Add `FieldSelection` properties:
   - `is_cartesian`;
   - `is_derived`;
   - `is_vector`.
5. Add factory and transformation methods:
   - `primary`;
   - `secondary`;
   - `magnitude`;
   - `intensity`;
   - `transverse`;
   - `with_part`;

### Tests

- Every parser accepts enums and supported strings.
- Invalid values have useful errors.
- Cartesian and derived selections validate correctly.
- Invalid component/complex-part combinations fail.
- Instances are immutable.
- Primary and secondary factory methods produce the expected values.

### Exit criteria

- New types are covered by unit tests.
- Public entry points still use existing types and strings.

## Phase 2: Centralize Representation and Storage Keys

### Objective

Eliminate scattered `Fxr`/`Fyi` mapping tables while preserving the existing
field-data dictionary format.

### Work

1. Add a focused module, preferably:

   ```text
   backend/numbat/field_selection.py
   ```

2. Place the only component/complex-part-to-storage-key mapping there:

   ```python
   storage_key(selection: FieldSelection) -> str
   selection_from_storage_key(
       field_type: FieldType,
       key: str,
   ) -> FieldSelection
   field_symbol(field_type: FieldType) -> str
   ```

3. Continue returning the current structure from
   `ModeInterpolator.interpolate_mode_i`:

   ```python
   {
       "Fxr": ...,
       "Fxi": ...,
       "Fyr": ...,
       "Fyi": ...,
       "Fzr": ...,
       "Fzi": ...,
       "Fa": ...,
   }
   ```

4. Convert callers to use:

   ```python
   fields[storage_key(selection)]
   ```

5. Replace internal use of `FieldTag.make_from_field_and_Fcode` with
   `selection_from_storage_key`.

### Tests

- Round-trip all valid selections through storage-key conversion.
- Assert exact existing key strings.
- Reject unsupported keys clearly.
- Test electric, magnetic, and acoustic display symbols separately.

### Exit criteria

- All storage-key mappings have one source of truth.
- Other mapping tables are removed or become thin compatibility wrappers.

## Phase 3: Introduce One Field-Resolution Policy

### Objective

Replace `FieldCode`'s implicit `force_AC` role with an explicit policy.

### Work

1. Add:

   ```python
   def resolve_field_type(
       requested: FieldType | str,
       *,
       simulation_is_acoustic: bool,
   ) -> FieldType:
       field_type = FieldType.parse(requested)
       return FieldType.AC if simulation_is_acoustic else field_type
   ```

2. Apply the resolver only at public boundaries:
   - `SimResult.plot_modes`;
   - `SimResult.plot_modes_1D`;
   - `SimResult.write_modes`;
   - `SimResult.write_modes_1D`;
   - direct `Mode` methods that can be invoked independently.
3. Keep lower-level helper functions typed to receive a resolved `FieldType`.
4. Retain `FieldCode` for one release as a deprecated wrapper that delegates
   to `resolve_field_type`.

### Tests

| Requested field | Simulation type | Expected resolved field |
|---|---|---|
| `EM_E` | EM | `EM_E` |
| `EM_H` | EM | `EM_H` |
| `AC` | AC | `AC` |
| `EM_E` | AC | `AC` |
| `EM_H` | AC | `AC` |

### Exit criteria

- Production code no longer needs `FieldCode`.
- Compatibility tests prove its legacy methods still work.

## Phase 4: Migrate Plotting to `FieldSelection`

### Objective

Remove mutable `FieldTag` from the plotting flow.

### Work

1. Replace `_field_tags_for_1D_plot` with:

   ```python
   def selections_for_1d_plot(
       field_type: FieldType,
       components: Sequence[FieldComponent],
       *,
       include_secondary: bool,
   ) -> list[FieldSelection]:
       ...
   ```

2. Preserve 1D selection order:

   ```text
   components=(X, Y, MAGNITUDE), include_secondary=True
   -> [X primary, Y primary, magnitude, X secondary, Y secondary]
   ```

3. Migrate scalar contours, vector plots, all-components plots, and
   single-component plots to accept `FieldSelection`.
4. Move presentation helpers out of `FieldTag` into plotting code:

   ```python
   tex_label(selection)
   line_color(selection)
   line_style(selection, *, show_magnitude)
   field_display_name(field_type)
   domain_display_name(field_type)
   ```

5. Add `FieldTag.as_selection()` as a compatibility adapter.

### Tests

- 1D selection order.
- Primary and secondary labels.
- Labels for `E`, `H`, and `u`.
- Styles for primary, secondary, magnitude, intensity, and vector quantities.
- Storage-key lookup for rendered scalar quantities.
- No secondary representation for `t`, `a`, or `i`.

### Exit criteria

- Plotting code uses `FieldSelection` directly.
- Plotting no longer accesses `FieldTag` private attributes.

## Phase 5: Upgrade Public Component Inputs

### Objective

Make enums the documented API while retaining existing strings.

### Work

1. Update public documentation and annotations to conceptually accept:

   ```python
   comps: Sequence[FieldComponent | str] = ()
   ```

2. Parse public component input once at the boundary:

   ```python
   components = tuple(FieldComponent.parse(component) for component in comps)
   ```

3. Preserve default components exactly:

   ```python
   (
       FieldComponent.X,
       FieldComponent.Y,
       FieldComponent.Z,
       FieldComponent.MAGNITUDE,
   )
   ```

4. Decide whether public intensity support is ready:
   - If yes, support `comps=("i",)`, calculate `Fi`, and document it.
   - If no, reject `i` at public plotting boundaries with a specific message
     while retaining the type for internal future use.
5. Update relevant docstrings and plotting tutorials to show enum-first usage,
   while documenting string compatibility.

### Tests

For each supported public plot API, test an old string call and an enum call
that resolve to identical selections, storage-key access, filenames, labels,
and plot metadata.

### Exit criteria

- New documentation uses enums first.
- Existing string inputs remain covered by regression tests.

## Phase 6: Deprecate Legacy Types

### Objective

Give downstream users a stable transition path before removal.

### Work

1. For one minor release:
   - keep `FieldCode` importable and functional;
   - keep `FieldTag` importable and functional;
   - prevent new production code from constructing either type;
   - emit deprecation warnings only for direct external use.
2. Use `DeprecationWarning` with an appropriate `stacklevel`, or a
   project-specific warning class if one exists.
3. Publish this migration table:

| Legacy API | Replacement |
|---|---|
| `FieldCode(field_type, force_AC)` | `resolve_field_type(field_type, simulation_is_acoustic=...)` |
| `FieldTag.make_from_field_and_component(ft, "x")` | `FieldSelection.primary(ft, FieldComponent.X)` |
| `tag.set_to_minor()` | `selection.with_part(ComplexPart.SECONDARY)` |
| `tag.component_as_F()` | `storage_key(selection)` |
| `tag.get_tex_plot_label()` | `tex_label(selection)` |
| `tag.linecolor()` | `line_color(selection)` |
| `tag.linestyle(comps)` | `line_style(selection, show_magnitude=...)` |

4. Remove legacy classes only after:
   - the stated deprecation window;
   - all production code uses the new API;
   - examples and documentation use the new API;
   - a released version has exposed warnings;
   - removal is covered by regression tests.

### Exit criteria

- The public migration window has elapsed.
- `FieldCode` and `FieldTag` can be removed without repository callers.

## Test Strategy

### Unit tests

Expand `tests/test_field_tags.py` during the transition. Rename it to
`tests/test_field_selection.py` when it primarily covers the new interface.

Cover:

- parser normalization;
- invalid values;
- selection invariants;
- primary/secondary mapping;
- all storage keys;
- TeX and plain-text labels;
- line styles and colours;
- legacy adapters;
- deprecation warnings.

### Integration tests

Add targeted tests for:

- electric EM plots;
- magnetic EM plots;
- acoustic displacement plots;
- default 1D component plots;
- 1D plots with secondary components;
- magnitude and vector-only requests;
- preserved plotting and writing filenames.

Prefer mocked interpolated data for presentation tests. Keep one end-to-end
solver-backed test per simulation family.

### Compatibility tests

Maintain parameterized old/new equivalence tests:

```python
old = result.plot_modes(field_type="EM_E", comps=("x", "a"))
new = result.plot_modes(
    field_type=FieldType.EM_E,
    comps=(FieldComponent.X, FieldComponent.MAGNITUDE),
)
```

Assert equivalent resolved selections, storage accesses, filenames, labels,
and plot metadata.

## Suggested Release Sequencing

| Release | Deliverables |
|---|---|
| Release A | Phases 0-2: new immutable types, parsers, central mapping; no public behaviour change |
| Release B | Phases 3-5: migrate production plotting and public parsing; document enum-first API |
| Release C | Phase 6 deprecation release: warnings for direct `FieldCode` and `FieldTag` use |
| Release D | Remove legacy classes after the documented deprecation window |
