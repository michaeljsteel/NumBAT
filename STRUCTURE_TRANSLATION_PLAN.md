# Structure translation implementation plan

Replace the separate EM and AC `set_xyshift_*` methods with one translation
specified when a structure is created. The translated geometry, mesh, solver
coordinates, plots, and exports should all use the same coordinate system.

## 1. Define the structure-level interface

Add an optional keyword to `nbapp.make_structure()` / `Structure`:

```python
wguide = nbapp.make_structure(
    "rectangular",
    domain_x=1500,
    domain_y=1500,
    inc_a_x=300,
    inc_a_y=280,
    translation=(200, -100),
    ...
)
```

- Default: `(0, 0)`.
- Meaning: **new coordinates = template coordinates + translation**.
- Units: the same as geometry dimensions—currently nm, changing to µm with
  the separate geometry-unit migration.
- Validate a pair of finite real numbers; reject malformed values explicitly.
- Store the translation on `Structure` as a common parameter, not a
  shape-specific one. Do not pass it through ordinary template parameter
  checks as an unrecognised geometry parameter.
- Apply the same translation to EM and AC; there will no longer be separate
  offsets.

## 2. Translate the generated Gmsh geometry

In [`Structure._build_mesh()`](backend/numbat/geometry/structure.py), after
`wg_geom.make_geometry()` returns:

1. Convert the physical translation to normalised Gmsh coordinates:
   ```python
   tx_normalised = tx / domain_x
   ty_normalised = ty / domain_x
   ```
2. Add a Gmsh translation to the generated `.geo` before running Gmsh.
3. Leave the zero-translation output unchanged.
4. Keep Gmsh's normalised coordinates and the existing metre scaling in Python
   and Fortran.

Putting this in the structure build rather than individual templates also
covers classes that override `make_geometry()`.

**First validate the Gmsh operation across supported templates.** A
point-translation probe worked for three examples, but it does not prove that
moving points after construction is safe for every Gmsh entity or external
script. In particular, check curved geometry, generated geometries, physical
groups, and any custom script that meshes or saves before the appended
operation. If point translation proves inadequate for a supported construction,
address it there rather than silently accepting a partially translated geometry.

The generated `.geo`, `.msh`, and `.mail` should all represent the translated
geometry. No Fortran interface change is planned.

## 3. Pass translation to wireframe hooks

Update the user hook to accept the structure translation:

```python
draw_mpl_frame(self, ax, styles, delta_xy)
```

The framework passes `delta_xy=(dx, dy)` in micrometres, matching the plotting
coordinates used by the hook. The hook applies the offset to its own geometry
before creating Matplotlib artists. Update built-ins, examples, and custom
waveguide documentation together. This keeps translation explicit and avoids
framework-side artist discovery and transform manipulation.

Test that the hook receives the correct offset and that representative built-in
and custom outlines align with the translated geometry.

## 4. Remove the old coordinate shifts

Once the mesh and wireframe share translated coordinates:

- Remove `set_xyshift_em()` and `set_xyshift_ac()` from
  [`structure.py`](backend/numbat/geometry/structure.py).
- Remove `get_xyshift()` from
  [`simresult.py`](backend/numbat/solvers/simresult.py).
- Remove the related arithmetic in:
  - [`mode_sampling.py`](backend/numbat/solvers/mode_sampling.py)
  - [`strain_sampling.py`](backend/numbat/solvers/strain_sampling.py)
  - [`modes.py`](backend/numbat/solvers/modes.py)
- Use translated FEM coordinates directly for 2D interpolation and 1D cuts.
  Convert units only where plots and cut inputs require it.
- Replace tests that encode the old shift behaviour with tests of the new
  shared-coordinate behaviour.

Leave `set_r0_offset()` and the width-reference API unchanged: they serve a
different purpose. Document that distinction.

## 5. Migration and compatibility

- Update the onion tutorial to provide one translation at construction.
- Document the sign, units, shared EM/AC coordinate system, and user-wireframe
  contract.
- Update the [micron migration plan](GEOMETRY_MICRONS_MIGRATION_PLAN.md) so it
  no longer proposes retaining the old setters.
- Check old saved simulations explicitly. A saved result with a nonzero old
  plotting shift has an **untranslated mesh**, so loading it as though it used
  the new model could produce misleading plots. Detect and report that case
  rather than silently reinterpret it; test zero-shift old results separately.
- Note in the migration guide that separate EM/AC plot origins no longer map
  directly to one structure translation.

## 6. Verification

### Unit and geometry tests

- Default and explicit zero translation; invalid inputs.
- Positive and negative translations on both axes, including a non-square
  domain to catch incorrect y scaling.
- Generated `.geo`, `.msh`, `.mail`, initial Python mesh, and solved mesh bounds.
- Built-in rectangular, circular, onion, rib/slot, and rounded geometries;
  external template and Python-generated geometry.
- Physical tags, material regions, and AC-to-EM mappings remain valid.

### Plot and sampling tests

- Wireframe hooks receive translation in micrometres and apply it to outlines.
- 2D mode and strain samples use the intended mesh points.
- x, y, and arbitrary-line cuts, including exported coordinate columns.
- Refractive-index/material profiles and mesh plots use the same coordinates.
- Test `set_r0_offset()` independently to ensure it is not applied twice.

### Real-solver checks

- Compare untranslated and translated versions of the same physical structure:
  EM effective indices, direct and EM-derived AC frequencies at nonzero `q_AC`,
  and a representative gain calculation.
- Compare within numerical tolerances rather than requiring identical node
  numbering: the Gmsh probe showed that translation can change remeshing details.
- Use matching Python/Fortran builds for the final integration run.

## Order of work

1. Add translation and Gmsh/mesh tests.
2. Implement the central wireframe transform and artist tests.
3. Remove downstream shift machinery and update sampling tests.
4. Update examples, docs, and saved-result handling.
5. Run representative real-solver and plotting checks.
