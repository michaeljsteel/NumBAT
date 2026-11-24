# Geometry dimensions in micrometres: migration plan

## Goal

Make geometry lengths supplied to `make_structure()` consistently
micrometres:

```python
# Current: nm
wguide = nbapp.make_structure("rectangular", 1500, 1500, 300, 280, ...)

# After migration: µm
wguide = nbapp.make_structure("rectangular", 1.5, 1.5, 0.3, 0.28, ...)
```

Keep the numerical solver and solved FEM meshes in metres, and keep Gmsh mesh
coordinates normalised by `domain_x`, as they are now.

This is an API-breaking change: the existing parameter names do not encode
their units, so old calls would otherwise run with a geometry 1000 times larger.

## Recommended approach

Make the switch explicitly rather than guessing units from the size of
`domain_x`. A large value could be either an old nanometre input or a
legitimately large micron-scale structure.

Make geometry classes themselves work in µm rather than converting inputs back
to nm internally. This gives user-defined geometry code, built-in geometry code,
and plot outlines the same convention.

### Compatibility decision

Before implementation, decide whether to provide an explicit temporary legacy
mode for old scripts or make this a clean breaking change with migration
instructions. A legacy mode is possible, but arbitrary user-defined geometry
parameters make automatic conversion tricky: NumBAT cannot know which custom
parameters are lengths and which are angles, counts, or ratios without
additional metadata. Prefer a clean break unless preserving old scripts is
important enough to justify defining that metadata.

## 1. Establish and document unit boundaries

| Quantity | Unit after migration |
|---|---|
| `domain_x`, `domain_y`, inclusion sizes, radii, gaps, offsets and other geometry lengths | µm |
| Geometry-class parameters and `draw_mpl_frame` coordinates | µm |
| `lc_bkg`, refinement factors and geometric ratios | Dimensionless; unchanged |
| Gmsh-generated coordinates and `.mail` coordinates | Normalised by `domain_x`; unchanged |
| `Simulation.d_in_m`, solved EM/AC FEM coordinates and solver inputs | m |
| Mode-plot coordinates and cut positions | µm, as today |
| Wavenumbers, frequencies and material constants | Existing units; unchanged |

This distinction matters because
[`FemMesh.build_from_gmsh_mail`](backend/numbat/mesh/fem.py) currently creates
an initial Python mesh in µm for plotting, whereas solved meshes are in metres.
The direct AC path converts that initial mesh to metres.

Do not automatically simplify every instance of `0.001`, `1000`, or `SI_um`:
many belong to unrelated physics or output conversions.

## 2. Change Python-to-solver scaling

In [`modecalcs.py`](backend/numbat/solvers/modecalcs.py):

- Change `d_in_m` from `domain_x * SI_nm` to `domain_x * SI_um`.
- Continue passing `d_in_m` to Fortran. The Fortran mesh loader already
  expects a scale in metres, so no Fortran unit change should be needed.
- Verify that every consumer of `d_in_m` continues to receive metres.

In [`fem.py`](backend/numbat/mesh/fem.py):

- Initial `.mail` mesh coordinates become `normalised_coordinate * domain_x`,
  producing µm directly. Remove the current `0.001` conversion.
- Retain the direct AC path's `* SI_um` conversion to metres.
- Update comments and docstrings to state units explicitly. Preserve the fact
  that both x and y normalised coordinates use `domain_x` as the scale.

This keeps the direct AC fix intact and avoids altering the plotting mesh's
current units.

## 3. Migrate geometry parameters and offsets

In [`structure.py`](backend/numbat/geometry/structure.py):

- Update the `make_structure`/`Structure` contract, parameter help, and examples
  to state µm.
- Ensure positional and keyword arguments use the same new convention,
  including the old `unitcell_x` and `unitcell_y` aliases if retained.
- Update the creation-time `translation=(x, y)` argument to use µm, and
  adjust its conversion to normalised Gmsh coordinates and wireframe plotting.
  The former separate EM/AC plotting-shift setters have been removed.
- Update examples that provide translation or derive it from `domain_x`.

The migration must cover more than standard inclusion dimensions: slab widths
and heights, coating thicknesses, gaps, separations, radii, offsets, and
shape-specific lengths all need the same treatment. Angles such as `alpha_tip`
and `rot_deg`, and mesh refinement parameters, must not be scaled.

## 4. Built-in geometry classes and outlines

The built-in classes pass geometry parameters through to Gmsh and use them to
draw outlines. In particular,
[`basic.py`](backend/numbat/mesh/builtin/basic.py) and other built-in modules
multiply lengths by `nmtoum` for Matplotlib.

- Remove the now-unnecessary nm-to-µm conversion from built-in
  `draw_mpl_frame()` implementations.
- Update [the shared helper](backend/numbat/mesh/builtin/_common.py), its
  imports, and comments.
- Check every built-in class: rectangular/circular, ribs, slots, onions,
  two-inclusion, and rounded shapes.
- Preserve dimensionless geometry calculations and angle units.
- Check generated geometries as well as `.geo`-based ones: rounded wedge
  variants construct normalised Gmsh coordinates directly in
  [`rounded.py`](backend/numbat/mesh/builtin/rounded.py).
- Check tolerances or comparisons involving dimensional values; hard-coded
  length tolerances may also need conversion.

An outline drawn for a 0.3 µm core should still line up with the field or
refractive-index plot of that core.

## 5. Gmsh templates

Most templates normalise dimensional parameters by `domain_x`, so geometry
formulas should generally remain dimensionless and unchanged. They are not
automatically safe just because ratios are used.

For each active template and its corresponding Python substitutions:

1. Inventory all dimensional defaults, including optional parameters that may
   not be supplied by a user.
2. Convert those defaults from nm to µm together. Otherwise, a substituted
   `domain_x=1.5` combined with an old default such as `a1=20` would describe a
   different geometry.
3. Rename misleading variables such as `dx_in_nm` and `dy_in_nm` to
   `dx_in_um` / `dy_in_um`, or unit-neutral names with clear comments.
4. Update the exact replacement strings in `apply_parameters()` at the same
   time. [`UserGeometryBase.make_geometry`](backend/numbat/geometry/usermesh.py)
   uses literal string replacement, so changing only the `.geo` file would
   silently leave defaults in place.
5. Check positional defaults, derived values, and any hard-coded dimensional
   literals separately from dimensionless constants, physical IDs, and
   refinement values.
6. Review `%f` substitution precision. It currently emits six decimal places;
   in µm this has a coarser physical resolution than it did in nm. Use
   sufficient significant digits for small dimensions and add a test that
   catches rounding to zero or losing meaningful precision.

There are 29 built-in `.geo` files in the tree; not all are necessarily active.
The inventory should distinguish active, shared, legacy, and example templates
rather than blindly rewriting them all.

## 6. User-defined waveguides and `.geo` files

The external geometry mechanism passes the parameter dictionary directly to
the user's class; see [`structure.py`](backend/numbat/geometry/structure.py).
NumBAT cannot automatically fix custom drawing or hard-coded dimensions in
external code.

Update the supplied example as the model for migration:

- [`elliptical_rod.py`](examples/user_waveguides/elliptical_rod.py): remove
  `/ 1000` from the outline and update substitution strings.
- [`elliptical_rod_msh_template.geo`](examples/user_waveguides/elliptical_rod_msh_template.geo):
  convert dimensional defaults and rename nm-labelled variables.
- [`waveguides.rst`](docs/source/waveguides.rst): explain the new contract,
  including that custom length parameters and plot outlines use µm, while
  generated mesh coordinates should retain the expected normalisation.

The migration guide should show both Python class and `.geo` changes. Tell users
to convert defaults as well as values provided by scripts, and to preserve
dimensionless parameters.

## 7. Examples, tests, and saved data

### Examples and documentation

- Convert geometry lengths in active tutorials, literature examples, and tests.
  At least 44 example/test files contain `make_structure()` calls.
- Review scans and derived quantities manually: changing an input array from nm
  to µm also affects plot labels, table columns, filenames, analytic formulas,
  and conversions used later in the script.
- Do not change wavelength or wavenumber calculations merely because they
  contain `1000` or `1e-9`.
- Update [`tutorial.rst`](docs/source/tutorial.rst),
  [`usage.rst`](docs/source/usage.rst), waveguide documentation, parameter help,
  and migration notes. Some wavelength text already appears out of date and
  should be made consistent while touching the units documentation.

### Saved simulations

Saved simulations use NumPy/pickle-based persistence in
[`modecalcs.py`](backend/numbat/solvers/modecalcs.py). Old saved objects may
contain a structure and user geometry parameters in nm, even though their solved
mesh coordinates are already in metres.

Before promising backward compatibility:

- Test loading an old saved EM and AC result, plotting it, and using its
  structure for a new calculation.
- Decide whether to support migration of old saved structures. Do not blindly
  multiply every numeric field: arrays, shifts, mesh coordinates, material data,
  and arbitrary custom parameters have different units.
- If reliable conversion cannot be established, document the limitation and
  emit a clear version/unit error where practical rather than silently
  producing a mis-scaled new calculation.

## 8. Verification

Build tests around physical equivalence, not just changed strings:

1. **Baseline before the switch:** record representative old-nm geometry and
   solver results, using the same mesh settings.
2. **Template tests:** compare normalised generated geometry for equivalent
   old/new dimensions, covering optional defaults, very small lengths, and
   user-defined geometry.
3. **Geometry coverage:** generate meshes for representative rectangular,
   circular, rib/slot, onion, rounded, and external templates; check bounds,
   connectivity/material regions, and ratios.
4. **Plotting:** check outline dimensions, material plot axes, mode-plot limits,
   and shared structure translation in µm.
5. **Mesh paths:** check the initial plotting mesh is in µm; solved EM mesh and
   both direct and EM-derived AC meshes are in metres. Extend
   [`test_ac_mesh_units.py`](tests/test_ac_mesh_units.py) with new-unit inputs.
6. **Solver integration:** compare EM effective indices and AC frequencies for
   equivalent physical structures, including AC-only and EM-derived AC at
   nonzero `q_AC`. Allow sensible numerical tolerances for eigensolver ordering
   and degeneracies.
7. **Migration checks:** run representative updated tutorials and test whichever
   legacy/saved-result policy is chosen.

## Suggested order

1. Decide the breaking-change/legacy policy and saved-result expectations.
2. Capture baseline cases and add equivalence tests.
3. Change the core unit boundaries and structure translation.
4. Migrate built-in geometry classes and templates together.
5. Migrate the user-waveguide example, documentation, and active examples.
6. Run mesh, plotting, and real-solver comparisons; address remaining cases
   exposed by the wider suite.
