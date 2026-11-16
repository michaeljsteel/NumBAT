# Acoustic Strain Post-Processing Plan

## Purpose

Add accurate acoustic strain post-processing, plotting, and export for
NumBAT acoustic modes.

The feature must calculate strain from the solved displacement on the native
six-node quadratic FEM elements before any interpolation to rectangular plot
grids or one-dimensional cuts. It must support:

- all six independent symmetric strain components;
- combined and individual two-dimensional plots;
- one-dimensional line-cut plots;
- two-dimensional and one-dimensional text exports;
- the existing plot trimming, colours, frame, title, and file-prefix options
  where applicable.

The intended physical definition is:

\[
S_{ij} = \frac{1}{2}
\left(
\frac{\partial u_i}{\partial x_j} +
\frac{\partial u_j}{\partial x_i}
\right).
\]

## Scope and Constraints

- Strain is an acoustic-only quantity. EM simulations must reject strain
  requests with a clear error.
- The source displacement is the native acoustic `fem_evecs` array, with
  shape `[component, local_node, mode, element]`.
- FEM elements have six local P2 nodes. Derivatives must use the P2 basis and
  the physical isoparametric Jacobian, not finite differences on an
  interpolated rectangle.
- Native strain can be discontinuous across material interfaces. Do not
  average values shared by neighbouring elements.
- Rectangular-grid plots and files, and line cuts, are presentation/output
  samples of the native strain field. They are not the source of the strain
  calculation.
- The obsolete rectangular-grid strain code in `backend/zz_spare.py` must
  not be restored. The broken current `Mode.plot_strain` hook must be
  replaced only after the new public API exists.

## Established Project Facts

### FEM Geometry and Sampling

- `FemMesh.m_elnd_to_mshpt` maps six local nodes to physical mesh nodes.
- `FemMesh.v_mshpt_xy` provides native physical x/y coordinates.
- `FemMesh.make_sub_triangulation()` splits each six-node P2 element into
  four display/interpolation triangles.
- `FemMesh.make_interpolator_for_grid()` already builds Matplotlib
  triangulation interpolation for rectangular grids and line cuts.
- `SimResult.get_mode_fields_on_fem_mesh()` demonstrates the existing
  component/node/mode/element ordering needed to extract element-local
  values.

### Longitudinal Convention and Voigt Order

The guided acoustic displacement convention is:

\[
u(x,y,z) = \tilde u(x,y)e^{iq_{\rm AC}z},
\qquad
\partial_z = i q_{\rm AC}.
\]

Use the established Voigt order:

\[
[S_{xx}, S_{yy}, S_{zz}, S_{yz}, S_{xz}, S_{xy}].
\]

The order and Fourier convention are consistent with
`backend/nbanalytic/elasticfieldconversions.py`.

For each element-local displacement field:

\[
\begin{aligned}
S_{xx} &= \partial_x u_x, \\
S_{yy} &= \partial_y u_y, \\
S_{zz} &= i q_{\rm AC} u_z, \\
S_{yz} &= \frac{1}{2}(\partial_y u_z + i q_{\rm AC} u_y), \\
S_{xz} &= \frac{1}{2}(\partial_x u_z + i q_{\rm AC} u_x), \\
S_{xy} &= \frac{1}{2}(\partial_y u_x + \partial_x u_y).
\end{aligned}
\]

Strain is dimensionless when displacement and native coordinates are both
expressed in SI units.

### Default Complex-Part Presentation

The acoustic displacement display-primary convention is:

```text
ux: real
uy: real
uz: imaginary
```

The corresponding default strain presentation components are:

```text
Re[Sxx]  Re[Syy]  Re[Szz]
Im[Syz]  Im[Sxz]  Re[Sxy]
```

This mapping must be explicit, not inferred from the existing Cartesian
displacement component rules:

| Component | Display-primary part |
|---|---|
| `Sxx` | real |
| `Syy` | real |
| `Szz` | real |
| `Syz` | imaginary |
| `Sxz` | imaginary |
| `Sxy` | real |

## Target Architecture

### New Typed Selection Model

Do not overload `FieldType`, `FieldComponent`, or `FieldSelection` for a
rank-two strain tensor. Introduce a dedicated model, initially in a new
module such as `numbat/strain_selection.py`:

```python
class StrainComponent(Enum):
    XX = "xx"
    YY = "yy"
    ZZ = "zz"
    YZ = "yz"
    XZ = "xz"
    XY = "xy"

@dataclass(frozen=True)
class StrainSelection:
    component: StrainComponent
    complex_part: ComplexPart
```

The module owns:

- component parsing;
- display-primary and secondary mappings;
- storage-key conversion;
- display label generation;
- deterministic component ordering.

Use unambiguous storage keys:

```text
Sxxr, Sxxi
Syyr, Syyi
Szzr, Szzi
Syzr, Syzi
Sxzr, Sxzi
Sxyr, Sxyi
```

### Native Strain Evaluation

Create `numbat/solvers/strain_sampling.py`, separate from
`mode_sampling.py`.

Its responsibilities:

1. Validate that the source simulation is acoustic.
2. Extract `ux`, `uy`, and `uz` element-local P2 nodal values for one mode.
3. Extract the matching six local x/y coordinates for every element.
4. Evaluate P2 basis-function derivatives in reference coordinates.
5. Construct the P2 isoparametric physical Jacobian at each requested native
   evaluation point.
6. Transform reference derivatives to physical x/y derivatives.
7. Form all six complex strain components, including the `i*q_AC`
   longitudinal terms.
8. Return an element-local native strain representation without averaging
   across elements.

The implementation must handle curvilinear P2 elements. It must not assume
that every six-node element is an affine triangle.

### Native Values and Sampling

The first native representation should retain values per element/local
evaluation point. It must retain duplicate values at shared nodes so a
material-interface jump is preserved.

For rectangular-grid and line-cut output:

1. Locate the containing native triangle/element.
2. Evaluate or interpolate the already-calculated native strain within that
   element.
3. Return a deterministic one-sided value if the sample lies exactly on a
   material interface.

Document this one-sided interface behavior in the public API and export
format documentation.

The existing Matplotlib triangulation support may be used for the sampling
stage, but only after native strain evaluation. If its duplicate-point or
interface behavior is unsuitable, implement a small strain-specific point
locator/evaluator rather than smoothing or averaging the strain.

### Rendering

Create a strain-specific renderer, for example:

```text
numbat/plotting/strain_rendering.py
numbat/plotting/strain_modes.py
```

Reuse shared low-level contour, colourbar, trim, frame, title, and output
helpers where their abstractions fit. Do not add quiver logic or
displacement-specific overview routing to strain rendering.

The initial combined overview has exactly six field panels:

```text
Re[Sxx]  Re[Syy]  Re[Szz]
Im[Syz]  Im[Sxz]  Re[Sxy]
```

Each component has an independent signed colour scale by default. A future
revision may add a metadata panel, tensor invariant, or user-controlled
shared colour scales, but none is required for the initial feature.

### Export

Create strain-specific serializers alongside `io/mode_fields.py`, for
example `io/strain_fields.py`.

#### Two-dimensional output

Write separate real and imaginary arrays for every strain component:

```text
AC_strain_mode_00_Sxx_re.txt
AC_strain_mode_00_Sxx_im.txt
...
AC_strain_mode_00_Sxy_re.txt
AC_strain_mode_00_Sxy_im.txt
```

Headers must include mode index, strain storage key, and output-grid
coordinate summary.

#### One-dimensional output

Write one file per cut, with all six complex components:

```text
x Sxxr Sxxi Syyr Syyi Szzr Szzi Syzr Syzi Sxzr Sxzi Sxyr Sxyi
```

The cut coordinate and filenames should follow the established displacement
field conventions.

### Public API

Add acoustic-specific equivalents of the existing field-mode APIs:

```python
ModeAC.plot_strain(...)
ModeAC.plot_strain_1D(...)
ModeAC.write_strain(...)
ModeAC.write_strain_1D(...)

ACSimResult.plot_strain_modes(...)
ACSimResult.plot_strain_modes_1D(...)
ACSimResult.write_strains(...)
ACSimResult.write_strains_1D(...)
```

The public options should parallel existing mode-profile methods where
meaningful: mode selection, plot resolution, cuts, prefixes/suffixes, axes,
colourbars, contour levels, trim, frame, title, and decorators.

Explicitly reject strain methods for EM results. Do not silently substitute
displacement fields.

## Implementation Phases

### Phase 1: Selection and Native P2 Gradient Foundation

**Status: Complete**

1. Add `StrainComponent`, `StrainSelection`, storage keys, labels, parsing,
   ordering, and display-primary mapping.
2. Add P2 reference basis functions and reference derivatives in Python.
3. Add P2 isoparametric Jacobian and physical-gradient transformations.
4. Add extraction of acoustic element-local displacement and coordinates.
5. Implement computation of all six element-local complex strain components.

Validation:

- an affine-element displacement has exact transverse derivatives;
- `q_AC != 0` gives correct longitudinal terms;
- all components follow the agreed Voigt order;
- strain values have dimensionless units;
- EM inputs are rejected.

### Phase 2: Native Mesh Inspection and Interface Policy

**Status: Complete**

1. Add an internal native strain data object with element-local values and
   geometry.
2. Add a mesh diagnostic renderer using triangle faces or subtriangles.
3. Confirm visually that material-interface jumps are preserved.
4. Define and test deterministic one-sided behavior for points exactly on an
   interface.

Validation:

- adjacent elements with distinct analytical gradients do not get averaged;
- native rendering shows the expected discontinuity;
- no duplicate-node ambiguity silently changes field values.

### Phase 3: Standard 2D Sampling and Six-Panel Plot

Status: complete.

1. [x] Add rectangular-grid sampling from the native strain evaluator.
2. [x] Add a dedicated six-panel 2D overview and individual-component plot
   support.
3. [x] Apply the agreed display-primary mapping and signed colour maps.
4. [x] Integrate public 2D plotting APIs.

Validation:

- [x] Plot values agree with analytical strain at selected grid points.
- [x] Panel order and default real/imaginary parts are correct.
- [x] Trim, frame, title, colourbar, prefix, and explicit component options
  work.
- [x] EM calls fail clearly.

### Phase 4: 2D Export

Status: complete.

1. [x] Implement strain 2D text serialization.
2. [x] Add stable filenames, headers, and real/imaginary output for all six
   components.
3. [x] Document the schema.

Each mode writes twelve files named
``AC_strain_mode_{mode:02d}_S{ij}_{re|im}.txt``. Each file contains one
``(x, y)``-indexed sampled grid and a header with the mode index, physical
strain storage key (for example ``Sxzr``), and sampled ``v_x``/``v_y`` bounds.

Validation:

- [x] All twelve output files are produced per mode.
- [x] Headers and array shapes match the selected sample grid.
- [x] Exported values agree with sampled plot data.

### Phase 5: 1D Cuts, Plots, and Export

1. Add direct native strain sampling for `x`, `y`, and arbitrary line cuts.
2. Add 1D component plots.
3. Add one-file-per-cut strain export.
4. Add public `ModeAC` and `ACSimResult` façade methods.

Validation:

- all cut directions agree with analytical values;
- each exported file contains x plus twelve strain columns in the documented
  order;
- line cuts crossing interfaces follow the documented one-sided rule.

### Phase 6: Legacy Cleanup and Documentation

1. Replace the broken current `Mode.plot_strain` implementation with a
   compatibility wrapper to the new API, or deprecate/remove it according to
   the public API decision made during implementation.
2. Remove or leave archived the obsolete rectangular-gradient code without
   making it reachable.
3. Add user documentation and examples for 2D plots, 1D cuts, and exports.
4. Add standalone import coverage for new solver, plotting, I/O, and
   selection modules.

## Non-Goals for the Initial Feature

- Stress, energy density, acoustic power flow, and strain invariants.
- Averaged continuous strain at material interfaces.
- Changes to the Fortran/f2py output interface.
- Reworking the existing displacement field export schema.
- A generic rank-two tensor framework beyond what strain immediately needs.

## Completion Criteria

The feature is complete when:

1. all six strain components are computed from native P2 FEM data;
2. no rectangular-grid finite difference is used to calculate strain;
3. native interface discontinuities are retained;
4. 2D and 1D plots and text exports are available through acoustic public
   APIs;
5. complex-part defaults match the agreed acoustic phase convention;
6. analytical, interface, rendering, and export regression tests pass; and
7. the former broken strain hook no longer calls a nonexistent method.
