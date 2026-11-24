"""Native P2 finite-element sampling for acoustic strain fields."""

from dataclasses import dataclass

import numpy as np
from matplotlib.tri import Triangulation

from numbat.mesh.fem import validate_1D_cut_format
from numbat.nbtypes import SI_um
from numbat.strain_selection import (
    StrainComponent,
    component_index,
    storage_key,
)

P2_REFERENCE_NODES = np.array(
    (
        (0.0, 0.0),
        (1.0, 0.0),
        (0.0, 1.0),
        (0.5, 0.0),
        (0.5, 0.5),
        (0.0, 0.5),
    )
)

P2_SUBTRIANGLES = np.array((
    (0, 3, 5),
    (1, 4, 3),
    (2, 5, 4),
    (3, 4, 5),
))


@dataclass(frozen=True)
class NativeStrainSamples:
    """Element-local complex strain values evaluated on native P2 elements."""

    reference_points: np.ndarray
    coordinates: np.ndarray
    values: np.ndarray
    element_material_indices: np.ndarray | None = None

    def __post_init__(self):
        if self.coordinates.shape[:2] != self.values.shape[:2]:
            raise ValueError("Native strain coordinates and values must align.")
        if self.coordinates.shape[2] != 2 or self.values.shape[2] != 6:
            raise ValueError(
                "Native strain samples require (element, point, xy) coordinates "
                "and six Voigt values."
            )

    @property
    def n_elements(self):
        return self.values.shape[0]

    def subtriangles(self):
        """Return display triangles indexed into flattened element-local nodes."""
        if not np.array_equal(self.reference_points, P2_REFERENCE_NODES):
            raise ValueError(
                "P2 subtriangles require strain samples at P2 reference nodes."
            )
        offsets = 6 * np.arange(self.n_elements)[:, np.newaxis, np.newaxis]
        return (offsets + P2_SUBTRIANGLES).reshape(-1, 3)

    def flattened_coordinates(self):
        """Return local P2-node coordinates without averaging duplicate nodes."""
        return self.coordinates.reshape(-1, 2)

    def display_values(self, selection):
        """Return flattened real or imaginary values for one strain selection."""
        values = self.values[:, :, component_index(selection.component)]
        return (
            values.real.reshape(-1)
            if storage_key(selection).endswith("r")
            else values.imag.reshape(-1)
        )


def _reference_points_array(reference_points):
    points = np.asarray(reference_points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError("reference_points must have shape (n_points, 2).")
    if np.any(points < 0) or np.any(points.sum(axis=1) > 1):
        raise ValueError("reference_points must lie within the unit triangle.")
    return points


def p2_basis(reference_points):
    """Evaluate the six P2 Lagrange basis functions on the unit triangle."""
    r, s = _reference_points_array(reference_points).T
    t = 1 - r - s
    return np.column_stack((
        t * (2 * t - 1),
        r * (2 * r - 1),
        s * (2 * s - 1),
        4 * t * r,
        4 * r * s,
        4 * s * t,
    ))


def p2_basis_gradients(reference_points):
    """Evaluate P2 basis gradients with columns ordered as (d/dr, d/ds)."""
    r, s = _reference_points_array(reference_points).T
    gradients = np.empty((len(r), 6, 2), dtype=float)
    gradients[:, 0, :] = (-3 + 4 * r + 4 * s)[:, np.newaxis]
    gradients[:, 1, 0] = -1 + 4 * r
    gradients[:, 1, 1] = 0
    gradients[:, 2, 0] = 0
    gradients[:, 2, 1] = -1 + 4 * s
    gradients[:, 3, 0] = 4 - 8 * r - 4 * s
    gradients[:, 3, 1] = -4 * r
    gradients[:, 4, 0] = 4 * s
    gradients[:, 4, 1] = 4 * r
    gradients[:, 5, 0] = -4 * s
    gradients[:, 5, 1] = 4 - 4 * r - 8 * s
    return gradients


def calculate_element_strain(
    displacements, element_coordinates, q_ac, reference_points=P2_REFERENCE_NODES
):
    """Evaluate all six complex strain components on P2 mesh elements.

    Args:
        displacements: Complex values with shape
            ``(n_elements, 6, 3)`` ordered as ``(ux, uy, uz)``.
        element_coordinates: Physical x/y coordinates with shape
            ``(n_elements, 6, 2)``.
        q_ac: Acoustic longitudinal wavevector in ``1/m``.
        reference_points: P2 reference-triangle points with shape
            ``(n_points, 2)``.

    Returns:
        Complex values with shape ``(n_elements, n_points, 6)`` in Voigt
        order ``(xx, yy, zz, yz, xz, xy)``.
    """
    displacements = np.asarray(displacements, dtype=np.complex128)
    element_coordinates = np.asarray(element_coordinates, dtype=float)
    points = _reference_points_array(reference_points)
    if displacements.ndim != 3 or displacements.shape[1:] != (6, 3):
        raise ValueError("displacements must have shape (n_elements, 6, 3).")
    if element_coordinates.shape != (*displacements.shape[:2], 2):
        raise ValueError(
            "element_coordinates must have shape (n_elements, 6, 2)."
        )

    basis = p2_basis(points)
    gradients_ref = p2_basis_gradients(points)
    jacobians = np.einsum("eac,pad->epcd", element_coordinates, gradients_ref)
    determinants = np.linalg.det(jacobians)
    if np.any(np.abs(determinants) < np.finfo(float).tiny):
        raise ValueError("Cannot calculate strain on an element with a zero Jacobian.")

    gradients_physical = np.einsum(
        "epkd,epdc->epkc",
        np.einsum("eak,pad->epkd", displacements, gradients_ref),
        np.linalg.inv(jacobians),
    )
    displacement_values = np.einsum("pa,eak->epk", basis, displacements)

    strain = np.empty((*displacement_values.shape[:2], 6), dtype=np.complex128)
    ux, uy, uz = (displacement_values[:, :, index] for index in range(3))
    dux_dx, dux_dy = gradients_physical[:, :, 0, 0], gradients_physical[:, :, 0, 1]
    duy_dx, duy_dy = gradients_physical[:, :, 1, 0], gradients_physical[:, :, 1, 1]
    duz_dx, duz_dy = gradients_physical[:, :, 2, 0], gradients_physical[:, :, 2, 1]

    strain[:, :, component_index(StrainComponent.XX)] = dux_dx
    strain[:, :, component_index(StrainComponent.YY)] = duy_dy
    strain[:, :, component_index(StrainComponent.ZZ)] = 1j * q_ac * uz
    strain[:, :, component_index(StrainComponent.YZ)] = 0.5 * (
        duz_dy + 1j * q_ac * uy
    )
    strain[:, :, component_index(StrainComponent.XZ)] = 0.5 * (
        duz_dx + 1j * q_ac * ux
    )
    strain[:, :, component_index(StrainComponent.XY)] = 0.5 * (dux_dy + duy_dx)
    return strain


def _element_coordinates(fem_mesh):
    node_indices = fem_mesh.element_node_connectivity.T - 1
    return np.transpose(fem_mesh.node_coordinates[:, node_indices], (1, 2, 0))


def sample_native_strain(samples, points, tolerance=1e-12):
    """Sample native strain using deterministic one-sided subtriangle values.

    Points on an element boundary select the first matching subtriangle in
    element/local-triangle order. Points outside the native mesh return NaN.
    """
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError("points must have shape (n_points, 2).")

    coordinates = samples.flattened_coordinates()
    triangles = samples.subtriangles()
    triangle_coordinates = coordinates[triangles]
    triangle_values = samples.values.reshape(-1, 6)[triangles]
    sampled = np.full((len(points), 6), np.nan + 1j * np.nan, dtype=np.complex128)

    for point_index, point in enumerate(points):
        for triangle, values in zip(triangle_coordinates, triangle_values):
            origin = triangle[0]
            transform = np.column_stack((triangle[1] - origin, triangle[2] - origin))
            determinant = np.linalg.det(transform)
            if abs(determinant) < np.finfo(float).tiny:
                continue
            second, third = np.linalg.solve(transform, point - origin)
            first = 1 - second - third
            barycentric = np.array((first, second, third))
            if np.all(barycentric >= -tolerance):
                sampled[point_index] = barycentric @ values
                break

    return sampled


def strain_cut_points(sim_result, s_cut, val1, val2, n_points):
    """Return a displayed cut axis and matching native-FEM evaluation points."""
    validate_1D_cut_format(s_cut, val1, val2)
    if n_points < 2:
        raise ValueError("n_points must be at least 2 for a strain cut.")

    x_min, x_max, y_min, y_max = sim_result.fem_mesh.get_xy_limits()
    match s_cut.lower():
        case "x":
            y_points = np.linspace(y_min, y_max, n_points)
            x_points = val1 * SI_um + np.zeros(n_points)
        case "y":
            x_points = np.linspace(x_min, x_max, n_points)
            y_points = val1 * SI_um + np.zeros(n_points)
        case "line":
            x_points = np.linspace(val1[0] * SI_um, val2[0] * SI_um, n_points)
            y_points = np.linspace(val1[1] * SI_um, val2[1] * SI_um, n_points)

    if s_cut.lower() == "x":
        cut_axis = y_points
    elif s_cut.lower() == "y":
        cut_axis = x_points
    else:
        cut_axis = np.hypot(
            x_points - x_points[0], y_points - y_points[0]
        )

    return cut_axis / SI_um, np.column_stack((x_points, y_points))


class NativeStrainInterpolator:
    """Fast element-local strain interpolator for dense plot or cut samples."""

    def __init__(self, samples):
        self.samples = samples
        self.local_triangles = samples.subtriangles()
        local_coordinates = samples.flattened_coordinates()
        unique_coordinates, inverse = np.unique(
            local_coordinates, axis=0, return_inverse=True
        )
        self._triangulation = Triangulation(
            unique_coordinates[:, 0],
            unique_coordinates[:, 1],
            inverse[self.local_triangles],
        )
        self._finder = self._triangulation.get_trifinder()
        self._triangle_coordinates = local_coordinates[self.local_triangles]
        self._triangle_values = samples.values.reshape(-1, 6)[
            self.local_triangles
        ]

    def __call__(self, points):
        """Return native one-sided strain samples, using NaN outside the mesh."""
        points = np.asarray(points, dtype=float)
        if points.ndim != 2 or points.shape[1] != 2:
            raise ValueError("points must have shape (n_points, 2).")

        triangle_indices = self._finder(points[:, 0], points[:, 1])
        sampled = np.full(
            (len(points), 6), np.nan + 1j * np.nan, dtype=np.complex128
        )
        valid = triangle_indices >= 0
        if not np.any(valid):
            return sampled

        triangle_coordinates = self._triangle_coordinates[triangle_indices[valid]]
        origins = triangle_coordinates[:, 0]
        transforms = np.stack((
            triangle_coordinates[:, 1] - origins,
            triangle_coordinates[:, 2] - origins,
        ), axis=-1)
        second_third = np.linalg.solve(
            transforms, (points[valid] - origins)[..., np.newaxis]
        ).squeeze(-1)
        barycentric = np.column_stack((
            1 - second_third.sum(axis=1),
            second_third,
        ))
        sampled[valid] = np.einsum(
            "ni,nij->nj",
            barycentric,
            self._triangle_values[triangle_indices[valid]],
        )
        return sampled


def make_native_strain_interpolator(samples):
    """Create a cached fast interpolator for native strain samples."""
    return NativeStrainInterpolator(samples)


def calculate_mode_strain(sim_result, mode_index, reference_points=P2_REFERENCE_NODES):
    """Calculate element-local native strain for one acoustic mode."""
    if not sim_result.is_AC():
        raise ValueError("Acoustic strain is only available for acoustic simulations.")

    displacements = np.transpose(
        sim_result.fem_evecs[:, :6, mode_index, :], (2, 1, 0)
    )
    coordinates = _element_coordinates(sim_result.fem_mesh)
    points = _reference_points_array(reference_points)
    values = calculate_element_strain(
        displacements, coordinates, sim_result.q_AC, points
    )
    return NativeStrainSamples(
        reference_points=points,
        coordinates=np.einsum("pa,eac->epc", p2_basis(points), coordinates),
        values=values,
        element_material_indices=getattr(
            sim_result.fem_mesh, "v_el_2_mat_idx", None
        ),
    )


def sample_mode_strain_cut(
    sim_result, mode_index, s_cut, val1, val2=None, n_points=501
):
    """Sample a native one-sided acoustic strain cut for one mode."""
    cut_axis, points = strain_cut_points(
        sim_result, s_cut, val1, val2, n_points
    )
    samples = calculate_mode_strain(sim_result, mode_index)
    return cut_axis, sample_native_strain(samples, points)
