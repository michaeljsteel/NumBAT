"""Native P2 finite-element sampling for acoustic strain fields."""

from dataclasses import dataclass

import numpy as np

from numbat.strain_selection import StrainComponent

_VOIGT_INDEX = {
    component: index for index, component in enumerate(StrainComponent)
}


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


@dataclass(frozen=True)
class NativeStrainSamples:
    """Element-local complex strain values evaluated on native P2 elements."""

    reference_points: np.ndarray
    coordinates: np.ndarray
    values: np.ndarray


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

    strain[:, :, _VOIGT_INDEX[StrainComponent.XX]] = dux_dx
    strain[:, :, _VOIGT_INDEX[StrainComponent.YY]] = duy_dy
    strain[:, :, _VOIGT_INDEX[StrainComponent.ZZ]] = 1j * q_ac * uz
    strain[:, :, _VOIGT_INDEX[StrainComponent.YZ]] = 0.5 * (
        duz_dy + 1j * q_ac * uy
    )
    strain[:, :, _VOIGT_INDEX[StrainComponent.XZ]] = 0.5 * (
        duz_dx + 1j * q_ac * ux
    )
    strain[:, :, _VOIGT_INDEX[StrainComponent.XY]] = 0.5 * (dux_dy + duy_dx)
    return strain


def _element_coordinates(fem_mesh):
    node_indices = fem_mesh.m_elnd_to_mshpt.T - 1
    return np.transpose(fem_mesh.v_mshpt_xy[:, node_indices], (1, 2, 0))


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
    )
