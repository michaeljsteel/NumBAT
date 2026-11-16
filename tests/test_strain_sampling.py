"""Analytical regression tests for native acoustic P2 strain sampling."""

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

sys.path.append(str(Path(__file__).resolve().parents[1] / "backend"))

from numbat.solvers.strain_sampling import (
    P2_REFERENCE_NODES,
    calculate_element_strain,
    calculate_mode_strain,
)
from numbat.strain_selection import (
    StrainComponent,
    StrainSelection,
    complex_storage_keys,
    selection_from_storage_key,
    storage_key,
)


def test_strain_selection_primary_storage_matches_acoustic_phase_convention():
    expected_primary_keys = {
        StrainComponent.XX: "Sxxr",
        StrainComponent.YY: "Syyr",
        StrainComponent.ZZ: "Szzr",
        StrainComponent.YZ: "Syzi",
        StrainComponent.XZ: "Sxzi",
        StrainComponent.XY: "Sxyr",
    }

    for component, expected_key in expected_primary_keys.items():
        selection = StrainSelection.primary(component)
        assert storage_key(selection) == expected_key
        assert selection_from_storage_key(expected_key) == selection
        assert complex_storage_keys(component) == (
            f"S{component.value}r",
            f"S{component.value}i",
        )


def test_native_p2_strain_matches_an_affine_displacement_on_curvilinear_element():
    coordinates = np.array((
        (
            (0.0, 0.0),
            (2.0, 0.0),
            (0.0, 1.5),
            (1.1, -0.1),
            (1.1, 0.8),
            (-0.1, 0.8),
        ),
    ))
    q_ac = 3.0
    gradients = np.array((
        (2.0 + 1.0j, -1.0 + 2.0j),
        (4.0 - 3.0j, 5.0 + 1.0j),
        (-2.0 + 4.0j, 3.0 - 2.0j),
    ))
    offsets = np.array((1.0 - 1.0j, -2.0 + 3.0j, 4.0 + 2.0j))
    displacements = offsets + np.einsum("eac,kc->eak", coordinates, gradients)

    strain = calculate_element_strain(displacements, coordinates, q_ac)
    points = np.einsum(
        "pa,eac->epc",
        np.array((
            (1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            (0.0, 1.0, 0.0, 0.0, 0.0, 0.0),
            (0.0, 0.0, 1.0, 0.0, 0.0, 0.0),
            (0.0, 0.0, 0.0, 1.0, 0.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 1.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, 1.0),
        )),
        coordinates,
    )
    expected = np.empty((1, len(P2_REFERENCE_NODES), 6), dtype=np.complex128)
    expected[:, :, 0] = gradients[0, 0]
    expected[:, :, 1] = gradients[1, 1]
    expected[:, :, 2] = 1j * q_ac * (
        offsets[2] + np.einsum("epc,c->ep", points, gradients[2])
    )
    expected[:, :, 3] = 0.5 * (
        gradients[2, 1]
        + 1j * q_ac * (offsets[1] + np.einsum("epc,c->ep", points, gradients[1]))
    )
    expected[:, :, 4] = 0.5 * (
        gradients[2, 0]
        + 1j * q_ac * (offsets[0] + np.einsum("epc,c->ep", points, gradients[0]))
    )
    expected[:, :, 5] = 0.5 * (gradients[0, 1] + gradients[1, 0])
    np.testing.assert_allclose(strain, expected)


def test_mode_strain_rejects_non_acoustic_simulations():
    sim_result = SimpleNamespace(is_AC=lambda: False)

    with pytest.raises(ValueError, match="Acoustic strain"):
        calculate_mode_strain(sim_result, 0)


def test_mode_strain_uses_acoustic_fem_storage_and_element_coordinates():
    coordinates = np.array((
        (0.0, 0.0),
        (1.0, 0.0),
        (0.0, 1.0),
        (0.5, 0.0),
        (0.5, 0.5),
        (0.0, 0.5),
    ))
    displacements = np.array(
        [[(2 * x + 3 * y, 5 * x + 7 * y, 11 * x + 13 * y) for x, y in coordinates]],
        dtype=np.complex128,
    )
    fem_evecs = np.zeros((3, 6, 1, 1), dtype=np.complex128)
    fem_evecs[:, :, 0, 0] = displacements[0].T
    fem_mesh = SimpleNamespace(
        m_elnd_to_mshpt=np.arange(1, 7, dtype=np.int64)[:, np.newaxis],
        v_mshpt_xy=coordinates.T,
    )
    sim_result = SimpleNamespace(
        is_AC=lambda: True,
        fem_evecs=fem_evecs,
        fem_mesh=fem_mesh,
        q_AC=0.0,
    )

    samples = calculate_mode_strain(sim_result, mode_index=0)

    expected = calculate_element_strain(displacements, coordinates[np.newaxis], 0.0)
    np.testing.assert_allclose(samples.values, expected)
    np.testing.assert_allclose(
        samples.coordinates, coordinates[np.newaxis]
    )
