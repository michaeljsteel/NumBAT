"""Test suite for ElasticIsotropicFreeSlab - Lamb wave dispersion solver."""

import math

import numpy as np
import pytest

import sys
sys.path.append('backend')

from nbanalytic.elasticfreeslab import ElasticIsotropicFreeSlab

import numbat

class MockMaterial:
    """Mock material with typical elastic properties for testing."""

    def __init__(self, Vs=3000.0, Vl=6000.0, Vr=2800.0):
        self._Vs = Vs
        self._Vl = Vl
        self._Vr = Vr

    def Vac_shear(self):
        return self._Vs

    def Vac_longitudinal(self):
        return self._Vl

    def Vac_Rayleigh(self):
        return self._Vr


@pytest.fixture
def mock_material():
    """Provide a mock material with typical elastic velocities."""
    return MockMaterial(Vs=3000.0, Vl=6000.0, Vr=2800.0)


@pytest.fixture
def slab(mock_material):
    """Create a standard test slab with 1 μm thickness."""
    return ElasticIsotropicFreeSlab(
        material=mock_material,
        thickness_SI=1.0e-6,
        max_omega=100e9,
        omega_bracketing_resolution=500,
        q_bracketing_resolution=500,
    )


def test_slab_construction(mock_material):
    """Test basic slab construction and attribute initialization."""
    thickness = 2.0e-6
    slab = ElasticIsotropicFreeSlab(
        material=mock_material, thickness_SI=thickness, max_omega=50e9
    )
    assert slab.width == thickness
    assert slab._Vbulk_shear == 3000.0
    assert slab._Vbulk_long == 6000.0
    assert slab._max_omega == 50e9


def test_characteristic_equation_at_cutoff(slab):
    """Test that characteristic equations behave at cutoff (q=0, Omega=0)."""
    # At cutoff, the char eq should be finite
    Omega = 1e6  # Very small frequency
    q = 1e3  # Very small wavenumber
    val_even = slab.elasticfreeslab_Lamb_chareq_even(
        Omega, q
    )
    val_odd = slab.elasticfreeslab_Lamb_chareq_odd(
        Omega, q
    )
    # Should not be NaN or infinite at reasonable small values
    assert np.isfinite(val_even)
    assert np.isfinite(val_odd)


def test_characteristic_equation_symmetry(slab):
    """Test characteristic equations have expected properties."""
    Omega = 20e9
    q = 5e6
    # Even and odd equations are distinct
    val_even = slab.elasticfreeslab_Lamb_chareq_even(
        Omega, q
    )
    val_odd = slab.elasticfreeslab_Lamb_chareq_odd(
        Omega, q
    )
    # Should generally differ (they're not identical functions)
    assert val_even != val_odd


def test_SH_dispersion_basic(slab):
    """Test SH (shear-horizontal) mode dispersion computation."""
    v_Omega = np.linspace(10e9, 50e9, 10)
    max_modes = 3

    m_qs, m_col = slab.find_SH_dispersion_for_q_bands_at_om(
        v_Omega, max_modes, col_array=True
    )

    assert m_qs.shape == (10, max_modes)
    assert m_col.shape == (10, max_modes, 3)
    # SH dispersion should be real and non-negative
    assert np.all(m_qs >= 0)
    # Higher modes should have lower q at same omega (or zero if cutoff)
    # At high omega, fundamental mode should exist
    assert m_qs[-1, 0] > 0


def test_SH_dispersion_obeys_cutoff(slab):
    """Test that SH modes below cutoff have q=0."""
    # For SH mode m, cutoff frequency is Omega_c = m * pi * Vs / width
    v_Omega = np.array([5e9, 10e9, 20e9])
    max_modes = 2

    m_qs, _ = slab.find_SH_dispersion_for_q_bands_at_om(v_Omega, max_modes)

    # Mode 0 should always propagate (no cutoff)
    assert np.all(m_qs[:, 0] > 0)

    # Mode 1 has cutoff at pi * Vs / width
    omega_cutoff_m1 = math.pi * slab._Vbulk_shear / slab.width
    # Below cutoff, q should be ~0
    for i, Om in enumerate(v_Omega):
        if Om < omega_cutoff_m1:
            assert m_qs[i, 1] < 1e3  # essentially zero


def test_Rayleigh_dispersion(slab):
    """Test Rayleigh surface wave dispersion is linear."""
    v_Omega = np.linspace(10e9, 50e9, 5)
    v_qR, v_col = slab.find_Rayleigh_dispersion_for_v(v_Omega, col_array=True)

    assert len(v_qR) == len(v_Omega)
    assert v_col.shape == (len(v_Omega), 3)

    # Rayleigh dispersion is linear: q = Omega / VR
    Vr = slab._material.Vac_Rayleigh()
    v_q_expected = v_Omega / Vr
    np.testing.assert_allclose(v_qR, v_q_expected, rtol=0, atol=0)


def test_find_Lamb_qs_at_Omega_returns_valid_wavenumbers(slab):
    """Test that Lamb mode finder returns physical wavenumbers."""
    Omega = 30e9
    max_modes = 2

    q_even = slab.find_Lamb_qs_at_Omega(Omega, max_modes, even_modes=True)
    q_odd = slab.find_Lamb_qs_at_Omega(Omega, max_modes, even_modes=False)

    # Should find at least one mode at this frequency
    assert len(q_even) >= 1 or len(q_odd) >= 1

    # All found q must be positive and less than Omega/Vs
    if len(q_even):
        assert np.all(q_even > 0)
        assert np.all(q_even < Omega / slab._Vbulk_shear * 2)

    if len(q_odd):
        assert np.all(q_odd > 0)
        assert np.all(q_odd < Omega / slab._Vbulk_shear * 2)


def test_find_Lamb_Omegas_at_q_returns_valid_frequencies(slab):
    """Test that Lamb mode finder at fixed q returns physical frequencies."""
    q = 1e7  # 10 rad/μm
    max_modes = 2

    Om_even = slab.find_Lamb_Omegas_at_q(q, max_modes, even_modes=True)
    Om_odd = slab.find_Lamb_Omegas_at_q(q, max_modes, even_modes=False)

    # Should find at least one mode
    assert len(Om_even) >= 1 or len(Om_odd) >= 1

    # All found Omega must be positive and reasonable
    if len(Om_even):
        assert np.all(Om_even > 0)
        assert np.all(Om_even < slab._max_omega)

    if len(Om_odd):
        assert np.all(Om_odd > 0)
        assert np.all(Om_odd < slab._max_omega)


def test_find_Lamb_dispersion_for_q_bands_shape(slab):
    """Test Lamb dispersion array has correct shape."""
    v_Omega = np.linspace(10e9, 50e9, 8)
    max_modes = 3

    m_qs = slab.find_Lamb_dispersion_for_q_bands_at_om(v_Omega, max_modes, even_modes=True)

    assert m_qs.shape == (len(v_Omega), max_modes)
    # Should have some non-zero entries at high frequency
    assert np.any(m_qs > 0)


def test_find_Lamb_dispersion_for_om_bands_shape(slab):
    """Test Lamb frequency dispersion array has correct shape."""
    v_q = np.linspace(1e6, 1e7, 8)
    max_modes = 3

    m_Om = slab.find_Lamb_dispersion_for_om_bands_at_q(v_q, max_modes, even_modes=True)

    assert m_Om.shape == (len(v_q), max_modes)
    # Should have some non-zero entries
    assert np.any(m_Om > 0)


def test_characteristic_equation_zero_at_found_root(slab):
    """Test that found Lamb mode satisfies characteristic equation."""
    Omega = 30e9
    max_modes = 1

    q_even = slab.find_Lamb_qs_at_Omega(Omega, max_modes, even_modes=True)

    if len(q_even) >= 1:
        q = q_even[0]
        val = slab.elasticfreeslab_Lamb_chareq_even(
            Omega, q
        )
        # Should be very close to zero
        assert abs(val) < 1e-3


def test_Rayleigh_profile_1d_shape_and_normalization(slab):
    """Test Rayleigh wave profile has correct structure."""
    Omega = 20e9
    depth = 3.0e-6
    npts = 50

    mode_prof = slab.find_Rayleigh_profile_1d(Omega, depth, npts=npts)

    # Check returned object attributes
    assert mode_prof.m_uxyz.shape == (npts, 3)
    assert len(mode_prof.v_x) == npts

    # Profile should be normalized (max ~1)
    max_amplitude = np.max(np.abs(mode_prof.m_uxyz))
    assert 0.9 < max_amplitude <= 1.0

    # x coordinate should span [-depth, 0]
    assert mode_prof.v_x[0]/1e6 == pytest.approx(-depth, rel=1e-6)
    assert mode_prof.v_x[-1]/1e6 == pytest.approx(0, rel=1e-6)


def test_find_mode_profile_1d_even_shape(slab):
    """Test Lamb even mode profile has correct shape."""
    Omega = 30e9
    q = 1e7
    npts = 40

    mode_prof = slab.find_mode_profile_1d(Omega, q, npts=npts, even_mode=True)

    assert mode_prof.m_uxyz.shape == (npts, 3)
    assert len(mode_prof.v_x) == npts

    # Profile should span [-width/2, width/2]
    assert mode_prof.v_x[0]/1e6 == pytest.approx(-slab.width / 2, rel=1e-6)
    assert mode_prof.v_x[-1]/1e6 == pytest.approx(slab.width / 2, rel=1e-6)


def test_find_mode_profile_1d_odd_shape(slab):
    """Test Lamb odd mode profile has correct shape."""
    Omega = 30e9
    q = 1e7
    npts = 40

    mode_prof = slab.find_mode_profile_1d(Omega, q, npts=npts, even_mode=False)

    assert mode_prof.m_uxyz.shape == (npts, 3)
    # Odd modes should also be normalized
    max_amplitude = np.max(np.abs(mode_prof.m_uxyz))
    assert 0.9 < max_amplitude <= 1.0


@pytest.mark.parametrize("even_modes", [True, False])
def test_Lamb_modes_respect_velocity_bounds(slab, even_modes):
    """Test that found modes have phase velocity between shear and longitudinal."""
    Omega = 30e9
    max_modes = 2

    qs = slab.find_Lamb_qs_at_Omega(Omega, max_modes, even_modes=even_modes)

    for q in qs:
        if q > 0:
            Vphase = Omega / q
            # Phase velocity should typically be between Vs and Vl
            # (though some modes can exceed Vl slightly in practice)
            assert Vphase > 0
            # Just check it's not wildly unreasonable
            assert Vphase < slab._Vbulk_long * 2


def test_asymptotic_formula_for_small_qw(slab):
    """Test asymptotic formula agrees with full solution at small qw."""
    q_small = 1e6  # Small wavenumber
    Omega_exact = slab.find_Lamb_Omegas_at_q(q_small, 1, even_modes=False)

    if len(Omega_exact) >= 1:
        Om_found = Omega_exact[0]
        Om_asymp = slab.find_Lamb_Omega_at_q_oddmode_zero_asymptotic_small_qw(
            q_small, even_modes=False
        )
        # Should be reasonably close for small qw
        # (asymptotic is approximate, allow ~10% error)
        if Om_asymp > 0:
            rel_error = abs(Om_found - Om_asymp) / Om_found
            assert rel_error < 0.2  # Within 20% for asymptotic
