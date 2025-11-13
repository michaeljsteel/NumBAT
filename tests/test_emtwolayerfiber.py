import sys
sys.path.append('backend')
import math

import numpy as np
import pytest

from nbanalytic.emtwolayerfiber import TwoLayerFiberEM
from nbanalytic.emconstants import EMPoln


def make_fiber(n_core=1.45, n_clad=1.44, arad=2.0e-6):
    return TwoLayerFiberEM(n_core=n_core, n_clad=n_clad, arad=arad)


def k_from_lambda(lam):
    return 2.0 * math.pi / lam


def test_vnumber_matches_formula():
    fib = make_fiber(arad=2.0e-6)
    lam = 1.55e-6
    k = k_from_lambda(lam)
    # Reference V-number
    V_ref = k * fib._rad_core * math.sqrt(fib._nco_sq - fib._ncl_sq)
    assert np.isclose(fib.Vnum_k(k), V_ref, rtol=0, atol=0.0)
    assert np.isclose(fib.Vnumb_lam(lam), V_ref, rtol=0, atol=0.0)


def test_he11_exists_and_within_bounds_singlemode():
    # Choose parameters giving V < 2.405 (single-mode regime)
    fib = make_fiber(arad=2.0e-6)
    lam = 1.55e-6
    k = k_from_lambda(lam)

    nmodes, v_neff = fib.find_neffs_for_k(k, EMPoln.HY, 1, 1, 3)
    assert nmodes >= 1
    neff = v_neff[0]
    # neff must lie strictly between cladding and core indices
    assert fib._n_clad < neff < fib._n_core


def test_hybrid_chareq_zero_at_root():
    fib = make_fiber(arad=2.0e-6)
    lam = 1.55e-6
    k = k_from_lambda(lam)

    nmodes, v_neff = fib.find_neffs_for_k(k, EMPoln.HY, 1, 1, 1)
    assert nmodes == 1
    neff = v_neff[0]
    # Characteristic equation should be ~0 at the found root
    val = fib.chareq_em_fib2_hy_m(neff, 1, k)
    assert abs(val) < 1e-6


@pytest.mark.parametrize(
    "arad, lam",
    [
        (2.0e-6, 1.55e-6),  # modest V, single-mode
        (4.0e-6, 1.55e-6),  # larger V, more modes possible
    ],
)
def test_all_found_neffs_within_bounds(arad, lam):
    fib = make_fiber(arad=arad)
    k = k_from_lambda(lam)
    nmodes, v_neff = fib.find_neffs_for_k(k, EMPoln.HY, 1, 3, 8)
    assert 0 <= nmodes <= 8
    if nmodes > 0:
        neffs = v_neff[:nmodes]
        assert np.all(neffs > fib._n_clad)
        assert np.all(neffs < fib._n_core)


def test_he11_tends_toward_ncore_for_large_V():
    # With larger radius (higher V), neff of HE11 approaches n_core
    fib = make_fiber(arad=10.0e-6)
    lam = 1.55e-6
    k = k_from_lambda(lam)

    neff = fib.find_neff_HE11_for_k(k)
    assert neff is not None
    # Expect HE11 to be close to n_core (loose tolerance)
    assert neff > fib._n_core - 1e-3
