"""
Fresh test suite for ElasticIsotropicCylinder and related helpers.

This suite targets the current API where characteristic equations are
instance methods on ElasticIsotropicCylinder, and the radial wavenumber
helper is get_q_transverse.
"""

import sys
sys.path.append('backend')

import numpy as np
import pytest
from numpy.testing import assert_allclose

from nbanalytic.elasticcylinder import (
    ElasticIsotropicCylinder,
    get_q_transverse,
)

# Material constants (silicon-like) and geometry
RHO = 2329.0
C11 = 165.7e9
C12 = 63.9e9
C44 = 79.6e9
RAD = 500e-9  # 500 nm

TWOPI = 2 * np.pi


@pytest.fixture
def cylinder() -> ElasticIsotropicCylinder:
    return ElasticIsotropicCylinder(RHO, C11, C12, C44, RAD)


class TestGetQTransverse:
    def test_real_when_propagating(self):
        Om = 1e9 * TWOPI
        q = 1e6
        V = 5000.0
        qt = get_q_transverse(Om, q, V)
        assert isinstance(qt, complex)
        assert qt.imag == 0.0
        assert qt.real > 0.0
        assert_allclose(qt, np.sqrt((Om/V)**2 - q**2))

    def test_imag_when_evanescent(self):
        Om = 1e9 * TWOPI
        q = 2e6
        V = 5000.0
        qt = get_q_transverse(Om, q, V)
        assert qt.real == 0.0
        assert qt.imag > 0.0
        assert_allclose(abs(qt), np.sqrt(q**2 - (Om/V)**2))

    def test_zero_at_cutoff(self):
        q = 1e6
        V = 5000.0
        Om = q * V
        qt = get_q_transverse(Om, q, V)
        assert qt == 0j or qt == 0.0 or (qt.real == 0.0 and qt.imag == 0.0)


class TestCharacteristic:
    def test_dispatcher_routes(self, cylinder: ElasticIsotropicCylinder):
        Om = 5e9 * TWOPI
        q = 1e6

        # m = -1 -> torsional (p0 1D)
        res_m1 = cylinder.chareq_elastic_cylinder(Om, -1, q)
        res_1d = cylinder.chareq_elastic_cylinder_p0_1D(Om, 0, q)
        assert_allclose(res_m1, res_1d)

        # m = 0 -> Pochhammer 2D
        res_m0 = cylinder.chareq_elastic_cylinder(Om, 0, q)
        res_2d = cylinder.chareq_elastic_cylinder_p0_2D(Om, 0, q)
        assert_allclose(res_m0, res_2d)

        # m > 0 -> general
        res_m2 = cylinder.chareq_elastic_cylinder(Om, 2, q)
        res_pp = cylinder.chareq_elastic_cylinder_ppos(Om, 2, q)
        assert_allclose(res_m2, res_pp)

    def test_dispersion_relation_wrapper(self, cylinder: ElasticIsotropicCylinder):
        q = 1e6
        nu = 3e9
        m = 0
        direct = cylinder.chareq_elastic_cylinder(TWOPI * nu, m, q)
        wrapped = cylinder.dispersion_relation_at_q_nu(q, nu, m)
        assert_allclose(wrapped, direct)


class TestModeFinding:
    def test_find_modes_basic(self, cylinder: ElasticIsotropicCylinder):
        q = 1e6
        Om_hi = TWOPI * 6e9
        m = 0
        nmodes, om = cylinder.find_Omega_at_q(q, Om_hi, m, 2)
        assert nmodes >= 0
        assert len(om) == 2
        if nmodes > 1:
            assert np.all(np.diff(om[:nmodes]) >= 0)

    @pytest.mark.parametrize('m', [-1, 0, 1])
    def test_find_modes_various_m(self, cylinder: ElasticIsotropicCylinder, m):
        q = 1e6
        Om_hi = TWOPI * 8e9
        nmodes, om = cylinder.find_Omega_at_q(q, Om_hi, m, 1)
        assert nmodes >= 0
        if nmodes > 0:
            assert om[0] > 0

class TestSpecifcValues:
    """Tests verifying specific known values of the solutions."""

    def test_specific_values_m_minus_one(self):
        cylinder = ElasticIsotropicCylinder(RHO, C11, C12, C44, RAD)
        
        q = 1e6
        Om_hi = TWOPI * 10e9
        m = -1
        
        nmodes, omegas = cylinder.find_Omega_at_q(q, Om_hi, m, 5)
        Exp_Oms = [ 6.03314086e+10, 5.84679647e+09]

        assert_allclose(omegas[:2], Exp_Oms)

    def test_specific_values_m_zero(self):
        cylinder = ElasticIsotropicCylinder(RHO, C11, C12, C44, RAD)
        
        q = 1e6
        Om_hi = TWOPI * 10e9
        m = 0
        
        nmodes, omegas = cylinder.find_Omega_at_q(q, Om_hi, m, 5)
        Exp_Oms = [8.35805679e+09, 3.43346555e+10, 4.58949076e+10 ]

        assert_allclose(omegas[:3], Exp_Oms)
        
    def test_specific_values_m_one(self):
        cylinder = ElasticIsotropicCylinder(RHO, C11, C12, C44, RAD)
        
        q = 1e6
        Om_hi = TWOPI * 10e9
        m = 1
        
        nmodes, omegas = cylinder.find_Omega_at_q(q, Om_hi, m, 5)
        Exp_Oms = [3.95276474e+09, 2.34053712e+10, 3.20899647e+10, 5.99645381e+10]

        assert_allclose(omegas[:4], Exp_Oms)
        
    def test_specific_values_m_two(self):
        cylinder = ElasticIsotropicCylinder(RHO, C11, C12, C44, RAD)
        
        q = 1e6
        Om_hi = TWOPI * 10e9
        m = 2
        
        nmodes, omegas = cylinder.find_Omega_at_q(q, Om_hi, m, 5)
        Exp_Oms = [8.43483342e+09, 2.82663196e+10, 3.70607583e+10, 4.75638438e+10,5.85002648e+09]
        assert_allclose(omegas, Exp_Oms)


if __name__ == '__main__':
    import pytest as _pytest
    _pytest.main([__file__, '-q'])
