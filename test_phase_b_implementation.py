"""
Test suite for Phase B: AC mesh trimming extraction.

Verifies that:
1. ACMeshTrimmer correctly extracts AC meshes from EM meshes
2. ACMeshTrimData stores trimming results properly
3. FemMesh.ac_build_from_em() uses ACMeshTrimmer correctly
4. Element conversion table is created and accessible
5. Backward compatibility is maintained
"""

import sys
from pathlib import Path

# Add backend to path
backend_path = Path(__file__).parent / "backend"
sys.path.insert(0, str(backend_path))

import numpy as np
import pytest

from numbat.mesh.fem import FemMesh
from numbat.mesh.ac_trimming import ACMeshTrimmer, ACMeshTrimData


# ============================================================================
# Mock Classes for Testing
# ============================================================================

class MockElasticProps:
    """Mock elastic properties for testing."""
    
    def __init__(self, elastic_material_indices=None):
        """Initialize with list of material indices that are elastic."""
        self.elastic_material_indices = elastic_material_indices or [2, 3]
        self.typ_el_AC = elastic_material_indices or [2, 3]
        self.n_mats_ac = len(self.elastic_material_indices)
        self.v_acoustic_mats = []  # Empty list for testing
    
    def is_elastic_material_index(self, mat_idx):
        """Return True if material index is elastic."""
        return mat_idx in self.elastic_material_indices
    
    def active_material_index(self, mat_idx):
        """Map original material index to active material index (1-based)."""
        # Simple mapping: elastic materials get sequential 1-based indices
        for i, idx in enumerate(sorted(self.elastic_material_indices)):
            if idx == mat_idx:
                return i + 1
        raise ValueError(f"Material {mat_idx} is not elastic")


class MockFemMesh:
    """Mock FEM mesh for testing."""
    
    def __init__(self, num_nodes=10, num_elements=4, elastic_material_indices=None):
        """Create a simple test mesh."""
        self.num_nodes = num_nodes
        self.num_elements = num_elements
        self.elastic_indices = elastic_material_indices or [2, 3]
        
        # Simple node coordinates arranged in a grid
        self._node_coordinates = np.zeros((2, num_nodes))
        for i in range(num_nodes):
            self._node_coordinates[0, i] = (i % 5) * 1.0
            self._node_coordinates[1, i] = (i // 5) * 1.0
        
        # Element connectivity: 6 nodes per element, 1-based indices
        self._element_connectivity = np.ones((6, num_elements), dtype=np.int64)
        for el in range(num_elements):
            base_idx = el * 2
            for node_offset in range(6):
                node_idx = (base_idx + node_offset) % num_nodes
                self._element_connectivity[node_offset, el] = node_idx + 1  # 1-based
        
        # Element materials: mix of elastic and non-elastic
        # Elements 0, 2 -> elastic (material 2, 3)
        # Elements 1, 3 -> vacuum (material 1)
        self._element_materials = np.array([2, 1, 3, 1], dtype=np.int64)
    
    def get_all_node_coordinates(self):
        return self._node_coordinates.copy()
    
    def get_all_element_connectivity(self):
        return self._element_connectivity.copy()
    
    def get_all_element_materials(self):
        return self._element_materials.copy()


# ============================================================================
# Test Cases
# ============================================================================

class TestACMeshTrimData:
    """Tests for ACMeshTrimData dataclass."""
    
    def test_trimdata_creation(self):
        """Verify ACMeshTrimData can be created with proper fields."""
        connectivity = np.array([[1, 2], [3, 4], [5, 6], [7, 8], [9, 10], [11, 12]])
        coordinates = np.array([[0, 1], [0.5, 1.5]])
        materials = [2, 3]
        el_conv = {0: 0, 1: 2}
        node_conv = {0: 0, 1: 1}
        
        trim_data = ACMeshTrimData(
            element_node_connectivity=connectivity,
            node_coordinates=coordinates,
            element_materials=materials,
            element_convert_table=el_conv,
            node_convert_table=node_conv,
        )
        
        assert trim_data.num_elements == 2
        assert trim_data.num_nodes == 2
        assert len(trim_data.element_materials) == 2
    
    def test_trimdata_properties(self):
        """Verify ACMeshTrimData properties work correctly."""
        connectivity = np.ones((6, 5))
        coordinates = np.ones((2, 3))
        
        trim_data = ACMeshTrimData(
            element_node_connectivity=connectivity,
            node_coordinates=coordinates,
            element_materials=[1, 2, 3, 4, 5],
            element_convert_table={},
            node_convert_table={},
        )
        
        assert trim_data.num_elements == 5
        assert trim_data.num_nodes == 3


class TestACMeshTrimmer:
    """Tests for ACMeshTrimmer utility class."""
    
    def test_trim_basic(self):
        """Test basic trimming of a simple mesh."""
        em_mesh = MockFemMesh(num_nodes=10, num_elements=4)
        elastic_props = MockElasticProps(elastic_material_indices=[2, 3])
        
        trim_data = ACMeshTrimmer.trim_em_mesh_to_ac(em_mesh, elastic_props)
        
        # Should have trimmed to 2 elastic elements (indices 0, 2)
        assert trim_data.num_elements == 2
        assert len(trim_data.element_materials) == 2
        assert len(trim_data.element_convert_table) == 2
        
        # Check conversion table maps correctly
        assert 0 in trim_data.element_convert_table
        assert 1 in trim_data.element_convert_table
        # AC element 0 should map to EM element 0 or 2
        assert trim_data.element_convert_table[0] in [0, 2]
        assert trim_data.element_convert_table[1] in [0, 2]
    
    def test_trim_preserves_node_uniqueness(self):
        """Verify trimming extracts unique nodes only."""
        em_mesh = MockFemMesh(num_nodes=10, num_elements=4)
        elastic_props = MockElasticProps(elastic_material_indices=[2, 3])
        
        trim_data = ACMeshTrimmer.trim_em_mesh_to_ac(em_mesh, elastic_props)
        
        # Node count should be less than 6 * num_elastic_elements due to sharing
        assert trim_data.num_nodes < 6 * trim_data.num_elements
        # Node count should be positive
        assert trim_data.num_nodes > 0
    
    def test_trim_no_elastic_raises_error(self):
        """Verify trimming raises error if no elastic elements."""
        em_mesh = MockFemMesh(num_nodes=10, num_elements=4)
        # Create elastic props with no overlapping materials
        elastic_props = MockElasticProps(elastic_material_indices=[99])
        
        with pytest.raises(ValueError, match="no elastic elements"):
            ACMeshTrimmer.trim_em_mesh_to_ac(em_mesh, elastic_props)
    
    def test_trim_connectivity_is_fortran_indexed(self):
        """Verify trimmed connectivity uses 1-based Fortran indexing."""
        em_mesh = MockFemMesh(num_nodes=10, num_elements=4)
        elastic_props = MockElasticProps(elastic_material_indices=[2, 3])
        
        trim_data = ACMeshTrimmer.trim_em_mesh_to_ac(em_mesh, elastic_props)
        
        # All node indices should be >= 1 (1-based Fortran)
        assert np.all(trim_data.element_node_connectivity >= 1)
        # All node indices should be <= num_nodes (valid range)
        assert np.all(trim_data.element_node_connectivity <= trim_data.num_nodes)
    
    def test_trim_coordinates_match_connectivity(self):
        """Verify trimmed coordinates are consistent with connectivity."""
        em_mesh = MockFemMesh(num_nodes=10, num_elements=4)
        elastic_props = MockElasticProps(elastic_material_indices=[2, 3])
        
        trim_data = ACMeshTrimmer.trim_em_mesh_to_ac(em_mesh, elastic_props)
        
        # Node coordinates shape should match num_nodes
        assert trim_data.node_coordinates.shape == (2, trim_data.num_nodes)
        # All coordinates should be finite
        assert np.all(np.isfinite(trim_data.node_coordinates))
    
    def test_trim_node_convert_table(self):
        """Verify node conversion table maps AC nodes to EM nodes."""
        em_mesh = MockFemMesh(num_nodes=10, num_elements=4)
        elastic_props = MockElasticProps(elastic_material_indices=[2, 3])
        
        trim_data = ACMeshTrimmer.trim_em_mesh_to_ac(em_mesh, elastic_props)
        
        # Node convert table should have entry for each AC node
        assert len(trim_data.node_convert_table) == trim_data.num_nodes
        # All values should be valid 0-based EM node indices
        for ac_idx, em_idx in trim_data.node_convert_table.items():
            assert 0 <= em_idx < em_mesh.num_nodes
            assert 0 <= ac_idx < trim_data.num_nodes


class MockWaveguide:
    """Reusable mock waveguide for FemMesh testing."""
    def __init__(self, elastic_props):
        self.elastic_props = elastic_props


class TestFemMeshACBuildFromEM:
    """Tests for FemMesh.ac_build_from_em() using ACMeshTrimmer."""
    
    def test_ac_build_from_em_basic(self):
        """Verify ac_build_from_em() populates FemMesh correctly."""
        em_mesh = MockFemMesh(num_nodes=10, num_elements=4)
        elastic_props = MockElasticProps(elastic_material_indices=[2, 3])
        
        ac_mesh = FemMesh()
        wg = MockWaveguide(elastic_props)
        ac_mesh.ac_build_from_em(wg, em_mesh)
        
        # Verify mesh properties are set
        assert ac_mesh.num_elements > 0
        assert ac_mesh.num_nodes > 0
        assert ac_mesh.num_elements <= 2  # Only 2 elastic elements
    
    def test_ac_build_sets_element_convert_table(self):
        """Verify element conversion table is set after trimming."""
        em_mesh = MockFemMesh(num_nodes=10, num_elements=4)
        elastic_props = MockElasticProps(elastic_material_indices=[2, 3])
        
        ac_mesh = FemMesh()
        wg = MockWaveguide(elastic_props)
        ac_mesh.ac_build_from_em(wg, em_mesh)
        
        # Element convert table should be populated
        convert_table = ac_mesh.get_element_convert_table()
        assert convert_table is not None
        assert len(convert_table) == ac_mesh.num_elements
        # All values should map to valid EM elements
        for ac_idx, em_idx in convert_table.items():
            assert 0 <= em_idx < em_mesh.num_elements
    
    def test_ac_build_maintains_coordinates(self):
        """Verify AC mesh coordinates are properly extracted."""
        em_mesh = MockFemMesh(num_nodes=10, num_elements=4)
        elastic_props = MockElasticProps(elastic_material_indices=[2, 3])
        
        ac_mesh = FemMesh()
        wg = MockWaveguide(elastic_props)
        ac_mesh.ac_build_from_em(wg, em_mesh)
        
        # Coordinates should be finite
        coords = ac_mesh.get_all_node_coordinates()
        assert np.all(np.isfinite(coords))
        # Should have correct shape
        assert coords.shape == (2, ac_mesh.num_nodes)
    
    def test_ac_build_deprecated_property_access(self):
        """Verify backward compatibility with deprecated el_convert_tbl property."""
        em_mesh = MockFemMesh(num_nodes=10, num_elements=4)
        elastic_props = MockElasticProps(elastic_material_indices=[2, 3])
        
        ac_mesh = FemMesh()
        wg = MockWaveguide(elastic_props)
        ac_mesh.ac_build_from_em(wg, em_mesh)
        
        # Old property should still work (with deprecation warning)
        import warnings
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            convert_tbl = ac_mesh.el_convert_tbl
            # Should have issued a deprecation warning
            assert len(w) > 0
            assert issubclass(w[-1].category, DeprecationWarning)
        
        # Content should match new accessor
        assert convert_tbl == ac_mesh.get_element_convert_table()


class TestPhaseABackwardCompatibility:
    """Verify Phase A accessors still work after Phase B changes."""
    
    def test_num_elements_works(self):
        """Verify num_elements() accessor works."""
        em_mesh = MockFemMesh(num_nodes=10, num_elements=4)
        elastic_props = MockElasticProps(elastic_material_indices=[2, 3])
        
        ac_mesh = FemMesh()
        wg = MockWaveguide(elastic_props)
        ac_mesh.ac_build_from_em(wg, em_mesh)
        
        # Should return positive integer
        n_el = ac_mesh.num_elements
        assert isinstance(n_el, int)
        assert n_el > 0
    
    def test_num_nodes_works(self):
        """Verify num_nodes() accessor works."""
        em_mesh = MockFemMesh(num_nodes=10, num_elements=4)
        elastic_props = MockElasticProps(elastic_material_indices=[2, 3])
        
        ac_mesh = FemMesh()
        wg = MockWaveguide(elastic_props)
        ac_mesh.ac_build_from_em(wg, em_mesh)
        
        # Should return positive integer
        n_pts = ac_mesh.num_nodes
        assert isinstance(n_pts, int)
        assert n_pts > 0
    
    def test_get_all_node_coordinates_works(self):
        """Verify get_all_node_coordinates() works after trimming."""
        em_mesh = MockFemMesh(num_nodes=10, num_elements=4)
        elastic_props = MockElasticProps(elastic_material_indices=[2, 3])
        
        ac_mesh = FemMesh()
        wg = MockWaveguide(elastic_props)
        ac_mesh.ac_build_from_em(wg, em_mesh)
        
        coords = ac_mesh.get_all_node_coordinates()
        assert coords.shape == (2, ac_mesh.num_nodes)
        assert np.all(np.isfinite(coords))
    
    def test_get_all_element_connectivity_works(self):
        """Verify get_all_element_connectivity() works after trimming."""
        em_mesh = MockFemMesh(num_nodes=10, num_elements=4)
        elastic_props = MockElasticProps(elastic_material_indices=[2, 3])
        
        ac_mesh = FemMesh()
        wg = MockWaveguide(elastic_props)
        ac_mesh.ac_build_from_em(wg, em_mesh)
        
        connectivity = ac_mesh.get_all_element_connectivity()
        assert connectivity.shape == (6, ac_mesh.num_elements)
        # Should be Fortran-indexed (>= 1)
        assert np.all(connectivity >= 1)
        # Should be valid node indices
        assert np.all(connectivity <= ac_mesh.num_nodes)


# ============================================================================
# Run Tests
# ============================================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
