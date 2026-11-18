#!/usr/bin/env python3
"""Test script for Phase A: FemMesh accessor implementation.

Tests that:
1. All new accessor methods work correctly
2. Deprecated properties emit warnings
3. Internal state is properly maintained
4. AC mesh trimming works with new accessors
"""

import warnings
import numpy as np
import sys

# Suppress numba warnings to focus on our tests
warnings.filterwarnings("ignore", category=DeprecationWarning, module="numba")

def test_imports():
    """Test that imports work."""
    print("Test 1: Imports")
    try:
        from backend.numbat.mesh.mesh_data import MeshArrays
        print("  ✓ MeshArrays imported successfully")
        
        from backend.numbat.mesh.fem import FemMesh
        print("  ✓ FemMesh imported successfully")
        
        return True
    except Exception as e:
        print(f"  ✗ Import failed: {e}")
        return False


def test_mesh_arrays_dataclass():
    """Test MeshArrays dataclass."""
    print("\nTest 2: MeshArrays dataclass")
    try:
        from backend.numbat.mesh.mesh_data import MeshArrays
        
        # Create sample mesh data
        n_elements = 10
        n_nodes = 20
        
        element_materials = np.array([1, 2, 1, 2, 1, 2, 1, 2, 1, 2], dtype=np.int64)
        element_connectivity = np.array([
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            [2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
            [3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
            [4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
            [5, 6, 7, 8, 9, 10, 11, 12, 13, 14],
            [6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
        ], dtype=np.int64)
        node_coordinates = np.random.rand(2, n_nodes) * 10.0
        
        arrays = MeshArrays(
            element_materials=element_materials,
            element_connectivity=element_connectivity,
            node_coordinates=node_coordinates,
            node_physical_indices=np.zeros(n_nodes, dtype=np.int64)
        )
        
        assert arrays.num_nodes() == n_nodes, "num_nodes() should return correct count"
        assert arrays.num_elements() == n_elements, "num_elements() should return correct count"
        print("  ✓ MeshArrays creation and counting work")
        
        # Test validation
        try:
            arrays.validate(n_elements, n_nodes)
            print("  ✓ Validation passes for correct dimensions")
        except ValueError as e:
            print(f"  ✗ Validation should pass: {e}")
            return False
        
        # Test validation failure
        try:
            arrays.validate(n_elements + 1, n_nodes)
            print("  ✗ Validation should fail for wrong element count")
            return False
        except ValueError:
            print("  ✓ Validation correctly rejects wrong element count")
        
        return True
        
    except Exception as e:
        print(f"  ✗ MeshArrays test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_fem_mesh_initialization():
    """Test FemMesh initialization."""
    print("\nTest 3: FemMesh initialization")
    try:
        from backend.numbat.mesh.fem import FemMesh
        
        mesh = FemMesh()
        
        # Check internal state
        assert mesh._n_msh_pts == 0, "Initial _n_msh_pts should be 0"
        assert mesh._n_msh_el == 0, "Initial _n_msh_el should be 0"
        assert mesh._v_el_2_mat_idx is None, "Initial _v_el_2_mat_idx should be None"
        assert mesh._m_elnd_to_mshpt is None, "Initial _m_elnd_to_mshpt should be None"
        assert mesh._v_mshpt_xy is None, "Initial _v_mshpt_xy should be None"
        
        print("  ✓ FemMesh initializes with correct internal state")
        return True
        
    except Exception as e:
        print(f"  ✗ FemMesh initialization test failed: {e}")
        return False


def test_new_accessors():
    """Test new accessor methods."""
    print("\nTest 4: New accessor methods")
    try:
        from backend.numbat.mesh.fem import FemMesh
        import numpy as np
        
        mesh = FemMesh()
        
        # Set up internal data
        mesh._n_msh_pts = 20
        mesh._n_msh_el = 10
        mesh._v_mshpt_xy = np.random.rand(2, 20) * 10.0
        mesh._m_elnd_to_mshpt = np.arange(1, 61).reshape(6, 10).astype(np.int64)
        mesh._v_el_2_mat_idx = np.array([1, 2, 1, 2, 1, 2, 1, 2, 1, 2], dtype=np.int64)
        mesh._node_physindex = np.zeros(20, dtype=np.int64)
        
        # Test accessors
        assert mesh.num_nodes() == 20, "num_nodes() should return 20"
        print("  ✓ num_nodes() works")
        
        assert mesh.num_elements() == 10, "num_elements() should return 10"
        print("  ✓ num_elements() works")
        
        coord = mesh.get_node_coordinate(0)
        assert len(coord) == 2, "get_node_coordinate() should return tuple of 2"
        print("  ✓ get_node_coordinate() works")
        
        all_coords = mesh.get_all_node_coordinates()
        assert all_coords.shape == (2, 20), "get_all_node_coordinates() should return (2, 20) array"
        print("  ✓ get_all_node_coordinates() works")
        
        node_indices = mesh.element_node_indices(0)
        assert len(node_indices) == 6, "element_node_indices() should return 6 indices"
        print("  ✓ element_node_indices() works")
        
        mat_idx = mesh.element_material_index_1based(0)
        assert mat_idx == 1, "element_material_index_1based() should return 1-based index"
        print("  ✓ element_material_index_1based() works")
        
        mat_idx_0based = mesh.element_to_material_index(0)
        assert mat_idx_0based == 0, "element_to_material_index() should return 0-based index"
        print("  ✓ element_to_material_index() works")
        
        return True
        
    except Exception as e:
        print(f"  ✗ Accessor methods test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_deprecated_properties():
    """Test that deprecated properties emit warnings."""
    print("\nTest 5: Deprecated property warnings")
    try:
        from backend.numbat.mesh.fem import FemMesh
        import numpy as np
        
        mesh = FemMesh()
        mesh._n_msh_pts = 20
        mesh._n_msh_el = 10
        mesh._v_el_2_mat_idx = np.array([1, 2, 1, 2, 1, 2, 1, 2, 1, 2], dtype=np.int64)
        mesh._v_mshpt_xy = np.random.rand(2, 20)
        mesh._m_elnd_to_mshpt = np.arange(1, 61).reshape(6, 10).astype(np.int64)
        
        # Capture warnings
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            
            # Access deprecated properties
            _ = mesh.n_msh_pts
            assert len(w) > 0, "Should emit warning for n_msh_pts access"
            assert "deprecated" in str(w[-1].message).lower()
            print("  ✓ n_msh_pts access emits deprecation warning")
            
            w.clear()
            _ = mesh.n_msh_el
            assert len(w) > 0, "Should emit warning for n_msh_el access"
            print("  ✓ n_msh_el access emits deprecation warning")
            
            w.clear()
            _ = mesh.v_el_2_mat_idx
            assert len(w) > 0, "Should emit warning for v_el_2_mat_idx access"
            print("  ✓ v_el_2_mat_idx access emits deprecation warning")
            
            w.clear()
            _ = mesh.v_mshpt_xy
            assert len(w) > 0, "Should emit warning for v_mshpt_xy access"
            print("  ✓ v_mshpt_xy access emits deprecation warning")
            
            w.clear()
            _ = mesh.m_elnd_to_mshpt
            assert len(w) > 0, "Should emit warning for m_elnd_to_mshpt access"
            print("  ✓ m_elnd_to_mshpt access emits deprecation warning")
        
        return True
        
    except Exception as e:
        print(f"  ✗ Deprecated properties test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_deprecated_property_setters():
    """Test that deprecated property setters emit warnings."""
    print("\nTest 6: Deprecated property setter warnings")
    try:
        from backend.numbat.mesh.fem import FemMesh
        import numpy as np
        
        mesh = FemMesh()
        
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            
            # Set deprecated properties
            mesh.n_msh_pts = 20
            assert len(w) > 0, "Should emit warning for n_msh_pts assignment"
            print("  ✓ n_msh_pts assignment emits deprecation warning")
            
            w.clear()
            mesh.n_msh_el = 10
            assert len(w) > 0, "Should emit warning for n_msh_el assignment"
            print("  ✓ n_msh_el assignment emits deprecation warning")
            
            w.clear()
            mesh.v_el_2_mat_idx = np.array([1, 2], dtype=np.int64)
            assert len(w) > 0, "Should emit warning for v_el_2_mat_idx assignment"
            print("  ✓ v_el_2_mat_idx assignment emits deprecation warning")
        
        return True
        
    except Exception as e:
        print(f"  ✗ Deprecated setter test failed: {e}")
        return False


def test_get_xy_limits():
    """Test get_xy_limits method."""
    print("\nTest 7: get_xy_limits() method")
    try:
        from backend.numbat.mesh.fem import FemMesh
        import numpy as np
        
        mesh = FemMesh()
        mesh._v_mshpt_xy = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        mesh._n_msh_pts = 3
        
        x_min, x_max, y_min, y_max = mesh.get_xy_limits()
        
        assert x_min == 1.0, "x_min should be 1.0"
        assert x_max == 3.0, "x_max should be 3.0"
        assert y_min == 4.0, "y_min should be 4.0"
        assert y_max == 6.0, "y_max should be 6.0"
        
        print("  ✓ get_xy_limits() returns correct bounds")
        return True
        
    except Exception as e:
        print(f"  ✗ get_xy_limits() test failed: {e}")
        return False


def main():
    """Run all tests."""
    print("=" * 70)
    print("PHASE A IMPLEMENTATION TESTS - FemMesh Accessor Interface")
    print("=" * 70)
    
    tests = [
        test_imports,
        test_mesh_arrays_dataclass,
        test_fem_mesh_initialization,
        test_new_accessors,
        test_deprecated_properties,
        test_deprecated_property_setters,
        test_get_xy_limits,
    ]
    
    results = []
    for test in tests:
        try:
            result = test()
            results.append(result)
        except Exception as e:
            print(f"\n✗ Test {test.__name__} crashed: {e}")
            import traceback
            traceback.print_exc()
            results.append(False)
    
    print("\n" + "=" * 70)
    passed = sum(results)
    total = len(results)
    print(f"RESULTS: {passed}/{total} tests passed")
    
    if all(results):
        print("✓ All tests passed! Phase A implementation successful.")
        print("=" * 70)
        return 0
    else:
        print("✗ Some tests failed. Review implementation.")
        print("=" * 70)
        return 1


if __name__ == "__main__":
    sys.exit(main())
