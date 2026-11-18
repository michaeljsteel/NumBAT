# Copyright (C) 2017-2025  Michael Steel, Bjorn Sturmberg, Kokou Dossou.

# NumBAT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

"""Data structures for FEM mesh metadata and arrays."""

from dataclasses import dataclass
from typing import Optional

import numpy as np


@dataclass
class MeshArrays:
    """Arrays returned from Fortran simulation or loaded from gmsh file.
    
    All arrays use Fortran indexing conventions where appropriate:
    - Element indices: 1-based
    - Node indices: 1-based
    - Material indices: 1-based
    
    Fortran equivalents:
    - MeshEM type: defined in meshprops.f90
    - MeshAC type: defined in meshprops.f90
    """
    
    element_materials: np.ndarray
    """Material index for each element. Shape (n_elements,), 1-based indices.
    
    Fortran equivalent:
    - MeshEM%v_elt_material or MeshAC%v_elt_material
    - type_el in calc_modes_em_impl.f90
    - Ranges: [1..n_mat_el] for EM, [1..n_mat_ac_el] for AC
    """
    
    element_connectivity: np.ndarray
    """Element-to-node mapping. Shape (6, n_elements), 1-based node indices (Fortran convention).
    
    Fortran equivalent:
    - MeshEM%m_elnd_to_mshpt or MeshAC%m_elnd_to_mshpt
    - m_elnd_to_mshpt in calc_modes_em_impl.f90
    - [i, j] = ith node (1..6) of element j
    - Nodes 1-3 are vertices, 4-6 are edge midpoints
    - Ranges: [1..n_msh_pts]
    """
    
    node_coordinates: np.ndarray
    """Physical (x, y) coordinates of all nodes. Shape (2, n_nodes), physical units.
    
    Fortran equivalent:
    - MeshEM%v_mshpt_xy or MeshAC%v_mshpt_xy
    - v_mshpt_xy in calc_modes_em_impl.f90
    - [:, i] = (x, y) coordinate of node i in physical units
    """
    
    node_physical_indices: Optional[np.ndarray] = None
    """Physical boundary line/surface indices for each node. Shape (n_nodes,).
    
    Used for EM meshes to identify boundary nodes; not typically set for AC meshes (which use Neumann BCs).
    If on a Gmsh Physical Line, contains the index of that line (1-based). Otherwise 0 (interior point).
    
    Fortran equivalent:
    - MeshEM%v_mshpt_physindex or MeshAC%v_mshpt_physindex
    - node_physindex in calc_modes_em_impl.f90
    - Ranges: [0..max_phys_index]
    """
    
    def num_nodes(self) -> int:
        """Return total number of nodes."""
        return self.node_coordinates.shape[1]
    
    def num_elements(self) -> int:
        """Return total number of elements."""
        return self.element_connectivity.shape[1]
    
    def validate(self, expected_num_elements: int, expected_num_nodes: int) -> None:
        """Validate array dimensions match expectations.
        
        Args:
            expected_num_elements: Expected number of elements.
            expected_num_nodes: Expected number of nodes.
            
        Raises:
            ValueError: If array dimensions don't match expectations.
        """
        actual_elements = self.num_elements()
        actual_nodes = self.num_nodes()
        
        if actual_elements != expected_num_elements:
            raise ValueError(
                f"Element count mismatch: expected {expected_num_elements}, got {actual_elements}"
            )
        if actual_nodes != expected_num_nodes:
            raise ValueError(
                f"Node count mismatch: expected {expected_num_nodes}, got {actual_nodes}"
            )
        if self.element_materials.shape[0] != expected_num_elements:
            raise ValueError(
                f"Material array length mismatch: expected {expected_num_elements}, "
                f"got {self.element_materials.shape[0]}"
            )
        if self.node_physical_indices is not None:
            if self.node_physical_indices.shape[0] != expected_num_nodes:
                raise ValueError(
                    f"Physical index array length mismatch: expected {expected_num_nodes}, "
                    f"got {self.node_physical_indices.shape[0]}"
                )
