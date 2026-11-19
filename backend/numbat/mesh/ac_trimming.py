"""AC mesh trimming utilities for extracting acoustic elements from EM meshes.

This module provides standalone utilities to trim an EM mesh down to only the
acoustic (elastic) regions, mapping element and node indices appropriately.

The trimming process:
1. Filters the EM mesh to keep only elastic material elements
2. Extracts unique nodes used by those elastic elements
3. Remaps element connectivity to use the new compressed node indices
4. Returns mapping tables for translating between EM and AC grids
"""

from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass
class ACMeshTrimData:
    """Result of trimming an EM mesh to acoustic (elastic) regions.
    
    Attributes:
        element_node_connectivity: (6, n_msh_el_AC) array, Fortran-indexed connectivity
            Maps each AC element to its 6 nodes (1-based indices into AC node array)
        node_coordinates: (2, n_msh_pts_AC) array of x, y coordinates
        element_materials: list of material indices (1-based) for each AC element
        element_convert_table: dict mapping AC element indices to original EM element indices
            Used to translate between trimmed AC mesh and full EM mesh for field extraction
        node_convert_table: dict mapping AC node indices (0-based) to original EM node indices (0-based)
            Used to translate AC node indices back to full EM mesh node indices
    """
    element_node_connectivity: np.ndarray
    node_coordinates: np.ndarray
    element_materials: list
    element_convert_table: dict
    node_convert_table: dict
    
    @property
    def num_elements(self) -> int:
        """Return number of acoustic elements."""
        return self.element_node_connectivity.shape[1]
    
    @property
    def num_nodes(self) -> int:
        """Return number of acoustic nodes."""
        return self.node_coordinates.shape[1]


class ACMeshTrimmer:
    """Utility class for trimming EM meshes to acoustic (elastic) regions."""
    
    @staticmethod
    def trim_em_mesh_to_ac(
        em_fem_mesh: Any, 
        elastic_props: Any
    ) -> ACMeshTrimData:
        """Extract acoustic mesh from an EM mesh by filtering to elastic materials.
        
        This method extracts only the elements with elastic material properties,
        remaps node indices, and returns the trimmed mesh data with conversion tables
        for translating between the full EM mesh and the trimmed AC mesh.
        
        Args:
            em_fem_mesh: FemMesh instance containing the full EM mesh
                Must have methods: num_elements(), get_all_element_connectivity(),
                get_all_node_coordinates(), get_all_element_materials()
            elastic_props: ElasticProps or similar interface providing:
                - is_elastic_material_index(mat_idx: int) -> bool
                - active_material_index(mat_idx: int) -> int (1-based material index)
        
        Returns:
            ACMeshTrimData: Trimmed mesh data with element/node mappings
        
        Raises:
            ValueError: If trimming results in no acoustic elements
        """
        # Extract full EM mesh data using accessors
        # Fortran: n_msh_el, m_elnd_to_mshpt, v_mshpt_xy, type_el
        num_elements_em = em_fem_mesh.num_elements
        element_node_connectivity_em = em_fem_mesh.get_all_element_connectivity()
        node_coordinates_em = em_fem_mesh.get_all_node_coordinates()
        element_materials_em = em_fem_mesh.get_all_element_materials()
        
        # Phase 1: Filter elements to keep only those with elastic materials
        # Fortran: type_el
        element_materials_ac = []  # Material index for each AC element
        # Fortran: m_elnd_to_mshpt
        element_node_connectivity_ac_tmp = np.zeros(
            np.shape(element_node_connectivity_em), dtype=np.int64
        )
        element_convert_table = {}  # Maps AC element index -> EM element index
        
        num_elements_ac = 0
        for el in range(num_elements_em):
            mat_idx = element_materials_em[el]
            
            if elastic_props.is_elastic_material_index(mat_idx):
                element_materials_ac.append(
                    elastic_props.active_material_index(mat_idx)
                )
                element_convert_table[num_elements_ac] = el
                element_node_connectivity_ac_tmp[:, num_elements_ac] = element_node_connectivity_em[:, el]
                num_elements_ac += 1
        
        if num_elements_ac == 0:
            raise ValueError(
                "AC mesh trimming produced no elastic elements. "
                "Check that elastic_props.is_elastic_material_index() identifies at least one material."
            )
        
        # Phase 2: Extract unique nodes used by AC elements
        ac_nodes_with_duplicates = []  # Collect all nodes (with duplicates)
        for el in range(num_elements_ac):
            ac_nodes_with_duplicates.extend(element_node_connectivity_ac_tmp[:6, el])
        
        # Remove duplicates and sort for consistent ordering
        # Note: nodes are 1-based Fortran indices at this point
        ac_nodes_1based = sorted(set(ac_nodes_with_duplicates))
        # Convert to 0-based Python indices for the mapping
        ac_nodes_0based = [n - 1 for n in ac_nodes_1based]
        num_nodes_ac = len(ac_nodes_0based)
        
        # Phase 3: Create node remapping table (1-based EM node idx from connectivity -> 0-based AC node idx)
        # Fortran: n_msh_pts (for num_nodes_ac)
        em_node_to_ac_node_map = {}  # Maps 1-based EM node idx (from connectivity) -> 0-based AC node idx
        node_convert_table = {}  # Maps 0-based AC node idx -> 0-based EM node idx
        for i, em_node_idx_0based in enumerate(ac_nodes_0based):
            em_node_idx_1based = em_node_idx_0based + 1
            em_node_to_ac_node_map[em_node_idx_1based] = i
            node_convert_table[i] = em_node_idx_0based
        
        # Phase 4: Remap element connectivity to use AC node indices
        # Fortran: m_elnd_to_mshpt
        element_node_connectivity_ac = []
        for i in range(6):
            el_tbl = []
            for el in range(num_elements_ac):
                el_tbl.append(em_node_to_ac_node_map[element_node_connectivity_ac_tmp[i][el]])
            element_node_connectivity_ac.append(el_tbl)
        
        # Convert to numpy array and adjust to Fortran indexing (1-based)
        element_node_connectivity_ac = np.array(element_node_connectivity_ac, dtype=np.int64) + 1
        
        # Phase 5: Extract coordinates for AC nodes
        # Fortran: v_mshpt_xy
        node_coordinates_ac = np.zeros((2, num_nodes_ac))
        for ac_node_idx, em_node_idx_0based in enumerate(ac_nodes_0based):
            node_coordinates_ac[:, ac_node_idx] = node_coordinates_em[:, em_node_idx_0based]
        
        return ACMeshTrimData(
            element_node_connectivity=element_node_connectivity_ac,
            node_coordinates=node_coordinates_ac,
            element_materials=element_materials_ac,
            element_convert_table=element_convert_table,
            node_convert_table=node_convert_table,
        )
