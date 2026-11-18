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
        n_msh_el_em = em_fem_mesh.num_elements()
        m_elnd_to_mshpt_em = em_fem_mesh.get_all_element_connectivity()
        v_mshpt_xy_em = em_fem_mesh.get_all_node_coordinates()
        type_el_em = em_fem_mesh.get_all_element_materials()
        
        # Phase 1: Filter elements to keep only those with elastic materials
        type_el_ac = []  # Material index for each AC element
        m_elnd_to_mshpt_ac_tmp = np.zeros(
            np.shape(m_elnd_to_mshpt_em), dtype=np.int64
        )
        el_convert_tbl = {}  # Maps AC element index -> EM element index
        
        n_msh_el_ac = 0
        for el in range(n_msh_el_em):
            mat_idx = type_el_em[el]
            
            if elastic_props.is_elastic_material_index(mat_idx):
                type_el_ac.append(
                    elastic_props.active_material_index(mat_idx)
                )
                el_convert_tbl[n_msh_el_ac] = el
                m_elnd_to_mshpt_ac_tmp[:, n_msh_el_ac] = m_elnd_to_mshpt_em[:, el]
                n_msh_el_ac += 1
        
        if n_msh_el_ac == 0:
            raise ValueError(
                "AC mesh trimming produced no elastic elements. "
                "Check that elastic_props.is_elastic_material_index() identifies at least one material."
            )
        
        # Phase 2: Extract unique nodes used by AC elements
        nodes_ac_mul = []  # Collect all nodes (with duplicates)
        for el in range(n_msh_el_ac):
            nodes_ac_mul.extend(m_elnd_to_mshpt_ac_tmp[:6, el])
        
        # Remove duplicates and sort for consistent ordering
        # Note: nodes are 1-based Fortran indices at this point
        nodes_ac_1based = sorted(set(nodes_ac_mul))
        # Convert to 0-based Python indices for the mapping
        nodes_ac = [n - 1 for n in nodes_ac_1based]
        n_msh_pts_ac = len(nodes_ac)
        
        # Phase 3: Create node remapping table (0-based EM node index -> 0-based AC node index)
        d_nodes_2_acnodes = {}  # Maps 1-based EM node idx (from connectivity) -> 0-based AC node idx
        node_convert_tbl = {}  # Maps 0-based AC node idx -> 0-based EM node idx
        for i, em_node_idx_0based in enumerate(nodes_ac):
            em_node_idx_1based = em_node_idx_0based + 1
            d_nodes_2_acnodes[em_node_idx_1based] = i
            node_convert_tbl[i] = em_node_idx_0based
        
        # Phase 4: Remap element connectivity to use AC node indices
        m_elnd_to_mshpt_ac = []
        for i in range(6):
            el_tbl = []
            for el in range(n_msh_el_ac):
                el_tbl.append(d_nodes_2_acnodes[m_elnd_to_mshpt_ac_tmp[i][el]])
            m_elnd_to_mshpt_ac.append(el_tbl)
        
        # Convert to numpy array and adjust to Fortran indexing (1-based)
        m_elnd_to_mshpt_ac = np.array(m_elnd_to_mshpt_ac, dtype=np.int64) + 1
        
        # Phase 5: Extract coordinates for AC nodes
        v_mshpt_xy_ac = np.zeros((2, n_msh_pts_ac))
        for ac_node_idx, em_node_idx_0based in enumerate(nodes_ac):
            v_mshpt_xy_ac[:, ac_node_idx] = v_mshpt_xy_em[:, em_node_idx_0based]
        
        return ACMeshTrimData(
            element_node_connectivity=m_elnd_to_mshpt_ac,
            node_coordinates=v_mshpt_xy_ac,
            element_materials=type_el_ac,
            element_convert_table=el_convert_tbl,
            node_convert_table=node_convert_tbl,
        )
