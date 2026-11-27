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

from numbat.reporting import OutputLevel, nbprint

import copy
import itertools
import warnings
from collections.abc import Sequence
from math import sqrt
from numbers import Real
from pathlib import Path
from typing import Any, Optional, Union

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

from .. import reporting, runtime
from numbat.mesh import mail
from numbat.mesh.ac_trimming import ACMeshTrimmer
from numbat.mesh.mesh_data import MeshArrays
from numbat.nbtypes import SI_to_gmpercc, SI_um
from numbat.utils.numeric import is_float_pair, is_real_number, np_min_max


def validate_1D_cut_format(
    s_cut: str, val1: Union[Real, Sequence[Real]], val2: Optional[Sequence[Real]]
) -> None:
    """Validate the direction and coordinates requested for a one-dimensional cut.

    Args:
        s_cut: Cut direction, either ``"x"``, ``"y"``, or ``"line"``.
        val1: Fixed coordinate for an axis-aligned cut, or the start point of a
            line cut.
        val2: End point of a line cut. It is required only when ``s_cut`` is
            ``"line"``.
    """

    if s_cut not in ('x', 'y', 'line'):
        reporting.report_and_exit(f'Invalid plot cut type. Must be "x", "y" or "line". Received {s_cut}.')

    if s_cut.lower() == "x":
        if not is_real_number(val1):
            reporting.report_and_exit('Invalid value for x cut plot: val1 must be a float within the x-domain.')

    elif s_cut.lower() == "y":
        if not is_real_number(val1):
            reporting.report_and_exit('Invalid value for y cut plot: val1 must be a float within the y-domain.')

    else:
        if not is_float_pair(val1) or not is_float_pair(val2):
            reporting.report_and_exit('Invalid value for line cut plot: val1 and val2 must be tuples of two floats.')


# def omake_interper_f_2D(tri_triang6p, finder, vx_out, vy_out, nx, ny):
#    return lambda femsol: matplotlib.tri.LinearTriInterpolator(
#    tri_triang6p, femsol, trifinder=finder)(vx_out, vy_out).reshape(
#        nx, ny)


# TODO: an attempt to allow matplotlib triangulation through multiprocessing
# ProcessPoolExecutor, because neither Triangulations or nested functions are
# pickleable
class TriInterpWrapper:
    """Callable wrapper around a linear FEM triangle interpolator."""

    def __init__(
        self,
        tri_triang6p: matplotlib.tri.Triangulation,
        finder: matplotlib.tri.TriFinder,
        vx_out: np.ndarray,
        vy_out: np.ndarray,
        nx: int = 0,
        ny: int = 0,
    ) -> None:
        """Store the triangulation and output grid used for interpolation.

        Args:
            tri_triang6p: Triangulation in six-node element ordering.
            finder: Triangle finder associated with the physical mesh.
            vx_out: Flattened output x coordinates.
            vy_out: Flattened output y coordinates.
            nx: Output grid width, or zero for a one-dimensional sample.
            ny: Output grid height, or zero for a one-dimensional sample.
        """
        self.tri_triang6p = tri_triang6p
        self.finder = finder
        self.vx_out = vx_out
        self.vy_out = vy_out
        self.nx = nx  # Nonzero if 2D context
        self.ny = ny  # Nonzero if 2D context

    def __call__(self, femsol: np.ndarray) -> np.ma.MaskedArray:
        """Interpolate a FEM solution onto the stored sampling coordinates.

        Args:
            femsol: Scalar field values in six-node element ordering.

        Returns:
            Interpolated values, reshaped to ``(nx, ny)`` for a two-dimensional
            sampling grid.
        """
        #print('\n\nCalling TriInterpWrapper:', type(self.tri_triang6p))
        #print(dir(self.tri_triang6p))
        #print(self.tri_triang6p._cpp_triangulation)

        #del self.tri_triang6p._cpp_triangulation
        #self.tri_triang6p._cpp_triangulation = None

        #mptri = matplotlib.tri.LinearTriInterpolator(
        #                                            self.tri_triang6p, femsol, trifinder=self.finder)(
        #                                                    self.vx_out, self.vy_out).reshape(self.nx, self.ny)

        mptri = matplotlib.tri.LinearTriInterpolator(
                                                    self.tri_triang6p, femsol, trifinder=self.finder)(
                                                            self.vx_out, self.vy_out)

        if self.nx:
            mptri = mptri.reshape(self.nx, self.ny)


        return mptri

def make_interper_f_2D(
    tri_triang6p: matplotlib.tri.Triangulation,
    finder: matplotlib.tri.TriFinder,
    vx_out: np.ndarray,
    vy_out: np.ndarray,
    nx: int,
    ny: int,
) -> TriInterpWrapper:
    """Create an interpolator for a flattened two-dimensional sampling grid.

    Args:
        tri_triang6p: Triangulation in six-node element ordering.
        finder: Triangle finder associated with the physical mesh.
        vx_out: Flattened x coordinates of the sampling grid.
        vy_out: Flattened y coordinates of the sampling grid.
        nx: Number of samples along the x axis.
        ny: Number of samples along the y axis.

    Returns:
        Callable interpolator that reshapes evaluated fields to ``(nx, ny)``.
    """

#    def mif2D(femsol):
#        return matplotlib.tri.LinearTriInterpolator(
#            tri_triang6p, femsol, trifinder=finder
#        )(vx_out, vy_out).reshape(nx, ny)
#
#    return mif2D

    return TriInterpWrapper(tri_triang6p, finder, vx_out, vy_out, nx, ny)



def make_interper_f_1D(
    tri_triang6p: matplotlib.tri.Triangulation,
    finder: matplotlib.tri.TriFinder,
    vx_out: np.ndarray,
    vy_out: np.ndarray,
) -> TriInterpWrapper:
    """Create an interpolator for samples along a one-dimensional line.

    Args:
        tri_triang6p: Triangulation in six-node element ordering.
        finder: Triangle finder associated with the physical mesh.
        vx_out: X coordinates of the sampling line.
        vy_out: Y coordinates of the sampling line.

    Returns:
        Callable interpolator that evaluates fields at the supplied points.
    """

    # def mif1D(femsol):
    #     return matplotlib.tri.LinearTriInterpolator(
    #         tri_triang6p, femsol, trifinder=finder
    #     )(vx_out, vy_out)

    # return mif1D
    return TriInterpWrapper(tri_triang6p, finder, vx_out, vy_out)




class FemMesh:
    """Mirrors the FEM mesh constructed in Fortran.

    Used by FEMScalarFieldPlotter to understand FEM grid to sample from
    Simulation (EM and AC) to communicate with Fortran side.

    The properties can be constructed eitehr by reading an existing .mail file or by
    taking the properties from the Fortran side after a simulation has been run.
    """


    def __init__(self) -> None:
        """Initialise an empty mesh ready for EM or acoustic mesh data."""
        ## filename of original Gmsh-derived .mail file
        self.mesh_mail_fname = ""

        # Number of mesh nodes (points) in the FEM calculation
        # (Different to original underlying mesh because of elt sub-triangulation)
        # Fortran: n_msh_pts
        self.num_nodes = 0

        # Number of elements in .msh mesh file
        # Fortran: n_msh_el
        self.num_elements = 0

        #
        # Quantities made by python
        #

        # Array[0:num_elements] of index of elt's material into list of active em or materials.
        #    Values in range:  em simulation - [1, optical_props.n_mats_em]
        #                      ac simulation - [1, elastic_props.n_mats_ac]
        # Fortran: type_el or v_el_2_mat_idx
        self.element_materials = None

        # An elastic el_conv_tbl_n (a dictionary!), completely diff to self.type_el
        # AC-specific: should be moved to ElasticProperties
        # Fortran: (no direct mapping)
        self._typ_el_AC = None

        #
        # Quantities made by fortran
        #
        # Map of each element to its 6 nodes by node index (1..6). shape = (6, num_elements])
        # Fortran: m_elnd_to_mshpt
        self.element_node_connectivity = None

        # Line or surface index of a node [1..num_gmsh_types],     shape = (num_nodes,1)
        # Fortran: type_nod
        self.node_physical_indices = None

        # physical scaled x-y, locations of every node             shape = (2, num_nodes)
        # Fortran: v_mshpt_xy
        self.node_coordinates = None

        # Nodes per each element (is always 6)
        self.nodes_per_element = 6

        self._ac_mesh_from_em = True  # Always True

        # Dict map of elastically active mesh elements into original em mesh elements.
        #    Length: num_elements(ac).   Values in [1..num_elements(em)]
        # AC-specific: should be moved to MaterialMeshMapping
        # Fortran: (no direct mapping)
        self.element_convert_table = None

        # measured from nbgmsh to have them early before we know node_coordinates
        self._extents = None

    # =========================================================================
    # DEPRECATED PROPERTIES - For backward compatibility only
    # Use the new direct attributes or accessor methods instead
    # =========================================================================

    @property
    def n_msh_pts(self) -> int:
        """Deprecated: Use num_nodes directly instead."""
        warnings.warn(
            "Direct access to n_msh_pts is deprecated. Use num_nodes directly instead.",
            DeprecationWarning,
            stacklevel=2
        )
        return self.num_nodes

    @n_msh_pts.setter
    def n_msh_pts(self, value: int) -> None:
        """Deprecated: Use num_nodes directly instead."""
        warnings.warn(
            "Direct assignment to n_msh_pts is deprecated.",
            DeprecationWarning,
            stacklevel=2
        )
        self.num_nodes = value

    @property
    def n_msh_el(self) -> int:
        """Deprecated: Use num_elements directly instead."""
        warnings.warn(
            "Direct access to n_msh_el is deprecated. Use num_elements directly instead.",
            DeprecationWarning,
            stacklevel=2
        )
        return self.num_elements

    @n_msh_el.setter
    def n_msh_el(self, value: int) -> None:
        """Deprecated: Use num_elements directly instead."""
        warnings.warn(
            "Direct assignment to n_msh_el is deprecated.",
            DeprecationWarning,
            stacklevel=2
        )
        self.num_elements = value

    @property
    def v_el_2_mat_idx(self) -> Optional[np.ndarray]:
        """Deprecated: Use element_materials directly instead."""
        warnings.warn(
            "Direct access to v_el_2_mat_idx is deprecated. "
            "Use element_materials directly instead.",
            DeprecationWarning,
            stacklevel=2
        )
        return self.element_materials

    @v_el_2_mat_idx.setter
    def v_el_2_mat_idx(self, value: Optional[np.ndarray]) -> None:
        """Deprecated: Use element_materials directly instead."""
        warnings.warn(
            "Direct assignment to v_el_2_mat_idx is deprecated.",
            DeprecationWarning,
            stacklevel=2
        )
        self.element_materials = value

    @property
    def m_elnd_to_mshpt(self) -> Optional[np.ndarray]:
        """Deprecated: Use element_node_connectivity directly instead."""
        warnings.warn(
            "Direct access to m_elnd_to_mshpt is deprecated. "
            "Use element_node_connectivity directly instead.",
            DeprecationWarning,
            stacklevel=2
        )
        return self.element_node_connectivity

    @m_elnd_to_mshpt.setter
    def m_elnd_to_mshpt(self, value: Optional[np.ndarray]) -> None:
        """Deprecated: Use element_node_connectivity directly instead."""
        warnings.warn(
            "Direct assignment to m_elnd_to_mshpt is deprecated.",
            DeprecationWarning,
            stacklevel=2
        )
        self.element_node_connectivity = value

    @property
    def v_mshpt_xy(self) -> Optional[np.ndarray]:
        """Deprecated: Use node_coordinates directly instead."""
        warnings.warn(
            "Direct access to v_mshpt_xy is deprecated. "
            "Use node_coordinates directly instead.",
            DeprecationWarning,
            stacklevel=2
        )
        return self.node_coordinates

    @v_mshpt_xy.setter
    def v_mshpt_xy(self, value: Optional[np.ndarray]) -> None:
        """Deprecated: Use node_coordinates directly instead."""
        warnings.warn(
            "Direct assignment to v_mshpt_xy is deprecated.",
            DeprecationWarning,
            stacklevel=2
        )
        self.node_coordinates = value

    @property
    def node_physindex(self) -> Optional[np.ndarray]:
        """Deprecated: Use node_physical_indices directly instead."""
        warnings.warn(
            "Direct access to node_physindex is deprecated. "
            "Use node_physical_indices directly instead.",
            DeprecationWarning,
            stacklevel=2
        )
        return self.node_physical_indices

    @node_physindex.setter
    def node_physindex(self, value: Optional[np.ndarray]) -> None:
        """Deprecated: Use node_physical_indices directly instead."""
        warnings.warn(
            "Direct assignment to node_physindex is deprecated.",
            DeprecationWarning,
            stacklevel=2
        )
        self.node_physical_indices = value

    @property
    def el_convert_tbl(self) -> Optional[dict]:
        """Deprecated: Use element_convert_table directly instead."""
        warnings.warn(
            "Direct access to el_convert_tbl is deprecated. "
            "Use element_convert_table directly instead.",
            DeprecationWarning,
            stacklevel=2
        )
        return self.element_convert_table

    @el_convert_tbl.setter
    def el_convert_tbl(self, value: Optional[dict]) -> None:
        """Deprecated: Use element_convert_table directly instead."""
        warnings.warn(
            "Direct assignment to el_convert_tbl is deprecated.",
            DeprecationWarning,
            stacklevel=2
        )
        self.element_convert_table = value

    @property
    def typ_el_AC(self) -> Optional[Any]:
        """Deprecated: This should be in ElasticProperties, not FemMesh."""
        warnings.warn(
            "Direct access to typ_el_AC is deprecated. "
            "This property will move to ElasticProperties in a future version.",
            DeprecationWarning,
            stacklevel=2
        )
        return self._typ_el_AC

    @typ_el_AC.setter
    def typ_el_AC(self, value: Optional[Any]) -> None:
        """Deprecated: This should be in ElasticProperties, not FemMesh."""
        warnings.warn(
            "Direct assignment to typ_el_AC is deprecated.",
            DeprecationWarning,
            stacklevel=2
        )
        self._typ_el_AC = value

    @property
    def extents(self) -> Optional[list]:
        """Deprecated: Use get_xy_limits() instead."""
        warnings.warn(
            "Direct access to extents is deprecated. Use get_xy_limits() instead.",
            DeprecationWarning,
            stacklevel=2
        )
        return self._extents

    @extents.setter
    def extents(self, value: Optional[list]) -> None:
        """Deprecated: Use get_xy_limits() instead."""
        warnings.warn(
            "Direct assignment to extents is deprecated.",
            DeprecationWarning,
            stacklevel=2
        )
        self._extents = value

    @property
    def ac_mesh_from_em(self) -> bool:
        """Deprecated: This property is always True and will be removed."""
        warnings.warn(
            "Access to ac_mesh_from_em is deprecated. This property will be removed.",
            DeprecationWarning,
            stacklevel=2
        )
        return self._ac_mesh_from_em

    @ac_mesh_from_em.setter
    def ac_mesh_from_em(self, value: bool) -> None:
        """Deprecated: This property is always True and will be removed."""
        warnings.warn(
            "Assignment to ac_mesh_from_em is deprecated.",
            DeprecationWarning,
            stacklevel=2
        )
        self._ac_mesh_from_em = value

    # =========================================================================
    # PUBLIC ACCESSORS - Recommended API
    # Note: num_nodes and num_elements are direct integer attributes.
    # For detailed access patterns, use the methods below.
    # =========================================================================

    def get_node_coordinate(self, node_idx: int) -> tuple[float, float]:
        """Get (x, y) coordinate of a single node.
        
        Args:
            node_idx: 0-based mesh node index.
            
        Returns:
            Tuple of (x, y) coordinates in physical units.
            
        Raises:
            ValueError: If node_idx is out of bounds.
        """
        if self.node_coordinates is None:
            raise RuntimeError("Mesh coordinates not loaded")
        if node_idx < 0 or node_idx >= self.num_nodes:
            raise ValueError(f"Node index {node_idx} out of bounds [0, {self.num_nodes})")
        return (float(self.node_coordinates[0, node_idx]), float(self.node_coordinates[1, node_idx]))

    def get_all_node_coordinates(self) -> np.ndarray:
        """Return array of all node coordinates.
        
        Returns:
            Array of shape (2, n_nodes) where [:, i] is (x, y) of node i.
            Returns a copy to prevent modification of internal data.
        """
        if self.node_coordinates is None:
            raise RuntimeError("Mesh coordinates not loaded")
        return self.node_coordinates.copy()

    def element_node_indices(self, elem_idx: int) -> np.ndarray:
        """Get 6 node indices for an element.
        
        Args:
            elem_idx: 0-based element index.
            
        Returns:
            Array of 6 node indices (1-indexed for Fortran compatibility).
            
        Raises:
            ValueError: If elem_idx is out of bounds.
        """
        if self.element_node_connectivity is None:
            raise RuntimeError("Element connectivity not loaded")
        if elem_idx < 0 or elem_idx >= self.num_elements:
            raise ValueError(f"Element index {elem_idx} out of bounds [0, {self.num_elements})")
        return self.element_node_connectivity[:, elem_idx]

    def get_all_element_connectivity(self) -> np.ndarray:
        """Return element-to-node mapping for all elements.
        
        Returns:
            Array of shape (6, n_elements) containing 1-based node indices.
            Returns a copy to prevent modification of internal data.
        """
        if self.element_node_connectivity is None:
            raise RuntimeError("Element connectivity not loaded")
        return self.element_node_connectivity.copy()

    def element_material_index_1based(self, elem_idx: int) -> int:
        """Get 1-based material index for an element (Fortran convention).
        
        Args:
            elem_idx: 0-based element index.
            
        Returns:
            1-based material index.
            
        Raises:
            ValueError: If elem_idx is out of bounds.
        """
        if self.element_materials is None:
            raise RuntimeError("Element material mapping not loaded")
        if elem_idx < 0 or elem_idx >= self.num_elements:
            raise ValueError(f"Element index {elem_idx} out of bounds [0, {self.num_elements})")
        return int(self.element_materials[elem_idx])

    def get_all_element_materials(self) -> np.ndarray:
        """Return material index for all elements (1-based, Fortran convention).
        
        Returns:
            Array of shape (n_elements,) with 1-based material indices.
            Returns a copy to prevent modification of internal data.
        """
        if self.element_materials is None:
            raise RuntimeError("Element material mapping not loaded")
        return self.element_materials.copy()

    def node_physical_index(self, node_idx: int) -> int:
        """Get physical boundary index for a node (EM meshes only).
        
        Args:
            node_idx: 0-based mesh node index.
            
        Returns:
            Physical boundary index (1-based, or 0 if not on boundary).
            
        Raises:
            ValueError: If node_idx is out of bounds.
        """
        if self.node_physical_indices is None:
            raise RuntimeError("Node physical indices not loaded (AC mesh?)")
        if node_idx < 0 or node_idx >= self.num_nodes:
            raise ValueError(f"Node index {node_idx} out of bounds [0, {self.num_nodes})")
        return int(self.node_physical_indices[node_idx])

    def get_all_node_physical_indices(self) -> Optional[np.ndarray]:
        """Return physical boundary indices for all nodes (EM meshes only).
        
        Returns:
            Array of shape (n_nodes,) with physical indices, or None if not set.
            Returns a copy to prevent modification of internal data.
        """
        if self.node_physical_indices is None:
            return None
        return self.node_physical_indices.copy()

    def get_element_convert_table(self) -> Optional[dict]:
        """Get mapping from AC element indices to original EM element indices.
        
        Only populated for AC meshes trimmed from EM meshes.
        
        Returns:
            Dict mapping AC element index (int) -> EM element index (int), or None if not an AC mesh.
            Used to translate between trimmed AC mesh and full EM mesh.
        """
        return self.element_convert_table

    def build_from_gmsh_mail(self, wguide: Any) -> None:
        """Load the electromagnetic mesh and physical coordinates from a waveguide.

        Args:
            wguide: Waveguide structure providing its ``.mail`` mesh filename,
                domain dimensions, and optical material properties.
        """
        # Takes list of all material refractive indices
        # Discards any that are zero

        # (MJS: Not sure about the counting from 1, possibly needed for fortran?)

        self.mesh_mail_fname = wguide.mesh_mail_fname
        mesh = mail.MailData(wguide.mesh_mail_fname)

        # mesh = wguide.get_mail_mesh_data()

        self.element_materials = mesh.v_elts_mat

        # Electromagnetic so the mesh properties come straight from the Mail file.
        self.num_nodes = mesh.n_msh_pts
        self.num_elements = mesh.n_msh_elts

        self.element_node_connectivity = mesh.v_elts[:, :6].T
        self.node_coordinates = np.zeros([2, self.num_nodes])

        # Messy: Mail file does not include the domain scaling which is in nm, not microns
        # Everything including yvalues is scaled by domain_x, not domain_y
        self.node_coordinates[0, :] = mesh.v_x * wguide.domain_x * 0.001
        self.node_coordinates[1, :] = (mesh.v_y * wguide.domain_x * 0.001)

        self._extents = [
            np.min(self.node_coordinates[0, :]),
            np.max(self.node_coordinates[0, :]),
            np.min(self.node_coordinates[1, :]),
            np.max(self.node_coordinates[1, :]),
        ]

        self._report_em_mesh_properties(wguide)

    def _report_em_mesh_properties(self, wguide):
        opt_props = wguide.optical_props
        nbprint(
            f"\n The final EM sim mesh has {self.num_nodes} nodes, {self.num_elements} elements and {opt_props.n_mats_em} element types (materials)."
        )

        #matvals = list(wguide.d_materials.values())[: opt_props.n_mats_em]
        matvals = wguide.get_optical_materials()

        #print(" The material index table is:", opt_props.el_conv_table_n, "\n")
        nbprint(f" There are {opt_props.n_mats_em} active materials:")
        for im, m in enumerate(matvals):
            # +1 because materials are reported by their Fortran index
            nbprint(f'  {m.material_name+",":20} n = {opt_props.v_refindexn[im]:.5f}, mat. index = {im+1}.')


    def store_fortran_em_mesh_properties(
        self,
        type_el: np.ndarray,
        node_physindex: np.ndarray,
        m_elnd_to_mshpt: np.ndarray,
        v_mshpt_xy: np.ndarray,
    ) -> None:
        """Store electromagnetic mesh arrays returned by the Fortran solver.

        Args:
            type_el: Per-element Fortran material indices (Fortran: type_el).
            node_physindex: Per-node physical-boundary indices (Fortran: type_nod).
            m_elnd_to_mshpt: Six-node element-to-mesh-node mapping (Fortran: m_elnd_to_mshpt).
            v_mshpt_xy: Physical x and y coordinates of mesh nodes (Fortran: v_mshpt_xy).
        """

        self.element_materials = type_el
        self.element_node_connectivity = m_elnd_to_mshpt
        self.node_coordinates = v_mshpt_xy
        self.node_physical_indices = node_physindex
        self.num_nodes = v_mshpt_xy.shape[1]
        self.num_elements = type_el.shape[0]

        # Set extents from loaded coordinates
        if v_mshpt_xy.shape[1] > 0:
            self._extents = [
                np.min(v_mshpt_xy[0, :]),
                np.max(v_mshpt_xy[0, :]),
                np.min(v_mshpt_xy[1, :]),
                np.max(v_mshpt_xy[1, :]),
            ]

    def store_fortran_ac_mesh_properties(
        self, type_el: np.ndarray, m_elnd_to_mshpt: np.ndarray, v_mshpt_xy: np.ndarray
    ) -> None:
        """Store acoustic mesh arrays returned by the Fortran solver.

        Args:
            type_el: Per-element Fortran material indices (Fortran: type_el).
            m_elnd_to_mshpt: Six-node element-to-mesh-node mapping (Fortran: m_elnd_to_mshpt).
            v_mshpt_xy: Physical x and y coordinates of mesh nodes (Fortran: v_mshpt_xy).
        """

        self.element_materials = type_el
        self.element_node_connectivity = m_elnd_to_mshpt
        self.node_coordinates = v_mshpt_xy
        self.num_nodes = v_mshpt_xy.shape[1]
        self.num_elements = type_el.shape[0]

        # Set extents from loaded coordinates
        if v_mshpt_xy.shape[1] > 0:
            self._extents = [
                np.min(v_mshpt_xy[0, :]),
                np.max(v_mshpt_xy[0, :]),
                np.min(v_mshpt_xy[1, :]),
                np.max(v_mshpt_xy[1, :]),
            ]

    def ac_build_from_em(self, wguide: Any, em_fem: "FemMesh") -> None:
        """Builds the FemMesh for an acoustic simulation based on an existing EM mesh
        and the known elastic properties.

        Called from ACSimulation.__init__()
        """

        el_props = wguide.elastic_props

        # Use ACMeshTrimmer to extract and remap AC mesh from EM mesh
        trim_data = ACMeshTrimmer.trim_em_mesh_to_ac(em_fem, el_props)

        # Store trimmed mesh data in FemMesh
        self.element_convert_table = trim_data.element_convert_table
        self.num_nodes = trim_data.num_nodes
        self.num_elements = trim_data.num_elements
        self.element_node_connectivity = trim_data.element_node_connectivity
        self.element_materials = trim_data.element_materials
        self.node_coordinates = trim_data.node_coordinates

        # AC FEM uses Neumann B.C.s so node_physical_indices is totally irrelevant!
        self.node_physical_indices = np.zeros(trim_data.num_nodes)

        # Set extents from loaded coordinates
        if trim_data.node_coordinates.shape[1] > 0:
            self._extents = [
                np.min(trim_data.node_coordinates[0, :]),
                np.max(trim_data.node_coordinates[0, :]),
                np.min(trim_data.node_coordinates[1, :]),
                np.max(trim_data.node_coordinates[1, :]),
            ]

        self._report_ac_mesh_properties(el_props)

    def ac_build_from_gmsh_mail(self, wguide: Any) -> None:
        """Builds the FemMesh for an acoustic simulation directly from the Mail mesh file.

        This creates the full EM mesh, then trims it for AC analysis without
        running the EM eigensolve. The resulting coordinates are in metres,
        as for an AC mesh built from an EM SimResult. Used when ACSimulation
        is called without an existing EM SimResult.

        Args:
            wguide: Waveguide structure providing its ``.mail`` mesh filename,
                domain dimensions, elastic material properties, and optical
                material properties (needed for mesh trimming).
        """
        # First, build the full EM mesh from the Mail file
        self.build_from_gmsh_mail(wguide)

        # Then trim it for AC using the elastic properties
        el_props = wguide.elastic_props
        trim_data = ACMeshTrimmer.trim_em_mesh_to_ac(self, el_props)

        # Replace EM mesh with trimmed AC mesh
        self.element_convert_table = trim_data.element_convert_table
        self.num_nodes = trim_data.num_nodes
        self.num_elements = trim_data.num_elements
        self.element_node_connectivity = trim_data.element_node_connectivity
        self.element_materials = trim_data.element_materials
        # The Mail-based Python mesh is in micrometres; the AC solver needs metres.
        self.node_coordinates = trim_data.node_coordinates * SI_um

        # AC FEM uses Neumann B.C.s so node_physical_indices is totally irrelevant!
        self.node_physical_indices = np.zeros(trim_data.num_nodes)

        # Set extents from loaded coordinates
        if self.node_coordinates.shape[1] > 0:
            self._extents = [
                np.min(self.node_coordinates[0, :]),
                np.max(self.node_coordinates[0, :]),
                np.min(self.node_coordinates[1, :]),
                np.max(self.node_coordinates[1, :]),
            ]

        self._report_ac_mesh_properties(el_props)

    def get_xy_limits(self) -> tuple[float, float, float, float]:
        """Return minimum and maximum physical coordinates of the mesh."""
        if self.node_coordinates is None:
            raise RuntimeError("Mesh not loaded")
        x_min, x_max = np_min_max(self.node_coordinates[0, :])
        y_min, y_max = np_min_max(self.node_coordinates[1, :])
        return x_min, x_max, y_min, y_max

    def _report_ac_mesh_properties(self, el_props):

        nbprint(
            f"\n The final elastic sim mesh has {self.num_nodes} nodes, {self.num_elements} mesh elements and {len(el_props.typ_el_AC)} element types (materials)."
        )

        nbprint(" The material index table is:", el_props.typ_el_AC, "\n")
        nbprint(f" There are {el_props.n_mats_ac} active materials:")
        for mat in el_props.v_acoustic_mats:
            vphases = mat.Vac_phase()
            nbprint(
                f'   {mat.material_name+",":20} rho = {mat.rho*SI_to_gmpercc:.3f} g/cc, v_i = [{vphases[0]:.4f}, {vphases[1]:.4f}, {vphases[2]:.4f}] km/s.'
            )
        # print("  elt2nodes index map", self.m_elnd_to_mshpt)


    def make_sub_triangulation(self) -> tuple[list[list[int]], list[list[int]]]:
        """Split every quadratic mesh element into four linear triangles.

        Returns:
            Triangle-node indices in six-node element ordering and global
            mesh-node ordering, respectively.
        """
        if self.element_node_connectivity is None:
            raise RuntimeError("Element connectivity not loaded")

        v_triang6p = []
        for idx in range(0, 6 * self.num_elements, 6):
            triangles = [
                [idx + 0, idx + 3, idx + 5],
                [idx + 1, idx + 4, idx + 3],
                [idx + 2, idx + 5, idx + 4],
                [idx + 3, idx + 4, idx + 5],
            ]
            v_triang6p.extend(triangles)

        # v_triang1p maps the same triangles to the actual nodes defined by Gmsh, of which there are only num_nodes
        tabnod_py = self.element_node_connectivity.T - 1  # shift fortran to python indexing

        v_triang1p = []
        for i_el in np.arange(self.num_elements):
            triangles = [
                [tabnod_py[i_el, 0], tabnod_py[i_el, 3], tabnod_py[i_el, 5]],
                [tabnod_py[i_el, 1], tabnod_py[i_el, 4], tabnod_py[i_el, 3]],
                [tabnod_py[i_el, 2], tabnod_py[i_el, 5], tabnod_py[i_el, 4]],
                [tabnod_py[i_el, 3], tabnod_py[i_el, 4], tabnod_py[i_el, 5]],
            ]
            v_triang1p.extend(triangles)

        return v_triang6p, v_triang1p

    def _save_triangulation_plots(self, triang1p, triang6p, v_mshpt_xy):
        fig, axs = plt.subplots(1, 1)
        axs[0].triplot(triang1p, linewidth=0.5)
        # axs[1].triplot(triang6p, linewidth=.5)
        # for ax in axs:
        axs[0].set_aspect(1.0)
        axs[0].scatter(v_mshpt_xy[0, :], v_mshpt_xy[1, :], s=2, c="red")

        axs[0].set_xlabel(r"$x$ [μm]")
        axs[0].set_ylabel(r"$y$ [μm]")

        pref = runtime.outprefix()
        fname = pref + f"-{'ac' if self.sim_result.is_AC else 'em'}_triplots.png"
        save_and_close_figure(fig, fname)

    def make_interpolator_for_grid(
        self, vx_out: np.ndarray, vy_out: np.ndarray, nx: int, ny: int
    ) -> TriInterpWrapper:
        """Constructs interpolator to map triangular fem grid scalar functions onto requested rectangular grid.

        Args:
            vx_out (array(float)): Flattened x-values of rectangular grid
            vy_out (array(float)): Flattened x-values of rectangular grid
            nx (int): x-dimension of rectangular grid
            ny (int): y-dimension of rectangular grid

        Returns:
            make_interpolator_for_grid: The interpolator ready for use
        """

        v_triang6p, v_triang1p = self.make_sub_triangulation()

        # This is for testing only. Normally turn off
        check_tris = False
        if check_tris:
            check_triangulation(self.node_coordinates[0, :], self.node_coordinates[1, :], v_triang1p)

        v_x6p, v_y6p = self.get_fullmesh_nodes_xy()

        # triangulations:  x and y coords of all points, list of triangles defined by triples of indices of the points
        tri_triang6p = matplotlib.tri.Triangulation(v_x6p, v_y6p, v_triang6p)
        tri_triang1p = matplotlib.tri.Triangulation(
            self.node_coordinates[0, :], self.node_coordinates[1, :], v_triang1p
        )

        draw_triangulation = False
        if draw_triangulation:
            pl_tri_triang6p = matplotlib.tri.Triangulation(
                v_x6p / SI_um, v_y6p / SI_um, v_triang6p
            )
            pl_tri_triang1p = matplotlib.tri.Triangulation(
                self.node_coordinates[0, :] / SI_um, self.node_coordinates[1, :] / SI_um, v_triang1p
            )

            self._save_triangulation_plots(
                pl_tri_triang1p, pl_tri_triang6p, self.node_coordinates / SI_um
            )

        # The trifinder only works with triang1p, not triang6p.
        # The latter is apparently an 'invalid triangulation'.
        finder = tri_triang1p.get_trifinder()

        if ny > 1:  # a 2D sampling, not a line cut
            interper_f = make_interper_f_2D(tri_triang6p, finder, vx_out, vy_out, nx, ny)
        else:
            interper_f = make_interper_f_1D(tri_triang6p, finder, vx_out, vy_out)

        return interper_f

    def rect_grid(self, n_points: int) -> tuple[np.ndarray, np.ndarray]:
        """Create an approximately square-pixel rectangular grid over the mesh.

        Args:
            n_points: Nominal linear resolution for the grid.

        Returns:
            X and y coordinate vectors spanning the mesh extents.
        """
        if self._extents is None:
            raise RuntimeError("Mesh extents not calculated")

        x_min, x_max = self._extents[0:2]
        y_min, y_max = self._extents[2:]

        area = abs((x_max - x_min) * (y_max - y_min))
        n_pts_x = int(n_points * abs(x_max - x_min) / np.sqrt(area))
        n_pts_y = int(n_points * abs(y_max - y_min) / np.sqrt(area))

        v_regx = np.linspace(x_min, x_max, n_pts_x)
        v_regy = np.linspace(y_min, y_max, n_pts_y)
        return v_regx, v_regy

    def element_to_material_index(self, i_el: int) -> int:
        """Return the zero-based material index for one mesh element.

        Args:
            i_el: Zero-based mesh-element index.

        Returns:
            Zero-based material index associated with the element.
        """
        # -1 because FEM material indices are fortran unit-indexed
        return int(self.element_material_index_1based(i_el) - 1)



    def get_fullmesh_nodes_xy(self) -> tuple[np.ndarray, np.ndarray]:
        """Returns vectors of x and y physical positions from the 6 nodes of each element
        in the standard order.

        Many nodes will appear multiple times because they are edges (twice) or vertex
        nodes (>=2).

        Used only by the raw mode plotter in plotmoderaw.py.
        """
        if self.element_node_connectivity is None or self.node_coordinates is None:
            raise RuntimeError("Mesh not loaded")

        v_x6p = np.zeros(6 * self.num_elements)
        v_y6p = np.zeros(6 * self.num_elements)

        tabnod_py = self.element_node_connectivity.T - 1  # shift fortran to python indexing

        i = 0
        for i_el in range(self.num_elements):
            for i_node in range(6):
                i_ex = tabnod_py[i_el, i_node]
                v_x6p[i], v_y6p[i] = self.node_coordinates[:, i_ex]
                i += 1

        return v_x6p, v_y6p




# Checks of mesh and triangles satisfy conditions for triangulation
# Quadratic algorithm. Use on the smallest grid possible
def check_triangulation(vx: np.ndarray, vy: np.ndarray, triangs: Sequence[Sequence[int]]) -> None:
    """Report duplicate points and duplicate triangles in a triangulation.

    Args:
        vx: X coordinates of triangulation nodes.
        vy: Y coordinates of triangulation nodes.
        triangs: Node-index triples defining the triangles to inspect.
    """
    # are points unique
    nbprint("\n\nChecking triangulation goodness")
    npts = len(vx)
    dsepmin = 1e6
    dsi = 0
    dsj = 0
    for i in range(npts):
        for j in range(i + 1, npts):
            dsep = sqrt((vx[i] - vx[j]) ** 2 + (vy[i] - vy[j]) ** 2)
            if dsep < dsepmin:
                dsepmin = dsep
                dsi = i
                dsj = j

    nbprint("  Closest space of triangle points was", dsepmin)
    if dsepmin < 1e-11:
        msg = f"Point collision at {dsi}, {dsj}: ({vx[dsi]},{vy[dsi]}) =  ({vx[dsj]},{vy[dsj]})."
        msg += "\nIt seems the mesh grid reordering has failed."
        reporting.report_and_exit(msg)

    # is list of triangles unique
    s_vtri = set()
    clean = True
    for tri in triangs:
        stri = str(tri)
        if stri in s_vtri:
            nbprint("        Double triangle at", stri, level=OutputLevel.ALWAYS)
            clean = False
        else:
            s_vtri.add(stri)
    if clean:
        nbprint("  No doubled triangles found")
    else:
        nbprint("  Found doubled triangles", level=OutputLevel.ALWAYS)