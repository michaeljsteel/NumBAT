#include "numbat_decl.h"

module class_Mesh

   use numbatmod
   use alloc

   private

   ! A class for holding the various mesh and node properties together.
   ! Also for exchanging mesh info between fortran and python
   type, public  :: MeshEM

      integer(8) :: n_msh_pts
      integer(8) :: n_msh_elts
      integer(8) :: n_elt_mats

      ! What kind of mesh point, size: [n_msh_pts]
      ! If on a Gmsh Physical Line, it is the index of that line.
      ! Else is zero
      ! If physical lines are only defined on outer boundary,
      ! it serves as a check for being an outer boundary mesh point
      ! ranges: [0..max_phys_index]
      integer(8), dimension(:), allocatable :: v_mshpt_physindex

      ! absolute coordinates of each node in the mesh  [2, n_msh_pts]
      double precision, dimension(:,:), allocatable :: v_mshpt_xy

      ! material index of each elt, size: [n_msh_elts]
      ! ranges:  [1..n_mat_el]
      integer(8), dimension(:), allocatable :: v_elt_material

      ! mapping of elt,node to index of mesh point, size: [6, n_msh_elts]
      ! .mail elt data, but transposed: tn[i,j] = ith node of elt j
      !  tn[1..3,j] = vertex nodes, [4..6, j] edge nodes
      ! ranges: [1..n_msh_pts]
      integer(8), dimension(:,:), allocatable :: m_elnd_to_mshpt

   contains

      procedure :: allocate => MeshEM_allocate
      !final :: MeshEM_destructor

      ! load mesh data from .mail file
      procedure :: load_mesh_tables => MeshEM_load_mesh_tables

      ! is mesh point lying on outer boundary (a GMsh PhysicalLine)
      procedure :: is_boundary_mesh_point => MeshEM_is_boundary_mesh_point
      procedure :: is_boundary_mesh_point_by_elt_node => MeshEM_is_boundary_mesh_point_by_elt_node

      ! find node indices and xy locations at a given element
      procedure :: find_nodes_for_elt => MeshEM_find_nodes_for_elt

      procedure :: node_phys_index_by_ref => MeshEM_node_phys_index_by_ref

      ! for returning to python front end
      procedure :: fill_python_arrays => MeshEM_fill_python_arrays




   end type MeshEM

! ---------------------------------------


! v_tags[1..14, 1..n_msh_elts] assigns a unique integer label to each _distinct_ entity
! some elements have the same value because they are the same physical mesh point
! accessed from a different element/node combo

! v_tags[1, :]     - tag of each face of the elt  (Not actually a DoF?)
! v_tags[2..4, :]  - tag of P2 edges 1..3 of the elt (common edges to two elements get the same tag)
! v_tags[5..7, :]  - tag of P3 vertices 1..3 of the elt (common vertices to two elements get the same tag)
! v_tags[8..13, :] - tag of P3 nodes 4..9 of the elt which surround each P2 edge (common vertices to two elements get the same tag)
! v_tags[14, :]    - tag of P3 interior node 10 of the elt


   type, public  :: MeshEntities

      integer(8) n_faces       ! num faces
      integer(8) n_edges       ! num physically distinct edges
      integer(8) n_msh_pts_p3  ! num physically distinct P3 points
      integer(8) n_entities    ! total num of physically distinct things (sum of last 3)

      ! v_tags[i, j]
      !  [N_ENTITY_PER_EL x n_msh_elts] = [14 x n_msh_elts]
      !  tn[1, j] = j  ! index of element
      !    -
      !  Each element has 1 face, 3 edges and 10 P3 nodes
      !  so v_tags has dimensions 14 x n_msh_elts
      integer(8), dimension(:,:), allocatable :: v_tags


      ! Positions indexed by tag
      ! [2 x n_entities ]
      ! [:, 1] - barycentre of the face
      ! [:, 2..4] - P2 edge positions cloned from MeshEM%v_mshpt_xy
      ! [:, 5..7] - P3 vertex positions cloned from MeshEM%v_mshpt_xy
      ! [:, 8..13] - P3 edge positions built from MeshEM%v_mshpt_xy at 1/3, 2/3 intervals
      ! [:, 14] - P3 interior at barycentre
      double precision, dimension(:,:), allocatable :: v_xy


      ! row 1: what physical type is corresponding node of entity. (essentially is boundary or not)
      ! row 2: what dimensionality is entity
      !   [2 x N_ENTITY_PER_EL]
      !   [1, :]  for face: interior (0)
      !   [1, :]  for edge: take from MeshEM%v_mshpt_physindex
      !   [1, :]  for P3 vertex: take from MeshEM%v_mshpt_physindex
      !   [1, :]  for P3 edges: take from P2 edge in MeshEM%v_mshpt_physindex
      !   [1, :]  for P3 interior point:  0

      !   [2, :]  dimension: face(2), P2 edge (1), P3 edges and interior (0)
      integer(8), dimension(:,:), allocatable :: v_ety_props

      ! Maps P2 edge nodes to surrounding vertices
      integer(8) edge_end_nodes(2,3)

      ! tmp workspace, emptied once constructed
      integer(8), dimension(:), allocatable :: visited

   contains

      procedure :: allocate => MeshEntities_allocate
      procedure :: build_mesh_tables => MeshEntities_build_mesh_tables
      procedure :: count_faces => MeshEntities_count_faces
      procedure :: count_edges => MeshEntities_count_edges
      procedure :: count_nodes_P3 => MeshEntities_count_nodes_P3
      procedure :: analyse_face_and_edges => MeshEntities_analyse_face_and_edges
      procedure :: analyse_p3_nodes => MeshEntities_analyse_p3_nodes


      procedure :: check_bdy_elements_are_consistent => MeshEntities_check_bdy_elements_are_consistent

   end type MeshEntities





   ! A class for holding the various mesh and node properties together.
   ! Also for exchanging mesh info between fortran and python
   type, public  :: MeshAC

      !character(FNAME_LENGTH) :: mesh_file

      integer(8) :: n_msh_pts
      integer(8) :: n_msh_elts
      integer(8) :: n_elt_mats

      ! What kind of mesh point [n_msh_pts]
      ! For EM,
      !   If on a Gmsh Physical Line, it is the index of that line.
      !   Else is zero
      !   If physical lines are only defined on outer boundary,
      !   it serves as a check for being an outer boundary mesh point
      ! For AC,
      !   Very commonly the outer boundary will be stripped away by vacuum and
      !    every element is zero
      integer(8), dimension(:), allocatable :: v_mshpt_physindex

      ! absolute coordinates of each node in the mesh  [2, n_msh_pts]
      double precision, dimension(:,:), allocatable :: v_mshpt_xy

      ! material of each elt, size: [n_msh_elts]
      ! ranges:  [1..n_mat_ac_el], these are indexed back into the full material list via
      ! elastic material translation table
      integer(8), dimension(:), allocatable :: v_elt_material

      ! mapping of elt,node to index of mesh point, [6 x n_msh_elts]
      ! .mail elt data, but transposed: tn[i,j] = ith node of elt j
      !  tn[1..3,j] = vertex nodes, [4..6, j] edge nodes

      integer(8), dimension(:,:), allocatable :: m_elnd_to_mshpt

   contains

      procedure :: allocate => MeshAC_allocate
      !final :: MeshEM_destructor

      procedure :: load_mesh_tables_from_scratch => MeshAC_load_mesh_tables_from_scratch
      procedure :: load_mesh_tables_from_py => MeshAC_load_mesh_tables_from_py

      procedure :: find_nodes_for_elt => MeshAC_find_nodes_for_elt
      procedure :: is_boundary_mesh_point => MeshAC_is_boundary_mesh_point

      !procedure :: fill_python_arrays => MeshEM_fill_python_arrays
      !procedure :: is_boundary_mesh_point_by_elt_node => MeshEM_is_boundary_mesh_point_by_elt_node
!
!      procedure :: node_phys_index_by_ref => MeshEM_node_phys_index_by_ref



   end type MeshAC


contains

#include "meshprops_impl_em.f90"


#include "meshentities_impl.f90"

#include "meshprops_impl_ac.f90"


end module class_Mesh
