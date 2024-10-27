#include "numbat_decl.h"

module class_MeshRaw

   use numbatmod
   use alloc

   private

   ! A class for holding the various mesh and node properties together.
   ! Also for exchanging mesh info between fortran and python
   type, public  :: MeshRaw

      !character(FNAME_LENGTH) :: mesh_file

      integer(8) :: n_msh_pts
      integer(8) :: n_msh_el
      integer(8) :: n_elt_mats

      ! what kind of mesh point [n_msh_pts]
      ! basically the Gmsh physindex.  Nonzero = boundary node
      integer(8), dimension(:), allocatable :: v_nd_physindex

      ! absolute coordinates of each node in the mesh  [2, n_msh_pts]
      double precision, dimension(:,:), allocatable :: v_nd_xy

      ! material of each elt [n_msh_el]
      integer(8), dimension(:), allocatable :: el_material

      ! mapping of elt,node to index of mesh point, [6 x n_msh_el]
      ! .mail elt data, but transposed: tn[i,j] = ith node of elt j
      !  tn[1..3,j] = vertex nodes, [4..6, j] edge nodes

      integer(8), dimension(:,:), allocatable :: elnd_to_mshpt

   contains

      procedure :: allocate => MeshRaw_allocate
      !final :: MeshRaw_destructor

      procedure :: find_nodes_for_elt => MeshRaw_find_nodes_for_elt
      procedure :: fill_python_arrays => MeshRaw_fill_python_arrays
      procedure :: is_boundary_node => MeshRaw_is_boundary_node
      procedure :: is_boundary_node_2 => MeshRaw_is_boundary_node_2

      procedure :: node_phys_index_by_ref => MeshRaw_node_phys_index_by_ref

      procedure :: construct_node_tables => MeshRaw_construct_node_tables


   end type MeshRaw

! ---------------------------------------


! v_tags[1..14, 1..n_msh_el] assigns a unique integer label to each entity
! v_tags[1, :]     - label of each face of the elt
! v_tags[2..4, :]  - label of P2 edges 1..3 of the elt (common edges to two elements get the same label)
! v_tags[5..7, :] - label of P3 vertices 1..3 of the elt (common vertices to two elements get the same label)
! v_tags[8..13, :] - label of P3 vertices 4..9 of the elt which surround each P2 edge (common vertices to two elements get the same label)
! v_tags[14, :] - label of P3 interior node 10 of the elt


   type, public  :: MeshEntities

      integer(8) n_faces       ! num faces
      integer(8) n_edges       ! num physically distinct edges
      integer(8) n_msh_pts_p3  ! num physically distinct P3 points
      integer(8) n_entities         ! total num of physically distinct things (sum of last 3)

      ! v_tags[i, j] 14 x n_msh_el
      !  tn[1, j] = j  ! index of element
      !    -
      !  Each element has 1 face, 3 edges and 10 P3 nodes
      !  so v_tags has dimensions 14 x n_msh_el
      integer(8), dimension(:,:), allocatable :: v_tags


      ! Positions indexed by tag
      ! [2 x n_ddl]
      ! [:, 1] - barycentre of the face
      ! [:, 2..4] - P2 edge positions cloned from MeshRaw%v_nd_xy
      ! [:, 5..7] - P3 vertex positions cloned from MeshRaw%v_nd_xy
      ! [:, 8..13] - P3 edge positions built from MeshRaw%v_nd_xy at 1/3, 2/3 intervals
      ! [:, 14] - P3 interior at barycentre
      double precision, dimension(:,:), allocatable :: v_xy


      !   [2 x n_ddl]
      !   [1, :]  for face: interior (0)
      !   [1, :]  for edge: take from MeshRaw%v_nd_physindex
      !   [1, :]  for P3 vertex: take from MeshRaw%v_nd_physindex
      !   [1, :]  for P3 edges: take from P2 edge in MeshRaw%v_nd_physindex
      !   [1, :]  for P3 interior point:  0

      !   [2, :]  dimension: face(2), P2 edge (1), P3 edges and interior (0)
      integer(8), dimension(:,:), allocatable :: v_ety_props

      ! Maps P2 edge nodes to surrounding vertices
      integer(8) edge_ends(2,3)

   contains

      procedure :: allocate => MeshEntities_allocate
      procedure :: build_mesh_tables => MeshEntities_build_mesh_tables
      procedure :: count_and_label_faces => MeshEntities_count_and_label_faces
      procedure :: count_and_label_edges => MeshEntities_count_and_label_edges
      procedure :: count_and_label_nodes_P3 => MeshEntities_count_and_label_nodes_P3
      procedure :: analyse_face_and_edges => MeshEntities_analyse_face_and_edges
      procedure :: analyse_p3_nodes => MeshEntities_analyse_p3_nodes


   end type MeshEntities





   ! A class for holding the various mesh and node properties together.
   ! Also for exchanging mesh info between fortran and python
   type, public  :: MeshRawAC

      !character(FNAME_LENGTH) :: mesh_file

      integer(8) :: n_msh_pts
      integer(8) :: n_msh_el
      integer(8) :: n_elt_mats

      ! what kind of mesh point [n_msh_pts]
      ! basically the Gmsh physindex.  Nonzero = boundary node
      integer(8), dimension(:), allocatable :: v_nd_physindex

      ! absolute coordinates of each node in the mesh  [2, n_msh_pts]
      double precision, dimension(:,:), allocatable :: v_nd_xy

      ! material of each elt [n_msh_el]
      integer(8), dimension(:), allocatable :: el_material

      ! mapping of elt,node to index of mesh point, [6 x n_msh_el]
      ! .mail elt data, but transposed: tn[i,j] = ith node of elt j
      !  tn[1..3,j] = vertex nodes, [4..6, j] edge nodes

      integer(8), dimension(:,:), allocatable :: elnd_to_mshpt

   contains

      procedure :: allocate => MeshRawAC_allocate
      !final :: MeshRaw_destructor

      !procedure :: find_nodes_for_elt => MeshRaw_find_nodes_for_elt
      !procedure :: fill_python_arrays => MeshRaw_fill_python_arrays
      !procedure :: is_boundary_node => MeshRaw_is_boundary_node
      !procedure :: is_boundary_node_2 => MeshRaw_is_boundary_node_2
!
!      procedure :: node_phys_index_by_ref => MeshRaw_node_phys_index_by_ref
!
      procedure :: construct_node_tables => MeshRawAC_construct_node_tables
      procedure :: load_node_tables_from_py => MeshRawAC_load_node_tables_from_py



   end type MeshRawAC


   type, public  :: MeshEntitiesAC

      integer(8) n_faces       ! num faces
      integer(8) n_edges       ! num physically distinct edges
      integer(8) n_msh_pts_p3  ! num physically distinct P3 points
      integer(8) n_entities         ! total num of physically distinct things (sum of last 3)

      ! v_tags[i, j] 14 x n_msh_el
      !  tn[1, j] = j  ! index of element
      !    -
      !  Each element has 1 face, 3 edges and 10 P3 nodes
      !  so v_tags has dimensions 14 x n_msh_el
      integer(8), dimension(:,:), allocatable :: v_tags


      ! Positions indexed by tag
      ! [2 x n_ddl]
      ! [:, 1] - barycentre of the face
      ! [:, 2..4] - P2 edge positions cloned from MeshRaw%v_nd_xy
      ! [:, 5..7] - P3 vertex positions cloned from MeshRaw%v_nd_xy
      ! [:, 8..13] - P3 edge positions built from MeshRaw%v_nd_xy at 1/3, 2/3 intervals
      ! [:, 14] - P3 interior at barycentre
      double precision, dimension(:,:), allocatable :: v_xy


      !   [2 x n_ddl]
      !   [1, :]  for face: interior (0)
      !   [1, :]  for edge: take from MeshRaw%node_phys_i
      !   [1, :]  for P3 vertex: take from MeshRaw%node_phys_i
      !   [1, :]  for P3 edges: take from P2 edge in MeshRaw%node_phys_i
      !   [1, :]  for P3 interior point:  0

      !   [2, :]  dimension: face(2), P2 edge (1), P3 edges and interior (0)
      integer(8), dimension(:,:), allocatable :: v_ety_props

      ! Maps P2 edge nodes to surrounding vertices
      integer(8) edge_ends(2,3)

   contains

      procedure :: allocate => MeshEntitiesAC_allocate
      ! procedure :: build_mesh_tables => MeshEntitiesAC_build_mesh_tables
      ! procedure :: count_and_label_faces => MeshEntitiesAC_count_and_label_faces
      ! procedure :: count_and_label_edges => MeshEntitiesAC_count_and_label_edges
      ! procedure :: count_and_label_nodes_P3 => MeshEntitiesAC_count_and_label_nodes_P3
      ! procedure :: analyse_face_and_edges => MeshEntitiesAC_analyse_face_and_edges
      ! procedure :: analyse_p3_nodes => MeshEntitiesAC_analyse_p3_nodes


   end type MeshEntitiesAC

contains

#include "meshprops_impl.f90"


#include "meshentities_impl.f90"

#include "meshprops_impl_ac.f90"


end module class_MeshRaw
