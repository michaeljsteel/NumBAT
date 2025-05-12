#include "numbat_decl.h"

 ! !  Storage locations in sequence
 ! !  - tab_N_E_F = table_N_E_F,   shape: 14 x n_msh_elts
 ! !  - table_edges     shape: 4 x n_msh_pts
 ! !

 ! !  V = number of vertices
 ! !  E = number of edges
 ! !  F = number of faces
 ! !  C = number of cells (3D, triangle)
 ! !
 ! !  From Euler's theorem on 3D graphs: V-E+F-C = 1 - (number of holes)
 ! !  n_msh_pts = (number of vertices) + (number of mid-edge point) = V + E;
 ! !

subroutine build_mesh_tables( n_msh_elts, n_msh_pts, nodes_per_el, n_ddl, &
   mesh_raw, entities, errco, emsg)

   use numbatmod
   use alloc
   use class_Mesh

   integer(8) n_msh_elts, n_msh_pts, nodes_per_el, n_ddl

   type(MeshEM) :: mesh_raw
   type(MeshEntities) :: entities

   integer(8), intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) :: emsg

   !  ----------------------------------------------
   integer(8) n_face, n_edge, n_msh_pts_p3
   integer(8) ui_out


   integer(8), dimension(:), allocatable :: visited

   !  ----------------------------------------------

   ui_out = stdout

   call integer_alloc_1d(visited, mesh_raw%n_ddl, 'visited', errco, emsg); RETONERROR(errco)


   !  Fills:  tab_N_E_F[1,:]
   call label_faces (mesh_raw%n_msh_elts, entities%m_elnd_to_mshpt)


   !  For P2 FEM n_msh_pts=N_Vertices+N_Edge
   !  note: each element has 1 face, 3 edges and 10 P3 nodes
   !  so table_N_E_F = tab_N_E_F has dimensions 14 x n_msh_elts

   !  each element is a face
   n_face = mesh_raw%n_msh_elts

   !  Fills: n_edge, table_edge[1..4,:], tab_N_E_F[2:4,:], visited[1:n_msh_pts]
   !  Todo!  move n_edge later in list as an out variable
   call label_edges (mesh_raw, entities, n_edge, visited)

   !  Fills: remainder of tab_N_E_F[5:,:], visited[1:n_msh_pts], n_msh_pts_3
   !  Todo: move n_msh_pts_p3 later
   call label_nodes_P3 (n_msh_elts, n_msh_pts, nodes_per_el, n_edge, n_msh_pts_p3, mesh_raw%m_elnd_to_mshpt, &
      entities%m_elnd_to_mshpt,  visited)


   !  TODO: what is signif of this quanitty?
   n_ddl = n_edge + n_face + n_msh_pts_p3

   write(*,*), 'dddls', n_ddl, n_edge , n_face , n_msh_pts_p3


   ! if (debug .eq. 1) then
   !    write(ui_out,*) "py_calc_modes.f: n_msh_pts, n_msh_elts = ", n_msh_pts, n_msh_elts
   !    write(ui_out,*) "py_calc_modes.f: n_msh_pts_p3 = ", n_msh_pts_p3
   !    write(ui_out,*) "py_calc_modes.f: n_vertex, n_edge, n_face,", " n_msh_elts = ", &
   !       (n_msh_pts - n_edge), n_edge, n_face, n_msh_elts
   !    write(ui_out,*) "py_calc_modes.f: 2D case of the Euler &
   !    & characteristic: V-E+F=1-(number of holes)"
   !    write(ui_out,*) "py_calc_modes.f: Euler characteristic: V - E + F &
   !    &= ", (n_msh_pts - n_edge) - n_edge + n_face
   ! endif


   !  Fills: entities%type_nod(1:2, 1:n_ddl), x_E_F(1:2, 1:n_ddl)
   !  Should be using c_dwork for x_E_F ?
   call type_node_edge_face (n_msh_elts, n_msh_pts, nodes_per_el, n_ddl, mesh_raw%type_nod, mesh_raw%m_elnd_to_mshpt, &
      entities%m_elnd_to_mshpt, visited , entities%type_nod, mesh_raw%v_mshpt_xy, entities%v_mshpt_xy )


   !  Fills: entities%type_nod(1:2, 1:n_ddl), x_E_F(1:2, 1:n_ddl)
   call get_coord_p3 (n_msh_elts, n_msh_pts, nodes_per_el, n_ddl, mesh_raw%m_elnd_to_mshpt, mesh_raw%type_nod, &
      entities%m_elnd_to_mshpt, entities%type_nod, mesh_raw%v_mshpt_xy, entities%v_mshpt_xy , visited)


end subroutine
