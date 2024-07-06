#include "numbat_decl.h"

subroutine build_mesh_tables( &
   n_msh_el, n_msh_pts, d_nodes_per_el, n_ddl, &
   type_nod, table_nod, mesh_xy, &
   type_N_E_F, table_N_E_F, xy_N_E_F, &
   debug, errco, emsg)

   use numbatmod
use alloc
   integer(8) n_msh_el, n_msh_pts, d_nodes_per_el, n_ddl
   integer(8), intent(out) :: type_nod(n_msh_pts)
   integer(8), intent(out) :: table_nod(d_nodes_per_el, n_msh_el)
   !integer(8), intent(out) :: type_el(n_msh_el)

   double precision, intent(out) :: mesh_xy(2,n_msh_pts)
   integer(8)  :: type_N_E_F(2, n_ddl)
   integer(8) :: table_N_E_F(14, n_msh_el)
   double precision :: xy_N_E_F(2, n_ddl)

   integer(8) debug
   integer(4), intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) :: emsg

   !  ----------------------------------------------
   integer(8) n_face, n_edge, n_msh_pts_p3
   integer(8) ui_out


   integer(8), dimension(:), allocatable :: visited

   !  ----------------------------------------------


   ui_out = stdout
   call integer_alloc_1d(visited, n_ddl, 'visited', errco, emsg); RETONERROR(errco)


   !  Fills:  table_edge_face[1,:]
   call list_face (n_msh_el, table_N_E_F)


   !  For P2 FEM n_msh_pts=N_Vertices+N_Edge
   !  note: each element has 1 face, 3 edges and 10 P3 nodes
   !  so table_N_E_F = table_edge_face has dimensions 14 x n_msh_el

   !  each element is a face
   n_face = n_msh_el

   !  Fills: n_edge, table_edge[1..4,:], table_edge_face[2:4,:], visited[1:n_msh_pts]
   !  Todo!  move n_edge later in list as an out variable
   call list_edge (n_msh_el, n_msh_pts, d_nodes_per_el, n_edge, type_nod, table_nod, &
      table_N_E_F, visited)

   !  Fills: remainder of table_edge_face[5:,:], visited[1:n_msh_pts], n_msh_pts_3
   !  Todo: move n_msh_pts_p3 later
   call list_node_P3 (n_msh_el, n_msh_pts, d_nodes_per_el, n_edge, n_msh_pts_p3, table_nod, &
      table_N_E_F,  visited)


   !  TODO: what is signif of this quanitty?
   n_ddl = n_edge + n_face + n_msh_pts_p3


   if (debug .eq. 1) then
      write(ui_out,*) "py_calc_modes.f: n_msh_pts, n_msh_el = ", n_msh_pts, n_msh_el
      write(ui_out,*) "py_calc_modes.f: n_msh_pts_p3 = ", n_msh_pts_p3
      write(ui_out,*) "py_calc_modes.f: n_vertex, n_edge, n_face,", " n_msh_el = ", &
         (n_msh_pts - n_edge), n_edge, n_face, n_msh_el
      write(ui_out,*) "py_calc_modes.f: 2D case of the Euler &
      & characteristic: V-E+F=1-(number of holes)"
      write(ui_out,*) "py_calc_modes.f: Euler characteristic: V - E + F &
      &= ", (n_msh_pts - n_edge) - n_edge + n_face
   endif


   !  Fills: type_N_E_F(1:2, 1:n_ddl), x_E_F(1:2, 1:n_ddl)
   !  Should be using c_dwork for x_E_F ?
   call type_node_edge_face (n_msh_el, n_msh_pts, nodes_per_el, n_ddl, type_nod, table_nod, &
      table_N_E_F, visited , type_N_E_F, mesh_xy, xy_N_E_F )


   !  Fills: type_N_E_F(1:2, 1:n_ddl), x_E_F(1:2, 1:n_ddl)
   call get_coord_p3 (n_msh_el, n_msh_pts, nodes_per_el, n_ddl, table_nod, type_nod, &
      table_N_E_F, type_N_E_F, mesh_xy, xy_N_E_F, visited)
end subroutine
