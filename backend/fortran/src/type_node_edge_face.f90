
subroutine type_node_edge_face (n_msh_el, n_msh_pts, nodes_per_el, n_ddl, &
   type_nod, table_nod, table_N_E_F, &
   visited, type_N_E_F, mesh_xy, x_E_F)

!!!!!!!!!!!!!!!!

   implicit none
   integer(8) n_msh_el, n_msh_pts, nodes_per_el, n_ddl
   integer(8) type_nod(n_msh_pts)
   integer(8) table_nod(nodes_per_el,n_msh_el), table_N_E_F(14,n_msh_el)
   integer(8) visited(n_ddl), type_N_E_F(2,n_ddl)
   double precision mesh_xy(2,n_msh_pts), x_E_F(2,n_ddl)
   double precision, parameter :: one_third = 1.d0/3.d0


   integer(8) i, j, j1
   integer(8) type_n(10)
!  integer(8) list_point_F(6,4)
!  integer(8), parameter :: nddl_0 =14

   double precision el_xy(2,6)


   if ( nodes_per_el .ne. 6 ) then
      write(*,*) "type_node_edge_face: problem nodes_per_el = ", nodes_per_el
      write(*,*) "type_node_edge_face: nodes_per_el should be equal to 6 !"
      write(*,*) "type_node_edge_face: Aborting..."
      stop
   endif


   !  Initialisation
   type_N_E_F = 0
   visited= 0

   !  do j=1,n_ddl   !OBSELETE
   !  type_N_E_F(1,j) = 0
   !  type_N_E_F(2,j) = 0
   !  enddo

   !  do i=1,n_ddl   !OBSELETE
   !  visited(i) = 0
   !  enddo


   do i=1,n_msh_el

      do j=1,nodes_per_el
         j1 = table_nod(j,i)
         type_n(j) = type_nod(j1)
         el_xy(1,j) = mesh_xy(1,j1)
         el_xy(2,j) = mesh_xy(2,j1)
      enddo

      !  an element is a face
      j=1
      j1 = table_N_E_F(j,i)

      !  centre of the elements
      !x_E_F(1,j1) = (el_xy(1,1) + el_xy(1,2) + el_xy(1,3))/3.0d0
      !x_E_F(2,j1) = (el_xy(2,1) + el_xy(2,2) + el_xy(2,3))/3.0d0

      x_E_F(:,j1) = (el_xy(:,1) + el_xy(:,2) + el_xy(:,3)) * one_third

      !  Topologically, a face is an interior domain
      type_N_E_F(1,j1) = 0

      !  Face => dimension two
      type_N_E_F(2,j1) = 2

      !  scan the 3 element edges
      do j=1,3
         j1 = table_N_E_F(j+1,i)
         !x_E_F(1,j1) = el_xy(1,j+3)
         !x_E_F(2,j1) = el_xy(2,j+3)

         x_E_F(:,j1) = el_xy(:,j+3)

         if (visited(j1) .eq. 0) then
            visited(j1) = 1
            type_N_E_F(1,j1) = type_n(j+3)

            !  Edge => dimension one
            type_N_E_F(2,j1) = 1
         endif

      enddo

   enddo

   return
end
