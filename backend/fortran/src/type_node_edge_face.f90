
subroutine type_node_edge_face (n_msh_el, n_msh_pts, nodes_per_el, n_ddl, &
   type_nod, table_nod, table_N_E_F, &
   visite, type_N_E_F, mesh_xy, x_E_F)
!
!***********************************************************************
!
   implicit none
   integer*8 n_msh_el, n_msh_pts, nodes_per_el, n_ddl
   integer*8 type_nod(n_msh_pts)
   integer*8 table_nod(nodes_per_el,n_msh_el), table_N_E_F(14,n_msh_el)
   integer*8 visite(n_ddl), type_N_E_F(2,n_ddl)
   double precision mesh_xy(2,n_msh_pts), x_E_F(2,n_ddl)


   integer*8 i, j, j1
   integer*8 type_n(10)
!      integer*8 list_point_F(6,4)
   integer*8 nddl_0
   parameter (nddl_0 = 14)
   double precision xel(2,6)


   if ( nodes_per_el .ne. 6 ) then
      write(*,*) "type_node_edge_face: problem nodes_per_el = ", nodes_per_el
      write(*,*) "type_node_edge_face: nodes_per_el should be equal to 6 !"
      write(*,*) "type_node_edge_face: Aborting..."
      stop
   endif


   !     Initialisation
   do j=1,n_ddl
      type_N_E_F(1,j) = 0
      type_N_E_F(2,j) = 0
   enddo

   do i=1,n_ddl
      visite(i) = 0
   enddo
!
   do i=1,n_msh_el

      do j=1,nodes_per_el
         j1 = table_nod(j,i)
         type_n(j) = type_nod(j1)
         xel(1,j) = mesh_xy(1,j1)
         xel(2,j) = mesh_xy(2,j1)
      enddo

      ! an element is a face
      j=1
      j1 = table_N_E_F(j,i)

      ! centre of the elements
      x_E_F(1,j1) = (xel(1,1) + xel(1,2) + xel(1,3))/3.0d0
      x_E_F(2,j1) = (xel(2,1) + xel(2,2) + xel(2,3))/3.0d0

      ! Topologically, a face is an interior domain
      type_N_E_F(1,j1) = 0

      ! Face => dimension two
      type_N_E_F(2,j1) = 2

      ! scan the 3 element edges
      do j=1,3
         j1 = table_N_E_F(j+1,i)
         x_E_F(1,j1) = xel(1,j+3)
         x_E_F(2,j1) = xel(2,j+3)

         if (visite(j1) .eq. 0) then
            visite(j1) = 1
            type_N_E_F(1,j1) = type_n(j+3)

            ! Edge => dimension one
            type_N_E_F(2,j1) = 1
         endif

      enddo

   enddo
!
   return
end
